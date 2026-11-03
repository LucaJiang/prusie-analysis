"""Strict population QC and PLINK conversions from the original PGEN."""
from __future__ import annotations
from pathlib import Path
from collections import Counter, defaultdict
import csv, gzip, json, os, shutil, time
import numpy as np
import pandas as pd
import pgenlib
import pysam
from .acquire import atomic_json, sha256
from .reference import _table, genotype_qc
from .strict_policy import POLICY, REFERENCE_STAGES, classify, ALLELE_SCOPE
from .strict_common import allele_detail, file_record, stages
from .strict_prepare import code_identity
from .runtime import paths, log_command, fasta_identity
class BedReader:
    """Independent SNP-major BED decoder, returning diploid BIM A1 counts."""
    def __init__(self, path, n, m):
        self.n, self.m = n, m
        self.stride = (n+3)//4
        self.handle = Path(path).open("rb")
        if self.handle.read(3) != b"\x6c\x1b\x01":
            raise ValueError("Invalid SNP-major BED header")
        if Path(path).stat().st_size != 3 + self.stride*m:
            raise ValueError("BED byte dimensions mismatch")

    def read(self, indices):
        indices = np.asarray(indices, dtype=np.int64)
        if np.any(indices < 0) or np.any(indices >= self.m):
            raise ValueError("BED variant index out of bounds")
        packed = np.empty((len(indices), self.stride), dtype=np.uint8)
        for j, index in enumerate(indices):
            self.handle.seek(3 + int(index)*self.stride)
            packed[j] = np.frombuffer(self.handle.read(self.stride), dtype=np.uint8)
        bits = ((packed[:,:,None] >> np.array([0,2,4,6], dtype=np.uint8)) & 3).reshape(len(indices), -1)[:, :self.n]
        return np.array([2,-9,1,0], dtype=np.int8)[bits]

    def close(self):
        self.handle.close()

class StrictReferenceReader:
    def __init__(self, prefix):
        self.prefix = Path(prefix)
        self.receipt = json.loads(Path(str(prefix)+".receipt.json").read_text())
        if self.receipt.get("state") != "READY" or self.receipt.get("policy") != POLICY:
            raise ValueError("Strict reference is not READY")
        self.bed = BedReader(str(prefix)+".bed", self.receipt["n_reference"], self.receipt["retained_variants"])
        self.variants = pysam.TabixFile(str(prefix)+".variants.tsv.gz")
        self.columns = self.variants.header[-1].lstrip("#").split("\t")
        self.qc = pysam.TabixFile(str(prefix)+".qc.tsv.gz")
        self.qc_columns = self.qc.header[-1].lstrip("#").split("\t")

    def fetch(self, chrom, start, end):
        return [dict(zip(self.columns, line.split("\t"), strict=True))
                for line in self.variants.fetch(str(chrom), max(0,start-1), end)]

    def read_alt_dosages(self, indices):
        a = self.bed.read(indices).T.astype(np.float64)
        a[a == -9] = np.nan
        return a

    def fetch_qc(self, chrom, start, end):
        return [dict(zip(self.qc_columns,line.split("\t"),strict=True))
                for line in self.qc.fetch(str(chrom),max(0,start-1),end)]

    def close(self):
        self.bed.close()
        self.variants.close()
        self.qc.close()

def idset(path):
    return set(Path(path).read_text().splitlines())

def prepare(run: Path, chrom: int):
    p = paths(run)
    root = p["reference"]
    work = root / "work" / f"chr{chrom}"
    work.mkdir(parents=True, exist_ok=True)
    sourceval = json.loads((p["index"] / "source_validation.json").read_text())
    for v in sourceval["source"].values():
        st = Path(v["path"]).stat()
        if (st.st_size, st.st_mtime_ns) != (v["size_bytes"], v["mtime_ns"]):
            raise ValueError("Original source changed since full validation")
    sample = {pop:_table(p["samples"] / f"{pop}.samples.tsv") for pop in ("FIN","EUR")}
    identity = dict(policy=POLICY, build="GRCh38", source=sourceval["source"],fasta=fasta_identity(run,"GRCh38"),
                    allele_scope=ALLELE_SCOPE,
                    sample_sha256={pop:sha256(p["samples"] / f"{pop}.samples.tsv") for pop in sample},
                    plink_sha256={k:sha256(p[k]) for k in ("plink1","plink2")},
                    code={name:d for name,d in code_identity().items() if name in ("strict_reference.py","strict_common.py","strict_policy.py","reference.py")})
    receipts = {pop:root / pop / f"chr{chrom}.receipt.json" for pop in sample}
    with pysam.FastaFile(str(p["fasta"])) as fa:
        sequence = fa.fetch(f"chr{chrom}").upper()
    sourceledger = work / "structural_qc.tsv.gz"
    candidates, basecounts = [], Counter()
    raw_details=Counter()
    with gzip.open(p["index"] / f"chr{chrom}.tsv.gz", "rt") as inp, gzip.open(sourceledger, "wt") as out:
        next(inp)
        out.write("chrom\tpos\toriginal_id\tref\talt\tsource_index\tfirst_failure\n")
        for line in inp:
            c, ps, vid, ref, alt, si, _ = line.rstrip().split("\t")
            pos = int(ps)
            raw_details[allele_detail(ref,alt)]+=1
            reason = "invalid_coordinate" if c != str(chrom) or not 1 <= pos <= len(sequence) else classify(ref, alt)
            if reason == "PASS" and sequence[pos-1] != ref:
                reason = "REF_mismatch"
            basecounts[reason] += 1
            out.write(f"{c}\t{pos}\t{vid}\t{ref}\t{alt}\t{si}\t{reason}\n")
            if reason in ("PASS", "palindrome"):
                candidates.append((pos, vid, ref, alt, int(si)))
    del sequence
    print(json.dumps(dict(chrom=chrom, structural_counts=basecounts)), flush=True)
    # Canonical duplicates are compared within each final population, before
    # genotype QC. Different alternate alleles at one position stay distinct.
    groups = defaultdict(list)
    previous_pos = None
    local = {}
    for i, row in enumerate(candidates):
        if row[0] != previous_pos:
            local = {}; previous_pos = row[0]
        key = (row[0], row[2], row[3])
        if key in local:
            groups[key] = [local[key]] if key not in groups else groups[key]
            groups[key].append(i)
        else:
            local[key] = i
    euridx = np.array([int(r["source_sample_index"]) for r in sample["EUR"]],dtype=np.uint32)
    finidx = np.array([int(r["source_sample_index"]) for r in sample["FIN"]],dtype=np.uint32)
    fincol = np.searchsorted(euridx,finidx)
    assert np.array_equal(euridx[fincol],finidx)
    source = pgenlib.PgenReader(os.fsencode(str(p["source"])+".pgen"), sample_subset=euridx,
                              allele_idx_offsets=np.load(p["index"] / "source_allele_offsets.npy", mmap_mode="r"))
    dup = {pop:{} for pop in sample}
    for key, indices in groups.items():
        values = np.empty((len(indices),len(euridx)),dtype=np.int8)
        source.read_list(np.array([candidates[i][4] for i in indices],dtype=np.uint32),values,allele_idx=1)
        for pop in sample:
            g = values if pop == "EUR" else values[:,fincol]
            identical = np.all(g == g[0])
            for i in indices[1:] if identical else indices:
                dup[pop][candidates[i][4]] = "identical_duplicate_collapsed" if identical else "conflicting_duplicate"
    allresults=[]
    for pop in ("FIN","EUR"):
        dest = root / pop
        dest.mkdir(exist_ok=True)
        prefix = dest / f"chr{chrom}"
        temp = work / pop
        temp.mkdir(exist_ok=True)
        selected = [r for r in candidates if r[4] not in dup[pop]]
        n, m = len(sample[pop]),len(selected)
        cp = temp / "candidate"
        palids = [f"s{r[4]}" for r in selected if classify(r[2],r[3]) == "palindrome"]
        (temp / "palindrome.ids").write_text("\n".join(palids)+"\n")
        with pgenlib.PgenWriter(os.fsencode(str(cp)+".pgen.part"),n,variant_ct=m,nonref_flags=False) as writer, \
                Path(str(cp)+".pvar.part").open("w") as pv:
            pv.write("##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\n")
            for start in range(0,m,8192):
                batch=selected[start:start+8192]
                values=np.empty((len(batch),len(euridx)),dtype=np.int8)
                source.read_list(np.array([r[4] for r in batch],dtype=np.uint32),values,allele_idx=1)
                writer.append_biallelic_batch(values if pop=="EUR" else np.ascontiguousarray(values[:,fincol]))
                pv.writelines(f"{chrom}\t{r[0]}\ts{r[4]}\t{r[2]}\t{r[3]}\n" for r in batch)
        for ext in ("pgen","pvar"):
            Path(str(cp)+"."+ext+".part").replace(Path(str(cp)+"."+ext))
        shutil.copyfile(p["samples"] / f"{pop}.psam",Path(str(cp)+".psam"))
        base=[p["plink2"],"--pfile",cp,"--threads","1","--memory","4096","--nonfounders","--no-pheno",
              "--exclude",temp/"palindrome.ids","--min-alleles","2","--max-alleles","2","--snps-only","just-acgt"]
        log_command(run,[*base,"--maf","0.05","--write-snplist","--out",temp/"maf"],temp/"maf.stdout.log")
        maf=idset(temp/"maf.snplist")
        log_command(run,[*base,"--extract",temp/"maf.snplist","--geno","0.02","--write-snplist","--out",temp/"geno"],temp/"geno.stdout.log")
        geno=idset(temp/"geno.snplist")
        log_command(run,[*base,"--extract",temp/"geno.snplist","--hardy","midp","--freq","counts",
                         "--missing","variant-only","--out",temp/"metrics"],temp/"metrics.stdout.log")
        log_command(run,[*base,"--extract",temp/"geno.snplist","--hwe","1e-6","0","midp",
                         "--make-bed","--out",temp/"plink2"],temp/"plink2.stdout.log")
        log_command(run,[p["plink1"],"--bfile",temp/"plink2","--keep-allele-order","--make-bed",
                         "--threads","1","--memory","4096","--out",str(prefix)+".part"],temp/"plink1.stdout.log")
        bim1=[l.split() for l in Path(str(prefix)+".part.bim").read_text().splitlines()]
        bim2=[l.split() for l in (temp/"plink2.bim").read_text().splitlines()]
        if bim1 != bim2:
            raise ValueError("PLINK1 changed allele or variant order")
        fam1=[l.split() for l in Path(str(prefix)+".part.fam").read_text().splitlines()]
        fam2=[l.split() for l in (temp/"plink2.fam").read_text().splitlines()]
        if fam1 != fam2 or [r[1] for r in fam1] != [r["sample_id"] for r in sample[pop]]:
            raise ValueError("PLINK sample order mismatch")
        by_source={r[4]:r for r in selected}
        final=[]
        for b in bim1:
            row=by_source[int(b[1][1:])]
            if (b[0],int(b[3]),b[4],b[5]) != (str(chrom),row[0],row[3],row[2]):
                raise ValueError("PVAR REF/ALT disagrees with PLINK2/1 A2/A1")
            final.append(row)
        hwe=pd.read_csv(temp/"metrics.hardy",sep="\t").set_index("ID")
        hcol="MIDP" if "MIDP" in hwe else "P" if "P" in hwe else "P_HWE"
        if hcol not in hwe:
            raise ValueError(f"Unexpected HWE report fields: {list(hwe)}")
        bed1=BedReader(str(prefix)+".part.bed",n,len(final))
        bed2=BedReader(temp/"plink2.bed",n,len(final))
        # pgenlib independently reads final BED as well; the handwritten
        # decoder's A1 interpretation is checked against ORIGINAL PGEN ALT.
        secondary=pgenlib.PgenReader(os.fsencode(str(prefix)+".part.bed"),raw_sample_ct=n)
        stats=[]
        for start in range(0,len(final),4096):
            batch=final[start:start+4096]
            original=np.empty((len(batch),len(euridx)),dtype=np.int8)
            source.read_list(np.array([r[4] for r in batch],dtype=np.uint32),original,allele_idx=1)
            expected=original if pop=="EUR" else original[:,fincol]
            indices=np.arange(start,start+len(batch),dtype=np.uint32)
            a,b=bed1.read(indices),bed2.read(indices)
            c=np.empty_like(a);secondary.read_list(indices,c,allele_idx=1)
            if not (np.array_equal(a,expected) and np.array_equal(b,expected) and np.array_equal(c,expected)):
                raise ValueError("Full source/PLINK2/PLINK1/decoder genotype concordance failed")
            af,mf,missing,_=genotype_qc(a,maf_min=.05,missing_max=.02)
            for k,row in enumerate(batch):
                pv=float(hwe.loc[f"s{row[4]}",hcol])
                if mf[k] < .05-1e-15 or missing[k]>.02+1e-15 or pv<1e-6*(1-1e-5):
                    raise ValueError("Final reference QC metric disagrees with PLINK decision")
                observed=a[k][a[k]>=0]
                if np.ptp(observed)==0:
                    raise ValueError("Unexpected zero variance after HWE QC; preserve and investigate")
                stats.append((float(af[k]),float(mf[k]),float(missing[k]),pv))
        bed1.close();bed2.close();secondary.close()
        passed={f"s{r[4]}" for r in final}
        if not passed<=geno<=maf:
            raise ValueError("Staged accepted sets are not nested")
        exclusions=dict(dup[pop])
        for row in selected:
            vid=f"s{row[4]}"
            exclusions[row[4]]=("reference_MAF_lt_0.05" if vid not in maf else
                "reference_missing_gt_0.02" if vid not in geno else
                "reference_HWE_midp_lt_1e-6" if vid not in passed else "PASS")
        counts=Counter()
        with gzip.open(sourceledger,"rt") as inp, pysam.BGZFile(str(prefix)+".qc.tsv.gz.part","w") as out:
            out.write(("#"+next(inp)).encode())
            for line in inp:
                row=line.rstrip().split("\t")
                if row[-1]=="PASS":row[-1]=exclusions[int(row[5])]
                counts[row[-1]]+=1
                out.write(("\t".join(row)+"\n").encode())
        columns=["chrom","pos","id","internal_id","ref","alt","source_index","variant_index","af","maf","missing_rate",
                 "hwe_midp","original_pos","original_ref","original_alt","plink2_A1","plink2_A2","final_A1","final_A2","counted_allele"]
        with pysam.BGZFile(str(prefix)+".variants.tsv.gz.part","w") as out:
            out.write(("#"+"\t".join(columns)+"\n").encode())
            for i,(row,st) in enumerate(zip(final,stats,strict=True)):
                pos,vid,ref,alt,si=row
                vals=[chrom,pos,vid,f"s{si}",ref,alt,si,i,*st,pos,ref,alt,alt,ref,alt,ref,alt]
                out.write(("\t".join(map(str,vals))+"\n").encode())
        for ext in ("bed","bim","fam"):
            Path(str(prefix)+".part."+ext).replace(Path(str(prefix)+"."+ext))
        for ext in ("variants.tsv.gz","qc.tsv.gz"):
            Path(str(prefix)+"."+ext+".part").replace(Path(str(prefix)+"."+ext))
        pysam.tabix_index(str(prefix)+".variants.tsv.gz",seq_col=0,start_col=1,end_col=1,force=True)
        pysam.tabix_index(str(prefix)+".qc.tsv.gz",seq_col=0,start_col=1,end_col=1,force=True)
        hwe_removed=hwe.loc[list(geno-passed)] if geno-passed else hwe.iloc[:0]
        strata={}
        # HWE impacts stratified by observed frequency in the same genotype
        # set; PLINK report genotype counts provide an independent denominator.
        for lo,hi in ((.05,.1),(.1,.2),(.2,.3),(.3,.4),(.4,.500001)):
            if all(c in hwe for c in ("HOM_A1_CT","HET_A1_CT","TWO_AX_CT")):
                ac=2*hwe.HOM_A1_CT+hwe.HET_A1_CT
                obs=hwe.HOM_A1_CT+hwe.HET_A1_CT+hwe.TWO_AX_CT
            elif all(c in hwe for c in ("HOM_REF_CT","HET_CT","HOM_ALT1_CT")):
                ac=2*hwe.HOM_ALT1_CT+hwe.HET_CT
                obs=hwe.HOM_REF_CT+hwe.HET_CT+hwe.HOM_ALT1_CT
            else:
                break
            f=np.minimum(ac/(2*obs),1-ac/(2*obs))
            mask=(f>=lo)&(f<hi)
            strata[f"{lo}-{hi}"]={"tested":int(mask.sum()),"removed":int((mask & ~hwe.index.isin(passed)).sum())}
        reports_dir=dest / f"chr{chrom}.plink"
        reports_dir.mkdir(exist_ok=True)
        for suffix in ('log','nosex'):
            artifact=Path(str(prefix)+'.part.'+suffix)
            if artifact.exists():artifact.replace(reports_dir/('plink1.final.'+suffix))
        for path in temp.iterdir():
            if path.suffix in (".log",".snplist"):
                shutil.copyfile(path,reports_dir/path.name)
            elif path.name.startswith("metrics.") and path.suffix not in (".log",):
                with path.open("rb") as inp,gzip.open(reports_dir/(path.name+".gz"),"wb") as out:shutil.copyfileobj(inp,out)
        result=dict(state="READY",policy=POLICY,build="GRCh38",release="IGSR NYGC 20220422_3202_phased_SNV_INDEL_SV",
                    population=pop,chromosome=chrom,n_reference=n,source_variants=sum(counts.values()),retained_variants=len(final),
                    first_failure_counts=dict(counts),raw_allele_representation_counts=dict(raw_details),stages=stages(counts,REFERENCE_STAGES),identity=identity,
                    allele_scope=ALLELE_SCOPE,
                    historical_site_cardinality='not asserted; no historical ancestry lookup or auxiliary VCF dependency',
                    hwe=dict(threshold=1e-6,midp=True,k=0,keep_fewhet=False,nonfounders_included=True,
                             tested=len(geno),removed=len(geno-passed),frequency_strata=strata,
                             report_columns=list(hwe.columns),phenotype="none; all selected unrelated samples"),
                    conversion=dict(all_retained_genotypes_compared=True,genotype_cells_compared=n*len(final),
                                    source_ALT_equals_plink2_A1_equals_final_A1=True,
                                    independent_BED_decoder_and_pgenlib_both_match=True,sample_order_verified=True),
                    files=[file_record(Path(str(prefix)+"."+ext)) for ext in ("bed","bim","fam","variants.tsv.gz","variants.tsv.gz.tbi","qc.tsv.gz","qc.tsv.gz.tbi")]
                          +[file_record(item) for item in sorted(reports_dir.iterdir()) if item.is_file()],
                    completed_epoch=time.time())
        atomic_json(str(prefix)+".receipt.json",result)
        allresults.append(result)
        print(json.dumps(dict(chrom=chrom,population=pop,n=n,p=len(final),counts=counts,hwe=result["hwe"])),flush=True)
    source.close()
    return allresults
