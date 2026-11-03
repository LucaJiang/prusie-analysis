"""Portable SNP-only cases with exact signed float64 LD and retained strict QC."""
from __future__ import annotations
from pathlib import Path
from collections import Counter, defaultdict
import csv, gzip, hashlib, json, math, os, shutil, time
import numpy as np
from .acquire import atomic_json, sha256
from .prepare import common_parameters, factor_from_dosages, orient_variant
from .strict_policy import POLICY, STUDY_STAGES, classify
from .strict_common import orient, allele_detail, stages

def code_identity():
    folder=Path(__file__).parent
    return {p.name:sha256(p) for p in sorted(folder.glob("*.py"))}
def row_block(path, p, start, stop):
    return np.memmap(path,mode="r",dtype="<f8",offset=8*start*p,shape=(stop-start,p))

def dense_ld(f, path, block=256):
    """Write all rows of the unmodified Gram product; bounded address space."""
    p=f.shape[1]
    needed=8*p*p
    if shutil.disk_usage(path.parent).free < 2*needed + 20*1024**3:
        raise OSError(28,"Insufficient reserve for complete dense LD and column layout")
    with path.open("wb") as out:
        for i in range(0,p,block):
            a=np.asarray(f[:,i:i+block].T@f,dtype="<f8",order="C")
            a.tofile(out)
        out.flush();os.fsync(out.fileno())
    return needed

def validate_matrix(path, f, block=193):
    """Independent saved-file residual against a newly decoded genotype factor.

    An exact real Gram matrix is PSD. The saved matrix's residual plus the
    standard dot-product roundoff bound yields a conservative spectral bound
    p*(max_residual + gamma_n*max_column_norm_squared). No eigendecomposition,
    positive-definite repair, or rejection of legal low rank is performed.
    """
    n,p=f.shape
    if path.stat().st_size!=8*p*p:raise ValueError("LD size mismatch")
    diagonal=symmetry=maximum=residual=0.
    for i in range(0,p,block):
        stop=min(p,i+block)
        a=row_block(path,p,i,stop)
        if not np.isfinite(a).all():raise ValueError("Nonfinite saved LD")
        diagonal=max(diagonal,float(np.max(np.abs(a[np.arange(stop-i),np.arange(i,stop)]-1))))
        maximum=max(maximum,float(np.max(np.abs(a))))
        expected=f[:,i:stop].T@f
        residual=max(residual,float(np.max(np.abs(a-expected))))
        del expected
        for j in range(0,p,block):
            b=row_block(path,p,j,min(p,j+block))
            symmetry=max(symmetry,float(np.max(np.abs(a[:,j:j+block]-b[:,i:stop].T))))
            del b
        del a
    gamma=n*np.finfo(float).eps/(1-n*np.finfo(float).eps)
    norm2=float(np.max(np.sum(f*f,axis=0)))
    if diagonal>1e-10 or symmetry>1e-12 or maximum>1+1e-10 or residual>16*gamma*norm2:
        raise ValueError(f"Saved Gram validation failed: {diagonal}, {symmetry}, {maximum}, {residual}")
    return dict(diagonal_max_abs_error=diagonal,symmetry_max_abs_error=symmetry,max_abs_r=maximum,
                saved_Gram_residual_max_abs=residual,bitwise_symmetric=symmetry==0,
                psd_method="redecoded genotype Gram certificate, exhaustive saved-matrix residual; no full eigendecomposition",
                psd_eigenvalue_lower_bound=-p*(residual+gamma*norm2),dot_product_gamma_n=gamma,
                psd_bound_target="symmetric part of saved R; symmetry separately checked, diagnostic only with no changed bytes",
                bound_limitation="conservative float64 forward-error bound; legitimate rank <= N_ref-1 is retained",
                all_entries_reopened=True,n_variants=p,n_reference=n)

def duplicates(entries):
    """Keep one exact oriented record; reject every conflicting group member."""
    signature=lambda e:(e["beta"],e["se"],e["pvalue"],e.get("mlogp",""),e["study_maf"],e["aligned_study_af"])
    if len({signature(e) for e in entries}) != 1:
        return [],[(e,"conflicting_duplicate") for e in entries]
    return entries[:1],[(e,"identical_duplicate_collapsed") for e in entries[1:]]

def prepare_case(rows, *, case_id, chrom, start, end, n, reference, sequence, output, settings, provenance,
                 significant_threshold, lead_id=None, plan_only=False):
    output.mkdir(parents=True,exist_ok=True)
    identity=dict(policy=POLICY,provenance=provenance,settings=common_parameters(settings,n,0),n=n,
                  region=[chrom,start,end],code=code_identity(),
                  rows_sha256=hashlib.sha256(json.dumps(rows,sort_keys=True).encode()).hexdigest())
    signature=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
    receipt=output/"preparation.json"
    result=dict(case_id=case_id,chrom=str(chrom),start=start,end=end,n_study=n,n_ref=reference.receipt["n_reference"],
                source_rows=len(rows),raw_allele_representation_counts=Counter(allele_detail(r["other"],r["effect"]) for r in rows),identity=identity,identity_sha256=signature,state="PREPARING",started_epoch=time.time())
    atomic_json(receipt,result)
    normalized=defaultdict(list);events=[]
    for i,row in enumerate(rows):
        event=dict(source_row=row.get("source_row",i),original_id=str(row["original_id"]),
                   original_pos=row["original_pos"],original_effect_allele=row["effect"],original_other_allele=row["other"],
                   source_mapping_status=row.get('mapping_status','not_required_same_build'),
                   reason="",significant=False,source_lead=bool(row.get("source_lead",str(row["original_id"])==str(lead_id))))
        events.append(event)
        # Significance is annotated independently of the removal sequence.
        try:
            pv=float(row["pvalue"]);lp=float(row.get("mlogp",0))
            event["significant"]=(0<=pv<significant_threshold)
            if pv==0 and "mlogp" in row:
                event["significant"]=lp>-math.log10(significant_threshold)
        except (ValueError,TypeError):pass
        try:
            opos=int(row["original_pos"])
            if opos<1 or str(row.get("original_chrom",chrom)) not in set(map(str,range(1,23))):raise ValueError("invalid_coordinate")
            kind=classify(row["other"],row["effect"])
            if kind!="PASS":raise ValueError(kind)
            beta,se,pv,af=(float(row[k]) for k in ("beta","se","pvalue","af"))
            if not all(math.isfinite(v) for v in (beta,se,pv,af)) or se<=0 or not 0<=pv<=1 or not 0<=af<=1 or not math.isfinite(beta/se):
                raise ValueError("invalid_statistics")
            if "mlogp" in row and (not math.isfinite(float(row["mlogp"])) or float(row["mlogp"])<0):
                raise ValueError("invalid_statistics")
            if min(af,1-af)<.05:raise ValueError("study_MAF_lt_0.05")
            if row.get("source_fasta_valid") is False:raise ValueError("source_FASTA_mismatch")
            if row.get("mapping_status","PASS")!="PASS":raise ValueError("mapping_failure")
            pos=int(row["pos"])
            if not 1<=pos<=len(sequence):raise ValueError("invalid_coordinate")
            try:
                pos,ref,alt,sign,orientation=orient_variant(pos,row.get("target_other",row["other"]),
                                                         row.get("target_effect",row["effect"]),sequence)
            except ValueError as exc:raise ValueError("allele_mismatch") from exc
            net_sign,net_orientation=orient(row['other'],row['effect'],ref,alt)
            if net_sign!=sign:raise ValueError('allele_mismatch')
            event.update(pos=pos,ref=ref,alt=alt,flip_sign=sign,orientation=orientation)
            entry=dict(event=event,pos=pos,ref=ref,alt=alt,z=(beta/se)*sign,beta=beta*sign,se=se,pvalue=pv,
                       mlogp=row.get("mlogp",""),study_maf=min(af,1-af),study_af=af,aligned_study_af=af if sign==1 else 1-af,
                       original_beta=beta,original_se=se,original_pos=opos,original_effect_allele=row["effect"],
                       original_other_allele=row["other"],mapping_strand=row.get("mapping_strand","+"),
                       original_id=str(row["original_id"]),flip_sign=sign,orientation=orientation,
                       original_to_final_orientation=net_orientation)
            normalized[(pos,ref,alt)].append(entry)
        except (ValueError,TypeError,KeyError) as exc:
            reason=str(exc)
            if reason not in STUDY_STAGES:raise
            event["reason"]=reason
    refrows=reference.fetch(chrom,start,end)
    bykey=defaultdict(list)
    for r in refrows:bykey[(int(r["pos"]),r["ref"],r["alt"])].append(r)
    qc=reference.fetch_qc(chrom,start,end)
    qkey=defaultdict(set);positions=set()
    for r in qc:
        qkey[(int(r["pos"]),r["ref"],r["alt"])].add(r["first_failure"]);positions.add(int(r["pos"]))
    aligned=[]
    for key,entries in sorted(normalized.items()):
        chosen,removed=duplicates(entries)
        for e,reason in removed:e["event"]["reason"]=reason
        for e in chosen:
            rr=bykey.get(key,[])
            if len(rr)==1:
                e["reference"]=rr[0];e["event"]["reason"]="PASS";aligned.append(e)
            elif len(rr)>1:
                e["event"]["reason"]="reference_conflicting_duplicate"
            else:
                reasons=qkey.get(key,set())
                known=[r for r in ("reference_MAF_lt_0.05","reference_missing_gt_0.02","reference_HWE_midp_lt_1e-6") if r in reasons]
                e["event"]["reason"]=(known[0] if len(known)==1 else "reference_conflicting_duplicate" if "conflicting_duplicate" in reasons else
                    "reference_allele_mismatch" if not reasons and key[0] in positions else "reference_QC_or_absent")
    counts=Counter(e["reason"] for e in events)
    sig=Counter(e["reason"] for e in events if e["significant"])
    fields=list(dict.fromkeys(k for e in events for k in e))
    with gzip.open(output/"study_variant_qc.tsv.gz.part","wt") as f:
        writer=csv.DictWriter(f,fields,delimiter="\t");writer.writeheader();writer.writerows(events)
    (output/"study_variant_qc.tsv.gz.part").replace(output/"study_variant_qc.tsv.gz")
    m=len(aligned)
    result.update(n_variants=m,first_failure_counts=counts,stages=stages(counts,STUDY_STAGES),
                  significant_signal_counts=sig,significance_threshold=significant_threshold,
                  source_lead_retained=any(e["source_lead"] and e["reason"]=="PASS" for e in events),
                  zero_p_source=sum(float(r["pvalue"])==0 for r in rows),zero_p_retained=sum(e["pvalue"]==0 for e in aligned),
                  orientation_counts=Counter(e["original_to_final_orientation"] for e in aligned),
                  orientation_after_coordinate_mapping_counts=Counter(e["orientation"] for e in aligned),
                  mapping_failure_detail_counts=Counter(e['source_mapping_status'] for e in events if e['reason']=='mapping_failure'),
                  strong_signal_stratum=any(e["pvalue"]<1e-5 for e in aligned),
                  full_tested_variants_retained_after_QC=True)
    if m<2:
        result.update(state="INPUT_INSUFFICIENT",reason="fewer than two aligned QC-passing variants")
        atomic_json(receipt,result);return result
    if plan_only:
        result.update(state="PLANNED", dense_matrix_bytes=8*m*m,
                      dense_output_upper_bytes=16*m*m,
                      final_order_sha256=hashlib.sha256(json.dumps([
                          (e['pos'], e['ref'], e['alt'], int(e['reference']['variant_index']))
                          for e in aligned]).encode()).hexdigest())
        atomic_json(receipt,result)
        return result
    indices=np.array([int(e["reference"]["variant_index"]) for e in aligned],dtype=np.uint32)
    g=reference.read_alt_dosages(indices)
    f,stats=factor_from_dosages(g,maf_min=.05,missing_max=.02);del g
    for i,e in enumerate(aligned):
        if abs(stats["af"][i]-float(e["reference"]["af"]))>1e-12:raise ValueError("Reference dosage/order disagreement")
    ids=[f"chr{chrom}_{e['pos']}_{e['ref']}_{e['alt']}" for e in aligned]
    z=np.array([e["z"] for e in aligned],dtype="<f8")
    z.tofile(output/"z.bin.part")
    (output/"variant_ids.txt.part").write_text("\n".join(ids)+"\n")
    columns=["variant_id","pos","ref","alt","counted_allele","beta","se","z","pvalue","mlogp","study_maf","study_af",
             "aligned_study_af","original_id","original_pos","original_effect_allele","original_other_allele",
             "original_beta","original_se","flip_sign","orientation","original_to_final_orientation","mapping_strand","reference_af","reference_maf",
             "reference_missing_rate","reference_hwe_midp","reference_index","source_genotype_index"]
    with gzip.open(output/"variants.tsv.gz.part","wt") as out:
        writer=csv.DictWriter(out,columns,delimiter="\t");writer.writeheader()
        for i,e in enumerate(aligned):
            rr=e["reference"]
            record={k:e[k] for k in columns if k in e}
            record.update(variant_id=ids[i],counted_allele=e["alt"],reference_af=stats["af"][i],reference_maf=stats["maf"][i],
                          reference_missing_rate=stats["missing_rate"][i],reference_hwe_midp=rr["hwe_midp"],reference_index=indices[i],source_genotype_index=rr["source_index"])
            writer.writerow(record)
    dense_ld(f,output/"R.bin.part");del f
    # Reopen saved bytes and independently redecode all final genotype columns.
    f,_=factor_from_dosages(reference.read_alt_dosages(indices),maf_min=.05,missing_max=.02)
    validation=validate_matrix(output/"R.bin.part",f);del f
    column="R.bin"
    if not validation["bitwise_symmetric"]:
        column="R_F.bin"
        with (output/(column+".part")).open("wb") as out:
            for j in range(0,m,256):
                block=np.empty((min(256,m-j),m),dtype="<f8")
                for i in range(0,m,256):
                    a=row_block(output/"R.bin.part",m,i,min(m,i+256));block[:,i:i+256]=a[:,j:j+256].T;del a
                block.tofile(out)
    for name in ("z.bin","variant_ids.txt","variants.tsv.gz","R.bin")+((column,) if column!="R.bin" else ()):
        (output/(name+".part")).replace(output/name)
    parameters=common_parameters(settings,n,m)
    case=dict(case_id=case_id,n_variants=m,z_path="z.bin",R_path="R.bin",R_column_major_path=column,
              variant_ids_path="variant_ids.txt",parameters=parameters,provenance=provenance,qc_policy=POLICY,
              study_n=n,reference_n=reference.receipt["n_reference"],population=reference.receipt["population"],
              genome_build="GRCh38",source_build=provenance["source_build"],region=[chrom,start,end],
              orientation_basis="Final BED BIM A1 = source PVAR ALT; z=(original effect/SE)*recorded sign exactly once",
              orientation_definition="orientation is after coordinate-mapping strand transform; original_to_final_orientation is the net original-source to final counted-allele transform",
              source_code_sha256=identity["code"],variant_provenance_path="variants.tsv.gz",preparation_path="preparation.json")
    atomic_json(output/"case.json",case)
    validation["portable_reader"]="Numerical arrays and variant order checked directly"
    assert len(ids)==len(set(ids))
    for e in aligned:
        assert e["z"]==e["original_beta"]/e["original_se"]*e["flip_sign"]
        assert classify(e["ref"],e["alt"])=="PASS"
    validation.update(z_and_order_verified=True,all_variant_alleles_and_QC_verified=True)
    atomic_json(output/"validation.json",validation)
    names=list(dict.fromkeys(["case.json","variant_ids.txt","z.bin","R.bin",column,"variants.tsv.gz","study_variant_qc.tsv.gz","validation.json"]))
    result.update(state="READY",validation=validation,case_path=str(output/"case.json"),
                  hashes={name:sha256(output/name) for name in names},sizes={name:(output/name).stat().st_size for name in names},finished_epoch=time.time())
    atomic_json(receipt,result)
    print(json.dumps(dict(case=case_id,state="READY",p=m,n_ref=result["n_ref"])),flush=True)
    return result
