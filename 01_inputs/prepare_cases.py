#!/usr/bin/env python3
"""Align all variants in selected windows and apply the manuscript's strict SNP QC."""
from pathlib import Path
import argparse,csv,gzip,json
import pyarrow.parquet as pq
import pysam
from _core.runtime import config
from _core.acquire import sha256
from _core.strict_prepare import prepare_case
from _core.strict_reference import StrictReferenceReader
from _core.strict_policy import classify

def gwas(cfg,chrom,case_ids=None):
    v=cfg['inputs'];root=Path(v['scratch'])
    with (root/'gwas_selection/gwas_regions.tsv').open() as stream:
        regions=[r for r in csv.DictReader(stream,delimiter='\t') if int(r['chrom'])==chrom and r['status']=='ELIGIBLE' and (not case_ids or r['case_id'] in case_ids)]
    reference=StrictReferenceReader(Path(v['reference_root'])/'FIN'/f'chr{chrom}')
    with pysam.FastaFile(v['fasta']) as fa:sequence=fa.fetch(f'chr{chrom}').upper()
    summary=Path(v['finngen_summary']);summaryhash=sha256(summary)
    manifest=list(csv.DictReader(Path(v['finngen_manifest']).open(),delimiter='\t'))
    metadata=next(r for r in manifest if r.get('phenocode',r.get('phenotype',''))=='T2D')
    n_cases=int(metadata['num_cases']);n_controls=int(metadata['num_controls']);n=n_cases+n_controls
    with pysam.TabixFile(str(summary)) as tab:
        header=next(h for h in tab.header if h.startswith('#chrom')).split('\t')
        for region in regions:
            start,end=int(region['start']),int(region['end']);rows=[]
            for i,line in enumerate(tab.fetch(str(chrom),start-1,end)):
                raw=dict(zip(header,line.split('\t'),strict=True))
                rows.append(dict(pos=raw['pos'],original_pos=raw['pos'],original_chrom=str(chrom),other=raw['ref'],effect=raw['alt'],
                    beta=raw['beta'],se=raw['sebeta'],af=raw['af_alt'],pvalue=raw['pval'],mlogp=raw['mlogp'],original_id=raw.get('rsids',''),source_row=i))
            lead=min(rows,key=lambda r:float(r['pvalue']) if float(r['pvalue'])>0 else -float(r['mlogp']))['original_id']
            lead_row=min(range(len(rows)),key=lambda i:float(rows[i]['pvalue']) if float(rows[i]['pvalue'])>0 else -float(rows[i]['mlogp']))
            for i,row in enumerate(rows):row['source_lead']=i==lead_row
            prov=dict(source='FinnGen',release='R13',endpoint='T2D',source_build='GRCh38',final_build='GRCh38',population='FIN',
                summary_sha256=summaryhash,source_n=n,n_cases=n_cases,n_controls=n_controls,
                n_convention='N_cases + N_controls; Gaussian RSS approximation to binary trait',effect_allele_contract='published ALT',
                reference_receipt_sha256=sha256(str(reference.prefix)+'.receipt.json'),region_census_sha256=sha256(root/'gwas_selection/gwas_regions.tsv'))
            prepare_case(rows,case_id=region['case_id'],chrom=str(chrom),start=start,end=end,n=n,reference=reference,
                sequence=sequence,output=root/'cases'/region['case_id'],settings=cfg['scientific_settings'],provenance=prov,significant_threshold=5e-8,lead_id=lead)
    reference.close()

def eqtl(cfg,chrom,case_ids=None):
    v=cfg['inputs'];root=Path(v['scratch']);tissue=Path(v['gtex_root'])/'tissues/Whole_Blood'
    extract=json.loads((root/'extraction'/f'chr{chrom}.json').read_text())
    with (root/'gtex_selection/frozen_genes.tsv').open() as stream:
        roster={r['gene_id']:r for r in csv.DictReader(stream,delimiter='\t') if int(r['chrom'])==chrom}
    with gzip.open(next((tissue/'egenes').glob('*.gz')),'rt') as stream:eg={r['gene_id']:r for r in csv.DictReader(stream,delimiter='\t')}
    # The released tissue covariate header is the analysis sample roster.
    covpath=next((tissue/'covariates').glob('*.txt'))
    with covpath.open() as stream:n=len(next(stream).rstrip('\n').split('\t'))-1
    if n!=800:raise ValueError('Expected the release-matched 800-subject Whole_Blood cohort')
    with pysam.FastaFile(v['fasta']) as fa:sequence=fa.fetch('chr'+str(chrom)).upper()
    ref=StrictReferenceReader(Path(v['reference_root'])/'EUR'/f'chr{chrom}')
    for entry in extract['genes']:
        gene=entry['gene_id']
        if case_ids and gene not in case_ids:continue
        row=roster[gene]
        if not entry['complete_nominal_coverage']:raise ValueError('Incomplete nominal coverage')
        p=Path(entry['path']);assert sha256(p)==entry['sha256']
        raw=pq.read_table(p,use_threads=False).to_pylist();rows=[];dist_errors=0;window_errors=0;ref_errors=0
        for r in raw:
            parts=r['variant_id'].split('_')
            if len(parts)!=5 or parts[0]!='chr'+str(chrom) or parts[4]!='b38':raise ValueError('Unexpected released variant schema')
            pos=int(parts[1]);other,effect=parts[2:4]
            dist_errors+=r['tss_distance']!=pos-int(row['tss']);window_errors+=not(int(row['start'])<=pos<=int(row['end']))
            valid_fasta=sequence[pos-1:pos-1+len(other)]==other
            if classify(other,effect)=='PASS' and not valid_fasta:ref_errors+=1
            rows.append(dict(original_id=r['variant_id'],pos=pos,original_pos=pos,original_chrom=str(chrom),other=other,effect=effect,
                beta=r['slope'],se=r['slope_se'],af=r['af'],pvalue=r['pval_nominal'],source_fasta_valid=valid_fasta,
                source_row=r['source_row_0based'],source_lead=r['variant_id']==eg[gene]['variant_id']))
        if dist_errors or window_errors:raise ValueError('TSS/cis-window disagreement')
        prov=dict(source='GTEx',release='v11',tissue='Whole_Blood',source_build='GRCh38',final_build='GRCh38',gene_id=gene,
            gene_name=row['gene_name'],biotype=row['biotype'],strand=row['strand'],tss=int(row['tss']),source_n=n,
            allpairs_sha256=extract['source_sha256'],nominal_gene_sha256=entry['sha256'],source_num_var=int(row['num_var']),
            complete_nominal_rows=len(rows),complete_nominal_coverage=True,
            source_float_precision='af/slope/slope_se float32 promoted exactly to float64; pval_nominal float64',
            tss_distance_errors=dist_errors,outside_cis_window=window_errors,nonpalindromic_SNV_REF_mismatch=ref_errors,
            reference_population='EUR',effect_allele_contract='released ALT; signed slope/slope_se, aligned exactly once to reference ALT',
            reference_receipt_sha256=sha256(str(ref.prefix)+'.receipt.json'),gene_roster_sha256=sha256(root/'gtex_selection/frozen_genes.tsv'))
        prepare_case(rows,case_id=gene,chrom=str(chrom),start=int(row['start']),end=int(row['end']),n=n,reference=ref,sequence=sequence,
            output=root/'cases'/gene,settings=cfg['scientific_settings'],provenance=prov,significant_threshold=float(eg[gene]['pval_nominal_threshold']),lead_id=eg[gene]['variant_id'])
    ref.close()

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--panel',choices=['GWAS','eQTL'],required=True);ap.add_argument('--chrom',type=int,required=True)
    ap.add_argument('--case-id',action='append');args=ap.parse_args()
    (gwas if args.panel=='GWAS' else eqtl)(config(args.config),args.chrom,args.case_id)
if __name__=='__main__':main()
