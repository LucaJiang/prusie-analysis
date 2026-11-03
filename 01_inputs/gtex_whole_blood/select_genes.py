#!/usr/bin/env python3
"""Select GTEx v11 Whole_Blood eGenes and inclusive strand-aware cis windows."""
from pathlib import Path
from collections import Counter
import argparse,csv,gzip,re,sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from _core.runtime import config

from _core.acquire import atomic_json

def write_tsv(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or list(dict.fromkeys(k for r in rows for k in r))
    tmp=path.with_name(path.name+'.tmp')
    with tmp.open('w') as stream:
        writer=csv.DictWriter(stream,fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    tmp.replace(path)

def select(v):
    root=Path(v['gtex_root']);tissue=root/'tissues/Whole_Blood'
    with (root/'annotations/gencode.v47.gene_tss.tsv').open() as stream:
        tss={r['gene_id']:r for r in csv.DictReader(stream,delimiter='\t')}
    gtf_genes={}
    with (root/'annotations/gencode.v47.genes.gtf').open() as stream:
        for line in stream:
            if line.startswith('#'):continue
            a=line.rstrip().split('\t')
            if a[2]!='gene':continue
            attrs=dict(re.findall(r'(\w+) "([^"]*)"',a[8]));gid=attrs['gene_id']
            gtf_genes[gid]=(a[0],int(a[3]),int(a[4]),a[6]);r=tss[gid]
            assert (r['chromosome'],int(r['gene_start_1based']),int(r['gene_end_1based']),r['strand'])==gtf_genes[gid]
            assert int(r['tss_1based'])==int(a[3] if a[6]=='+' else a[4])
    assert len(gtf_genes)==len(tss)
    with gzip.open(next((tissue/'egenes').glob('*.gz')),'rt') as stream:
        egenes=list(csv.DictReader(stream,delimiter='\t'))
    assert len({r['gene_id'] for r in egenes})==len(egenes)
    roster=[];reasons=Counter()
    for source_index,r in enumerate(egenes):
        q=float(r['qval']);chrom=r['gene_chr'].removeprefix('chr');status='ELIGIBLE'
        if not np.isfinite(q):status='NONFINITE_Q'
        elif q>.05:status='NOT_EGENE_Q_GT_0.05'
        elif chrom not in set(map(str,range(1,23))):status='OUTSIDE_AUTOSOMES'
        elif r['gene_id'] not in tss:status='MISSING_GENCODE47'
        annotation=tss.get(r['gene_id']);pos=int(annotation['tss_1based']) if annotation else None
        start=max(1,pos-1000000) if pos else None;end=pos+1000000 if pos else None
        if annotation:
            assert (r['gene_chr'],int(r['gene_start']),int(r['gene_end']),r['strand'])==gtf_genes[r['gene_id']],r['gene_id']
        if status=='ELIGIBLE' and chrom=='6' and start<=34000000 and end>=25000000:status='WHOLE_WINDOW_HLA'
        reasons[status]+=1
        roster.append(dict(source_index=source_index,gene_id=r['gene_id'],gene_name=r['gene_name'],
            biotype=r.get('biotype',annotation.get('gene_type','') if annotation else ''),chrom=chrom,
            gene_start=r['gene_start'],gene_end=r['gene_end'],strand=r['strand'],tss=pos,start=start,end=end,
            qval=r['qval'],num_var=int(r['num_var']),status=status))
    roster.sort(key=lambda r:(int(r['chrom']) if r['chrom'].isdigit() else 99,r['tss'] or 0,r['gene_id']))
    eligible=[r for r in roster if r['status']=='ELIGIBLE']
    for i,r in enumerate(eligible):r['frozen_index']=i
    out=Path(v['scratch'])/'gtex_selection';write_tsv(out/'gene_census.tsv',roster);write_tsv(out/'frozen_genes.tsv',eligible)
    atomic_json(out/'selection.json',dict(selected=len(eligible),counts=dict(reasons)))
    return eligible

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True)
    args=ap.parse_args();select(config(args.config)['inputs'])
if __name__=='__main__':main()
