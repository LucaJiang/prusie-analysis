#!/usr/bin/env python3
"""Generate the release-matched strand-aware GENCODE 47 gene/TSS table."""
from pathlib import Path
import argparse,re,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from _core.runtime import config
from _core.finngen import write_tsv

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True)
    args=ap.parse_args();root=Path(config(args.config)['inputs']['gtex_root'])/'annotations'
    src=root/'gencode.v47.genes.gtf';dest=root/'gencode.v47.gene_tss.tsv'
    if dest.exists():raise FileExistsError('Use the existing table, or choose a new annotation directory')
    rows=[];seen=set()
    with src.open() as f:
        for line in f:
            if line.startswith('#'):continue
            fields=line.rstrip('\n').split('\t');assert len(fields)==9
            if fields[2]!='gene':continue
            attributes=dict(re.findall(r'(\w+) "([^"]*)"',fields[8]))
            gene=attributes['gene_id'];assert gene not in seen;seen.add(gene)
            assert fields[6] in ('+','-')
            rows.append({'gene_id':gene,'gene_name':attributes.get('gene_name',''),'chromosome':fields[0],'gene_start_1based':fields[3],
                'gene_end_1based':fields[4],'strand':fields[6],'tss_1based':fields[3] if fields[6]=='+' else fields[4],
                'gene_type':attributes.get('gene_type','')})
    assert rows
    write_tsv(dest,rows,list(rows[0]))

if __name__=='__main__':main()
