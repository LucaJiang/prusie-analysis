#!/usr/bin/env python3
"""Extract complete nominal gene rows in original chromosome-file order."""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,csv,sys
import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from _core.runtime import config
from _core.acquire import atomic_json,sha256

def extract(v,chrom,case_ids=None):
    root=Path(v['scratch']);source=Path(v['gtex_root'])/'tissues/Whole_Blood/all_associations'/f'Whole_Blood.v11.allpairs.chr{chrom}.parquet'
    with (root/'gtex_selection/frozen_genes.tsv').open() as stream:
        roster=[r for r in csv.DictReader(stream,delimiter='\t') if int(r['chrom'])==chrom and (not case_ids or r['gene_id'] in case_ids)]
    selected={r['gene_id']:r for r in roster};folder=root/'nominal'/f'chr{chrom}';folder.mkdir(parents=True,exist_ok=True)
    pf=pq.ParquetFile(source);total=0;counts=Counter();parts=defaultdict(list);values=pa.array(list(selected))
    for rg in range(pf.num_row_groups):
        table=pf.read_row_group(rg,use_threads=False);rgrows=len(table)
        table=table.append_column('source_row_0based',pa.array(np.arange(total,total+rgrows,dtype=np.int64)))
        fullcounts=Counter(table['gene_id'].to_pylist());counts.update(fullcounts)
        table=table.filter(pc.is_in(table['gene_id'],value_set=values));frame=table.to_pandas(use_threads=False)
        for gene,frame_gene in frame.groupby('gene_id',sort=False):
            path=folder/gene/f'rg{rg:03}.parquet';path.parent.mkdir(exist_ok=True)
            tab=pa.Table.from_pandas(frame_gene,schema=table.schema,preserve_index=False)
            pq.write_table(tab,path,compression='zstd');parts[gene].append(path)
        total+=rgrows;del frame,table
    assert total==pf.metadata.num_rows
    genes=[]
    for row in roster:
        gene=row['gene_id'];path=folder/gene/'complete.parquet';coverage=int(row['num_var'])==counts[gene]
        if parts[gene]:
            tables=[pq.read_table(p,use_threads=False) for p in parts[gene]];table=pa.concat_tables(tables)
            assert len(table)==counts[gene];pq.write_table(table,path,compression='zstd');del tables,table
        else:pq.write_table(pa.Table.from_batches([],schema=pf.schema_arrow),path)
        genes.append(dict(gene_id=gene,published_num_var=int(row['num_var']),nominal_rows=counts[gene],
                          complete_nominal_coverage=coverage,path=str(path),sha256=sha256(path)))
    atomic_json(root/'extraction'/f'chr{chrom}.json',dict(chrom=chrom,source_sha256=sha256(source),genes=genes))
    return genes

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--chrom',type=int,required=True);ap.add_argument('--case-id',action='append')
    args=ap.parse_args();extract(config(args.config)['inputs'],args.chrom,args.case_id)
if __name__=='__main__':main()
