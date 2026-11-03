#!/usr/bin/env python3
"""Prepare the exact FIN/EUR sample selection and strict SNP reference."""
from pathlib import Path
import argparse,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from _core.runtime import config,paths
from _core.reference import select_samples,census_pvar,verify_local_sources
from _core.strict_reference import prepare

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,required=True)
    ap.add_argument('--stage',choices=['samples','index','validate','prepare'],required=True)
    ap.add_argument('--chrom',type=int)
    args=ap.parse_args();cfg=config(args.config);v=cfg['inputs'];p=paths(args.config)
    if args.stage=='samples':
        select_samples(Path(str(p['source'])+'.psam'),Path(v['pedigree']),Path(v['corrected_psam']),
                       Path(v['kinship_exclusions']),Path(v['kinship_pairs']),p['samples'])
    elif args.stage=='index':census_pvar(Path(str(p['source'])+'.pvar'),p['index'])
    elif args.stage=='validate':verify_local_sources(p['source'],p['index'])
    else:
        if args.chrom is None:ap.error('--chrom is required for prepare')
        prepare(args.config,args.chrom)
if __name__=='__main__':main()
