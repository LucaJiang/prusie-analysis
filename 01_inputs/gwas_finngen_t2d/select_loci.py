#!/usr/bin/env python3
"""Select all autosomal FinnGen R13 T2D loci before examining fit outcomes."""
from pathlib import Path
import argparse,sys
import pysam
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from _core.runtime import config
from _core.finngen import census

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True)
    args=ap.parse_args();v=config(args.config)['inputs']
    with pysam.FastaFile(v['fasta']) as fa:
        lengths={str(c):fa.get_reference_length('chr'+str(c)) for c in range(1,23)}
    census(Path(v['finngen_summary']),lengths,Path(v['scratch'])/'gwas_selection')
if __name__=='__main__':main()
