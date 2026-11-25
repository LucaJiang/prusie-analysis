# Prepare the manuscript inputs

The scripts in `01_inputs/` process FinnGen R13 T2D summary statistics, GTEx v11
Whole Blood all-pairs *cis*-eQTL associations, and GRCh38 high-coverage 1000
Genomes reference genotypes into aligned z-score vectors and dense float64 LD
matrices in the configured scratch directory.

Copy `config/example.yaml` to `config/local.yaml` and set the local input and
output paths (relative paths resolve against the configuration file directory).
Install the Python dependencies in `environment/requirements-inputs.txt` and
provide PLINK 1.9 and PLINK 2 executables.

```sh
python 01_inputs/ld_plink/prepare_reference.py --config config/local.yaml --stage samples
python 01_inputs/ld_plink/prepare_reference.py --config config/local.yaml --stage index
python 01_inputs/ld_plink/prepare_reference.py --config config/local.yaml --stage validate
python 01_inputs/ld_plink/prepare_reference.py --config config/local.yaml --stage prepare --chrom 1
python 01_inputs/gwas_finngen_t2d/select_loci.py --config config/local.yaml
# Generate the GENCODE v47 TSS table if not already present:
python 01_inputs/gtex_whole_blood/annotation.py --config config/local.yaml
python 01_inputs/gtex_whole_blood/select_genes.py --config config/local.yaml
python 01_inputs/gtex_whole_blood/extract_genes.py --config config/local.yaml --chrom 1
python 01_inputs/prepare_cases.py --config config/local.yaml --panel GWAS --chrom 1
python 01_inputs/prepare_cases.py --config config/local.yaml --panel eQTL --chrom 1
python 01_inputs/ld_plink/reconstruct_ld.py --config config/local.yaml --case-id T2D_chr1_14521917_17521917
```

Repeat `--stage prepare`, `extract_genes.py`, and `prepare_cases.py` across
chromosomes 1–22, then run `reconstruct_ld.py` for each case. The benchmark
harness reads `<scratch>/plink_ld/<case_id>/case.json`. Passing `--case-id` to
`extract_genes.py` or `prepare_cases.py` restricts execution to specific named
cases.

## Input layout and quality-control rules

`gtex_root` expects the directory structure below. `annotation.py` extracts
1-based strand-aware transcription start sites (gene start on `+`, gene end on
`-`) for `gene` features in the GENCODE v47 GTF, and `select_genes.py` verifies
the table against the GTF.

```text
gtex_root/
  annotations/gencode.v47.genes.gtf
  annotations/gencode.v47.gene_tss.tsv
  tissues/Whole_Blood/
    egenes/*.gz
    covariates/*.txt
    all_associations/Whole_Blood.v11.allpairs.chr1.parquet
```

- **Reference panel selection and QC:** Sample selection preserves the source
  PSAM order and applies the release's first- and second-degree KING exclusions
  (kinship $> 0.0884$), pedigree checks, and official IGSR population labels,
  retaining 99 FIN and 525 EUR individuals. Variant QC in PLINK requires
  $\text{MAF} \ge 0.05$, missingness $\le 0.02$, and exact HWE mid-$P \ge
  10^{-6}$, excluding palindromic SNVs (A/T and C/G).
- **Association QC and allele alignment:** Association QC retains biallelic
  non-palindromic A/C/G/T SNVs with study $\text{MAF} \ge 0.05$. Identical
  duplicate records are collapsed and conflicting duplicates are removed.
  Variants are matched by position and alleles, aligning the signed effect size
  to the counted reference ALT allele. All QC-passing variants within a selected
  locus or *cis* window are retained, and full float64 LD matrices are written
  without pruning or regularization.
- **Locus and eGene selection:** The GWAS selector orders genome-wide
  significant signals ($P < 5 \times 10^{-8}$) by published $P$ value (using
  `mlogp` on underflow), constructs $\pm 1.5\text{ Mb}$ lead windows, merges
  overlapping or adjacent intervals, recursively shrinks merged regions
  exceeding $10\text{ Mb}$ by 10%, and excludes windows overlapping the HLA
  region (`chr6:25,000,000-34,000,000`). The GTEx selector retains autosomal
  eGenes with $q \le 0.05$ and inclusive $\text{TSS} \pm 1\text{ Mb}$ windows
  outside the HLA region. Analysis sample sizes are $n = 486,019$ for GWAS and
  $n = 800$ for Whole Blood.
- **LD matrix construction:** `prepare_cases.py` aligns the association
  statistics and writes the initial case files; `reconstruct_ld.py` runs PLINK
  1.9 (`--r square bin`) to produce the signed float64 correlation matrix
  `R.bin` in canonical variant order.
