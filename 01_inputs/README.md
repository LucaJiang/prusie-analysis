# Prepare the manuscript inputs

These scripts use the released FinnGen R13 T2D statistics, GTEx v11 Whole_Blood
all-pairs associations and GRCh38 high-coverage 1000 Genomes reference genotypes.
They write private association and LD files to the configured scratch directory.
Do not add those outputs to this repository.

Copy `config/example.yaml` to a local configuration and set each input path. Paths
are resolved relative to that configuration. Install the Level B dependencies in
`environment/requirements-inputs.txt`. PLINK 1.9 and PLINK 2 executables are external
inputs; their binaries are not distributed here.

```sh
python 01_inputs/ld_plink/prepare_reference.py --config config/local.yaml --stage samples
python 01_inputs/ld_plink/prepare_reference.py --config config/local.yaml --stage index
python 01_inputs/ld_plink/prepare_reference.py --config config/local.yaml --stage validate
python 01_inputs/ld_plink/prepare_reference.py --config config/local.yaml --stage prepare --chrom 1
python 01_inputs/gwas_finngen_t2d/select_loci.py --config config/local.yaml
# If the derived GENCODE TSS table is absent:
python 01_inputs/gtex_whole_blood/annotation.py --config config/local.yaml
python 01_inputs/gtex_whole_blood/select_genes.py --config config/local.yaml
python 01_inputs/gtex_whole_blood/extract_genes.py --config config/local.yaml --chrom 1
python 01_inputs/prepare_cases.py --config config/local.yaml --panel GWAS --chrom 1
python 01_inputs/prepare_cases.py --config config/local.yaml --panel eQTL --chrom 1
python 01_inputs/ld_plink/reconstruct_ld.py --config config/local.yaml --case-id T2D_chr1_14521917_17521917
```

Repeat chromosome preparation, extraction and alignment for chromosomes 1–22;
run PLINK reconstruction for each resulting case. The benchmark harness consumes
`<scratch>/plink_ld/<case_id>/case.json`. `--case-id` on extraction and alignment
restricts execution to named cases without changing their variants or QC.

`gtex_root` has the released-data layout below. Its annotation table is obtained
from the release-matched GENCODE 47 GTF: keep gene features, preserve versioned gene
IDs, and use the gene start on the positive strand or gene end on the negative
strand as the 1-based TSS. `select_genes.py` verifies every table entry against the GTF.

```text
gtex_root/
  annotations/gencode.v47.genes.gtf
  annotations/gencode.v47.gene_tss.tsv
  tissues/Whole_Blood/
    egenes/*.gz
    covariates/*.txt
    all_associations/Whole_Blood.v11.allpairs.chr1.parquet
```

Reference sample selection preserves the source order and checks the release's
first- and second-degree KING exclusions (kinship >0.0884), pedigree edges and
population metadata. It retains 99 FIN and 525 EUR individuals. Reference QC uses
PLINK's MAF ≥0.05, missingness ≤0.02 and exact HWE mid-P ≥10⁻⁶, excluding palindromes.
The original PGEN/PVAR and their lossless `.zst` copies are required by the source
validation step. The FASTA must contain the corresponding `chr1`–`chr22` sequences.

Association QC retains current-record biallelic A/C/G/T SNVs, excludes all A/T and
C/G pairs, and requires study MAF ≥0.05. Exact duplicates collapse to one record;
all conflicting duplicates are removed. Positions and alleles determine matching,
and the signed effect/SE is aligned once to the counted ALT allele. No P-value
filter is applied within a selected locus or cis window. Full numerical LD is
written without pruning, regularization, rounding or rank repair.

The GWAS locus selector preserves the published-P ordering (using released mlogp
for underflow), ±1.5 Mb lead windows, overlap/adjacency merging, recursive 10%
shrinking of regions exceeding 10 Mb, and final whole-window HLA exclusion.
GTEx selection uses q≤0.05, autosomes, inclusive TSS ±1 Mb and the same whole-window
HLA exclusion. Complete nominal row counts and signed float32 slope/SE promotion
are checked before alignment. Analysis n is 486,019 for GWAS and 800 for Whole_Blood;
the external-reference sample count is not used as analysis n.

`prepare_cases.py` first writes the historical normalized-genotype Gram matrix.
`reconstruct_ld.py` then applies the exact PLINK 1.9 commands used for the timed
campaign, restoring the frozen canonical order and counted alleles. These two LD
calculations can differ in floating-point entries. Numerical comparisons and the
local execution scope are recorded in the repository README.

Reference sample selection was also rerun from the official sample, pedigree and
kinship metadata; the FIN/EUR sample files matched their frozen hashes exactly.
The GENCODE TSS table was regenerated from the GTF and compared with its frozen
source table.

The full reference PVAR census, PGEN validation and chromosome-wide reference QC
were **not re-executed** for this repository. Their sanitized source is provided;
the four-case check reused the previously prepared strict FIN/EUR reference and
re-extracted the association records, aligned variants and reconstructed their LD.
Case JSON contains new provenance and source-code hashes, so its byte digest is
expected to differ from the historical manifest. The original frozen JSON digests
are checked separately; z and variant-order digests are compared directly.
