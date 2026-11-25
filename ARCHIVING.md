# Manuscript benchmark archive

Benchmark archive version **1.0.1** (Git tag `analysis-v1.0.1`) captures the
analysis scripts and measurement records accompanying the manuscript evaluation
of **prusie 0.3.0** and **susieR 0.16.6**.

Version 1.0.1 updates the reproduction documentation, protocol descriptions, and result metadata. It carries forward the measurement records and numerical benchmark results from version 1.0.0.

## Reproduction workflows

- **Offline table and figure regeneration (Level A):** `make check` recomputes
  all numerical summaries from `03_results/` and verifies the SHA256 checksums
  of the archived measurement records and execution scripts. `make all` writes
  the summary tables to `build/tables/` and redraws all manuscript and
  supplementary figures. This workflow requires only the Level A Python
  environment (`environment/requirements-level-a.txt`) plus Cairo and `librsvg`
  for rendering `architecture.pdf`.
- **End-to-end benchmark execution from source data (Level B):** See
  [README.md](README.md), [environment/](environment/README.md),
  [01_inputs/](01_inputs/README.md), [02_benchmark/](02_benchmark/README.md),
  and [supplementary/](supplementary/README.md) for instructions on obtaining
  the external GWAS/eQTL summary statistics and reference genotypes, preparing
  aligned z-scores and PLINK LD matrices, and running the paired benchmark
  harness.

## Reporting definitions and execution archive

The primary manuscript speed-up is defined as the ratio of summed per-case
median complete-call runtimes ($\sum \text{median}(T_{\text{susieR}}) /
\sum \text{median}(T_{\text{prusie}})$), counting each analysis case once. The
cumulative ratio across all retained calls and the median of per-case speed-up
ratios are reported as companion summaries. Full definitions are given in
[03_results/README.md](03_results/README.md) and
[benchmark_manuscript_protocol.json](03_results/benchmark_manuscript_protocol.json).

`02_benchmark/archive/` preserves the scripts, patches, and environment records
from the original cluster execution. In the top-level `02_benchmark/archive/`
directory, cluster and host identifiers were redacted in eight files after the
campaign; their original pre-redaction bytes are stored in
[02_benchmark/archive/as_run/](02_benchmark/archive/as_run/), including
[02_benchmark/archive/as_run/protocol.json](02_benchmark/archive/as_run/protocol.json).
`02_benchmark/verify_archive.py` verifies all archived files against
`checksums.json`. For portable benchmark execution on new systems, use the
configurable drivers in `02_benchmark/` with [config/example.yaml](config/example.yaml).

## macOS validation environment

Create an isolated Python 3.12 virtual environment using the pinned package
versions in `environment/requirements-level-a.txt`. On Apple Silicon macOS with
Cairo and `librsvg` installed via Homebrew (`/opt/homebrew/lib`), pass
`DYLD_FALLBACK_LIBRARY_PATH` directly to the Python command invoked by `make`
(since macOS system `make` strips dynamic-linker variables from its own
environment):

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg \
  make all PYTHON='DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib .venv/bin/python'
```

## Citation metadata and licensing

[CITATION.cff](CITATION.cff) provides the archive title, author list and
affiliations, abstract, keywords, version, and **GPL-3.0-or-later** license in
Citation File Format (CFF 1.2.0) for GitHub and Zenodo indexing. Archived
`susieR` files retain their BSD-3-Clause notices; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

For manuscript reproducibility, cite the version-specific Zenodo archive
corresponding to release tag `analysis-v1.0.1`.
