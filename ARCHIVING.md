# Manuscript benchmark archive

The benchmark archive version is **1.0.0**, with Git tag **analysis-v1.0.0**. This version identifies the analysis repository and its measurement records. The benchmarked software remains **prusie 0.3.0** and **susieR 0.16.6**; this archive does not release a new version of prusie.

## Reproduction scope

`make check` recomputes numerical summaries from the archived records and authenticates the frozen records and original execution scripts. `make all` additionally generates tables and redraws figures. These commands need the Level A Python environment and Cairo/librsvg for the architecture figure; they do not obtain association statistics, construct LD, install prusie or R, or rerun any fits. See [README.md](README.md), [environment/](environment/README.md), [01_inputs/](01_inputs/README.md), [02_benchmark/](02_benchmark/README.md), and [supplementary/](supplementary/README.md) for the separate externally sourced input-preparation and benchmark workflows.

The repository documents earlier full locus/gene selection checks and end-to-end checks on the first two GWAS loci and first two eQTL genes. These selected checks do not constitute an independent rerun of the complete benchmark campaign. The locally prepared eQTL LD census also has a documented numerical provenance limitation: its raw LD digests match the historical timing matrices in 0/14,157 cases. The census counts and historical runtime associations must be interpreted with that limitation.

## Current reporting definitions and frozen history

The primary manuscript speed-up is the ratio of the sums of per-case median complete-call times, with each case counted once. The all-call cumulative ratio and median of per-case ratios are separate statistics. The older `timing.headline` wording in [02_benchmark/protocol.json](02_benchmark/protocol.json) describes the all-call cumulative ratio; it does not define the current manuscript headline. The original protocol and historical execution archive are retained without editing. Current definitions are in [03_results/README.md](03_results/README.md) and [benchmark_manuscript_protocol.json](03_results/benchmark_manuscript_protocol.json).

The historical archive README links to the former `02_benchmark/benchmark_protocol.json` location. The restored original frozen protocol is now [02_benchmark/archive/as_run/protocol.json](02_benchmark/archive/as_run/protocol.json). Eight restored original copies under `archive/as_run/` match the original checksum manifest; unchanged originals remain in the parent archive. `02_benchmark/verify_archive.py` authenticates both locations. The historical files themselves are preserved.

Relative links in the restored historical `as_run/README.txt` reflect the original execution layout, before restoration into a subdirectory and removal of unreported pre-campaign qualification/smoke records. For current navigation, use [02_benchmark/README.md](02_benchmark/README.md), the parent [archive/README.txt](02_benchmark/archive/README.txt), and the restored protocol link above.

The execution archive intentionally retains original cluster paths, login and compute host names, scheduler user/group identifiers, allocation account names, and job records needed for provenance. These are historical execution details, not portable configuration or authentication credentials. Portable scripts use [config/example.yaml](config/example.yaml). No original association statistics, individual-level genotypes, or LD matrices are distributed. Pre-rendered manuscript figure PDFs are retained as reference artifacts; local redraws can differ in rendering.

## macOS validation environment

Use a separate Level A environment with the exact versions in `environment/requirements-level-a.txt`; the benchmark environment uses a different NumPy version. On Apple Silicon macOS with Cairo and librsvg already installed by Homebrew in `/opt/homebrew/lib`, the system `make` can remove a `DYLD_FALLBACK_LIBRARY_PATH` supplied only to its own environment. Set that variable in the Python command that `make` executes: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg make all PYTHON='DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib .venv/bin/python'`. This library-search setting changes no benchmark inputs, computations, or archived measurements. Other installations should use their actual library location.

The 23 November 2026 archival preparation used CPython 3.12.2 from the `py312` conda environment to create an isolated Level A virtual environment with the exact pinned dependencies. `make check` passed, matching all five regenerated summary/claim files, validating the full LD census, checking 20 frozen source digests, and authenticating 56 original execution-archive files (eight via restored copies). `make all` passed with the macOS library-search setting above, generating tables and running all eight figure scripts. The redrawn PDFs were kept as local validation artifacts; the tracked reference PDFs and all archived measurements were left byte-identical. `cffconvert --validate` passed for `CITATION.cff` schema 1.2.0. No benchmark fits were rerun during this preparation.

## Citation metadata and licensing

[CITATION.cff](CITATION.cff) supplies the title, ordered authors, manuscript-verified affiliations, description, keywords, archive version, and **GPL-3.0-or-later** license. The authors and affiliations follow the accompanying manuscript metadata; they are not inferred from the GitHub contributor list. No unverified ORCID, funding information, or manuscript DOI is included. The archived susieR files retain their BSD-3-Clause notices; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

The Zenodo GitHub integration reads `CITATION.cff`. A `.zenodo.json` file is unnecessary for this archive's metadata and is intentionally absent: Zenodo would otherwise ignore `CITATION.cff` and use only `.zenodo.json`. See Zenodo's [metadata guidance](https://help.zenodo.org/docs/github/describe-software/) and [supported CFF fields](https://help.zenodo.org/docs/github/describe-software/citation-file/).

After publication, cite the version-specific DOI assigned to `analysis-v1.0.0` for manuscript reproducibility. The concept DOI, if supplied by Zenodo, refers to all versions. A DOI badge and citation identifiers can subsequently be added to the default branch without creating another release. Obtain the exact archived Git commit from the release and record it alongside the DOI.
