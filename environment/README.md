# Environments

## Level A: offline summaries, tables, and figures

Level A uses CPython 3.12 with `environment/requirements-level-a.txt`. Rendering
`architecture.svg` to `architecture.pdf` also uses the system libraries Cairo
and `librsvg` (`libcairo2` and `librsvg2-2` on Debian/Ubuntu). Keep the Level A
summary environment (NumPy 2.2.6) separate from the Level B benchmark
environment (NumPy 1.26.4), as tie-breaking and least-squares rounding can
differ across major NumPy versions.

## Level B: input preparation and paired benchmark execution

Input preparation (`01_inputs/`) uses `environment/requirements-inputs.txt`,
PLINK 1.9, and PLINK 2.

The paired benchmark (`02_benchmark/`) uses CPython 3.12, NumPy 1.26.4, PyYAML
6.0.3, and the `prusie` 0.3.0 Linux x86_64 wheel:

```text
prusie-0.3.0-cp312-cp312-manylinux_2_34_x86_64.whl
sha256 85285c26c1f6ae3caaed99c4d7766ede7f0083ea98eb2d26b7763d8654840b59
```

Download the wheel from the
[`prusie` v0.3.0 release](https://github.com/LucaJiang/prusie/releases/tag/v0.3.0),
verify its SHA256 digest, and install it with `pip`. The fitting-stage timing
campaign used an instrumented build of the same source commit
(`6d60d5e5723461c3e1fc5a46fd0a0e7b6378cfa7`, wheel SHA256
`ea54f79bae02240370811345ac34d751b6c58f4affd54477b31d649f90fd0af8`) created
with `02_benchmark/archive/native_timing.patch`.

The R reference environment uses **R 4.4.0** and **susieR 0.16.6** at commit
`8e56a8e038e989856d106d9ca5175cc664fea9d2`, together with `jsonlite`:

```r
lib <- Sys.getenv("ANALYSIS_R_LIBRARY")
stopifnot(nzchar(lib))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
install.packages(c("jsonlite", "remotes"), lib = lib)
remotes::install_github("stephenslab/susieR@8e56a8e038e989856d106d9ca5175cc664fea9d2",
                        lib = lib, upgrade = "never")
stopifnot(as.character(packageVersion("susieR")) == "0.16.6")
```

In the manuscript benchmark campaign, R linked against OpenBLAS 0.3.20 and
Python/NumPy used OpenBLAS 0.3.23.dev (recorded in
`02_benchmark/archive/environment.json`). Pin all numerical libraries to a
single thread for paired benchmark runs. A C compiler, R development headers,
and `libdl` are required when compiling `02_benchmark/blas_probe.c`.
