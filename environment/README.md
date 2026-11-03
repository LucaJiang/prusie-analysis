# Environments

Level A uses CPython 3.12 with `requirements-level-a.txt`. It needs no prusie or R.
Rendering `architecture.svg` also uses system libraries Cairo and librsvg
(`libcairo2` and `librsvg2-2` on Debian/Ubuntu). No Python binding is required.
The summary environment uses NumPy 2.2.6; the benchmark environment used NumPy
1.26.4. Keep these environments separate: sorting ties and least-squares rounding
can differ between NumPy versions.

Input preparation uses `requirements-inputs.txt`, PLINK 1.9 and PLINK 2.
The portable complete-call benchmark uses CPython 3.12, NumPy 1.26.4,
PyYAML 6.0.3 and the prusie 0.3.0 Linux CPython 3.12 release wheel:

```text
prusie-0.3.0-cp312-cp312-manylinux_2_34_x86_64.whl
sha256 85285c26c1f6ae3caaed99c4d7766ede7f0083ea98eb2d26b7763d8654840b59
```

Obtain it from the [0.3.0 release](https://github.com/LucaJiang/prusie/releases/tag/v0.3.0),
verify its SHA256, then install the file with pip. The original measured wheel
has SHA256 `ea54f79bae02240370811345ac34d751b6c58f4affd54477b31d649f90fd0af8`;
it adds stage timers and is not byte-identical to the release wheel.

The R reference is **R 4.4.0**, **susieR 0.16.6** at commit
`8e56a8e038e989856d106d9ca5175cc664fea9d2`, plus jsonlite and susieR's declared
dependencies. Install the source from that exact commit, not an unpinned HEAD:

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

This installs the pinned statistical source; dependency versions and BLAS must
also be recorded for a new timing experiment. The original R used OpenBLAS
0.3.20; Python used OpenBLAS 0.3.23.dev. The archived environment records contain
these identities. Set every numerical library to one thread for paired fits.
A C compiler, R headers and `libdl` are needed for the optional BLAS probe.
