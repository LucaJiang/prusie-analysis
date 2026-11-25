# Paired R/Python benchmark

`campaign.py` executes the paired benchmark comparing `prusie` 0.3.0 and
`susieR` 0.16.6 (commit `8e56a8e038e989856d106d9ca5175cc664fea9d2`). Each case
and backend runs in a fresh process with one warmup call and three retained
timed calls. Input arrays are preloaded before timing, and garbage collection
and result serialization occur outside the timed span. Backends run sequentially
within each case following a seeded balanced schedule (`seed = 20261028`).

After preparing the inputs with [`01_inputs/`](../01_inputs/README.md) and
configuring the `benchmark` section in [`config/example.yaml`](../config/example.yaml),
run:

```sh
python 02_benchmark/make_manifest.py --config config/local.yaml
python 02_benchmark/campaign.py --config config/local.yaml
```

`make_manifest.py` orders the 14,357 primary cases and appends the six nested
variant windows ($p \in \{500, 1000, 2000, 4000, 10000, 19221\}$ of
`T2D_chr10_61277948_71054950`), verifying the z-score and variant-order SHA256
digests against `03_results/primary_case_manifest.tsv`. To run a subset of
locally prepared cases, build the manifest with `make_manifest.py --allow-missing`
and pass `--case-keys` (for example,
`--case-keys GWAS_00000,GWAS_00001,eQTL_00000,eQTL_00001`) to `campaign.py`.

## Manifest format and configuration

The case manifest is a JSON list where each entry specifies `case_key`, `panel`
(`GWAS`, `eQTL`, or `windows`), and `case_path`, plus `window` and `kind` for
nested-window entries. Both drivers use the fitting parameters stored in each
`case.json` and export the final `alpha` and `pip` arrays, per-call runtimes,
convergence status, process peak RSS, and maximum absolute R/Python `alpha`
difference.

Configuration options in `benchmark` include `python`, `rscript`, optional
`r_build`, `cases_manifest`, `output_dir`, `r_libs`, `r_libs_user`, optional
`primary_case_manifest`, and optional `blas_probe`. If `blas_probe` is omitted,
the harness compiles `blas_probe.c` via `R CMD SHLIB` (requiring a C compiler,
R headers, and `libdl`).

## Fitting-stage timing instrumentation

By default (`benchmark.instrumented: false`), the portable harness times
complete `susie_rss` calls using the standard `prusie` 0.3.0 release wheel.
To record same-call fitting-stage times in addition to complete-call times,
build an instrumented `prusie` wheel with
[`archive/native_timing.patch`](archive/native_timing.patch) and set
`benchmark.instrumented: true`. When enabled, `stage_timer.py` records the
native solver and result-construction spans in Python, while `fit_timer.R`
verifies the installed `susieR` function bodies against the pinned sources in
`archive/R_source/` and records the `susie_workhorse` fitting span in R.

## Protocol and execution archive

- [`protocol.json`](protocol.json) records the benchmark protocol parameters. As
  described in [`ARCHIVING.md`](../ARCHIVING.md) and
  [`03_results/README.md`](../03_results/README.md), the primary manuscript
  speed-up is the ratio of summed per-case median complete-call times, whereas
  the legacy `timing.headline` field in `protocol.json` refers to the cumulative
  all-call ratio.
- [`archive/`](archive/README.txt) and [`archive/as_run/`](archive/as_run/README.md)
  preserve the exact scripts, patches, and environment manifests from the
  manuscript cluster campaign, authenticated by `verify_archive.py` against
  `checksums.json`.
- [`slurm.sh`](slurm.sh) provides a single-node Slurm job template for cluster
  execution.
