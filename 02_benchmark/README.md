# Paired R/Python benchmark

`campaign.py` runs the frozen z/LD cases with prusie 0.3.0 and susieR 0.16.6
(commit `8e56a8e038e989856d106d9ca5175cc664fea9d2`). Each case/backend uses a
fresh process, one warmup and three retained complete calls. Inputs are loaded
before the clock; garbage collection and export are outside it. Backends run
sequentially within each case. The original seeded balanced order is frozen
before selecting smoke cases. Independent cases may run concurrently.

Prepare cases with [01_inputs](../01_inputs/README.md), fill the `benchmark`
section in [example.yaml](../config/example.yaml), then run:

```sh
python 02_benchmark/make_manifest.py --config config/local.yaml
python 02_benchmark/campaign.py --config config/local.yaml
```

The manifest builder uses the published primary order and appends the six fixed
nested windows. It checks z and variant-order bytes against their frozen hashes
and preserves the full case parameters. For a four-case smoke with incomplete
local inputs, use `make_manifest.py --allow-missing`, then select
`--case-keys GWAS_00000,GWAS_00001,eQTL_00000,eQTL_00001` in `campaign.py`.
Backend order is still frozen from the complete manifest before that selection.

The case manifest is a JSON list. Each entry needs `case_key`, `panel`
(`GWAS`, `eQTL` or `windows`) and `case_path`. Paths resolve relative to the
manifest. A nested-window entry additionally supplies `window` and `kind`.
Case files retain all original fitting parameters; neither driver filters,
reorders nor changes associations. The harness writes final full alpha and
PIP arrays, per-call timing and convergence, process high-water RSS and maximum
absolute R/Python alpha differences. Alpha rows keep their original order.
No component is excluded or permuted.

The `python`, `rscript`, optional `r_build`, `cases_manifest`, `output_dir`,
`r_libs`, `r_libs_user`, optional `primary_case_manifest` and optional `blas_probe` paths come from configuration.
Relative paths resolve against that configuration. Output must be a new directory.
Without a supplied probe, the harness builds `blas_probe.c` using `R CMD SHLIB`.
The probe requires a C compiler, R development headers and Linux `libdl`.
For a custom Rscript wrapper, set `r_build` to the matching `R` build executable.
`r_libs` and `r_libs_user` may be lists of library directories. The default
portable mode measures complete calls with the frozen released wheel. It does
not report a fitting-stage time as if that instrumented wheel were installed.

The manuscript fitting-stage measurements used a timing-instrumented derivative
of the frozen wheel, plus the guarded R timer in `fit_timer.R`. Enable
`benchmark.instrumented: true` only with a separate build made using the exact
[archived patch](archive/native_timing.patch). `stage_timer.py` requires the
two native timing keys added by that patch and rejects an ordinary release
wheel. No statistical operation is changed by the timing helpers. The R helper
checks the installed function bodies against the pinned source excerpts.
The patch build and fitting-stage instrumentation are not re-executed here.

[protocol.json](protocol.json) records the portable measurement definitions.
For the evaluated cases, the headline speed-up is the sum of all retained R
times divided by the corresponding prusie total. Ratios of case medians and
geometric mean ratios are separate summaries. The manuscript includes all
fixed cases, including iteration-limit fits and alpha discrepancies.

[archive/](archive/README.txt) preserves the published execution archive verbatim,
including the original `checksums.json`. Eight published copies had subsequently
redacted names, as that README explains. Their original bytes, recovered from
Git history and verified against the retained hashes, are in
[archive/as_run/](archive/as_run/README.md); unchanged files stay in the parent
archive. `archive/as_run/protocol.json` preserves the original frozen protocol.
The archive contains cluster-specific paths and operational controls for historical
traceability. The portable drivers omit those controls and do not submit jobs.
[slurm.sh](slurm.sh) is an editable single-node template; it has not been submitted.

Local smoke execution and numerical differences are documented in the repository
README. They are descriptive reproduction checks on four fixed cases; they do not
replace the published complete-panel measurements.
