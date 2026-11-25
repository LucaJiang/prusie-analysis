# Supplementary analyses

The offline summary scripts run as part of the Level A workflow using the
tables in `03_results/`:

```sh
python 04_figures_tables/derived_analyses.py
python 04_figures_tables/manuscript_numbers.py --output build/manuscript_numbers.json
```

`derived_analyses.py` writes `03_results/derived_analyses.json` by default
(or to `ANALYSIS_OUTPUT` when set, as in `make check` and `make tables`).

## Producers requiring local case inputs or saved fits

The scripts in `supplementary/` regenerate the supplementary tables in
`03_results/` from local case directories and saved benchmark outputs. Set
`ANALYSIS_CONFIG=config/local.yaml` with the `supplementary` paths configured.

| Script | Description |
| --- | --- |
| `ld_class_census.py PANEL SHARD NSHARDS OUT.jsonl` | Count exact signed float64 LD column equivalence classes ($G$) for each primary case in `GWAS` or `eQTL` |
| `export_ld_census.py` | Merge census JSONL shards into the seven-column `ld_class_census.tsv` table and verify panel coverage |
| `cs_membership.py SHARD NSHARDS OUT.jsonl` | Derive credible sets from saved final `alpha` arrays and local LD matrices using the `susieR` 0.16.6 rule |
| `summarize_cs.py` | Merge credible-set JSONL shards into `cs_membership.tsv` and `summary.json` |
| `export_threads.py` | Export multithreaded `prusie` 0.3.0 timing records to `thread_scaling.tsv` |
| `primary_case_manifest.py` | Extract case coordinates, sample sizes, and input SHA256 digests into `primary_case_manifest.tsv` |

### LD column-class census

`ld_class_census.py` groups columns of the float64 LD matrix that are identical
up to sign by canonicalizing the sign of the first nonzero entry, normalizing
signed zeros, screening with a 16-byte BLAKE2b hash, and verifying exact array
equality. To recompute the census across all primary cases and export the TSV
table:

```sh
python supplementary/ld_class_census.py GWAS 0 1 build/gwas_census.jsonl
python supplementary/ld_class_census.py eQTL 0 1 build/eqtl_census.jsonl
python supplementary/export_ld_census.py \
  --primary-table 03_results/benchmark_primary.tsv \
  --output build/ld_class_census.tsv build/gwas_census.jsonl build/eqtl_census.jsonl
```

For sharded execution, run shard indices `0` through `NSHARDS - 1` with distinct
output files per shard. Because `03_results/ld_class_census.tsv` is included in
the repository, offline table and figure generation (`make check`, `make all`)
uses the committed table directly.

### Credible-set membership and thread scaling

`cs_membership.py` applies the `susieR` 0.16.6 credible-set construction rule
(stable decreasing sort of component probabilities, long-double cumulative
sum to coverage $\ge 0.95$, row deduplication, and minimum absolute pairwise LD
purity $\ge 0.5$) identically to the saved `alpha` arrays of both backends.

The archived scripts in
[`02_benchmark/archive/thread_measurement/`](../02_benchmark/archive/thread_measurement/)
record the worker harness and balanced round schedule for the multithreaded
benchmark across 182 primary cases with $p \ge 1,024$ (one warmup call and eight
paired rounds at 1, 2, 4, and 8 threads).

### Analysis report figures

`04_figures_tables/make_report_figure.py` and
`04_figures_tables/make_report_figure_compact.py` compose the schematic diagram
and the pre-rendered report screenshots in `03_results/report_snippets/` (from
the `susieR` `FinemappingConvergence` dataset).
