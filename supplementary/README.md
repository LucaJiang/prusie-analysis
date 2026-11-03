# Supplementary analyses

The offline summaries are part of Level A:

```sh
python 04_figures_tables/derived_analyses.py
python 04_figures_tables/manuscript_numbers.py --output build/manuscript_numbers.json
```

The first command writes `03_results/derived_analyses.json`; `make check` uses a
temporary output directory to protect the frozen copy. The inputs are the published
per-case tables, LD-column census and thread-scaling table. No fit is run.

The following producers require local source inputs or saved fits. Set
`ANALYSIS_CONFIG=config/local.yaml`; their paths are in the configuration's
`supplementary` section.

| Script | Calculation | Local execution |
| --- | --- | --- |
| `ld_class_census.py GWAS 0 1 census.jsonl` | Exact signed column classes; all GWAS, every twentieth eQTL | First two cases of each panel recalculated; available published counts matched |
| `cs_membership.py 0 1 cs.jsonl` | Credible sets from both saved final alpha arrays and local LD, component IDs preserved | First two cases of each panel recalculated; counts and membership matched |
| `summarize_cs.py` | Merge complete CS shards | Not re-executed |
| `export_threads.py` | Export applicable C4 thread measurements to TSV | Full 546-row TSV regenerated exactly |
| `primary_case_manifest.py` | Extract metadata and input hashes from the frozen full case manifest | Full 14,357-row manifest regenerated exactly |

The LD census compares exact float64 rows up to sign; it counts column classes,
not matrix rank. CS extraction retains the original stable order, long-double
cumulative coverage, deduplication and minimum-purity rule. The full saved alpha
arrays and LD matrices needed for these producers are not distributed here.

The [archived thread measurement scripts](../02_benchmark/archive/thread_measurement/)
record the original fresh-process worker, paired-round calculations and schedule.
The manuscript thread table includes 182 cases with p≥1,024, one warmup and eight
paired rounds at 1, 2, 4 and 8 threads. Those thread measurements were not repeated.
The portable benchmark manifest includes all six first-p nested windows; no new
window measurements are implied by their offline regeneration.

The report figure scripts compose the published `03_results/report_snippets/`
screenshots from the susieR convergence toy example. They redraw the figures
without generating a new report or fitting a model.
