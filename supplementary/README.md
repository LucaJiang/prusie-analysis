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
| `ld_class_census.py PANEL SHARD NSHARDS OUT.jsonl` | Exact signed column classes for every primary case in `GWAS` or `eQTL`; no convergence filtering | Complete 200-GWAS and 14,157-eQTL census supplied; first two cases of each panel checked with this producer |
| `export_ld_census.py` | Merge complete JSONL shards into the public seven-column TSV, checking coverage and identities | Full-table export reproduced from completed census records |
| `cs_membership.py 0 1 cs.jsonl` | Credible sets from both saved final alpha arrays and local LD, component IDs preserved | First two cases of each panel recalculated; counts and membership matched |
| `summarize_cs.py` | Merge complete CS shards | Not re-executed |
| `export_threads.py` | Export applicable C4 thread measurements to TSV | Full 546-row TSV regenerated exactly |
| `primary_case_manifest.py` | Extract metadata and input hashes from the frozen full case manifest | Full 14,357-row manifest regenerated exactly |

The LD census compares exact float64 rows up to sign; it counts column classes,
not matrix rank. It canonicalizes the first nonzero sign and normalizes signed
zeros. Hashes screen candidates; final grouping requires exact array equality.
No approximate grouping or numerical tolerance is used. All 224 nonconverged
eQTL cases remain in the full primary panel.

The census uses locally prepared LD (99 FIN individuals for GWAS, 525 EUR for
eQTL). Local raw LD hashes match Burgundy timing LD hashes in 0/14,157 eQTL
cases. Exact classes in these local matrices therefore need not be the classes
in independently reconstructed timing matrices. Runtime associations use the
existing fitting-stage records and are descriptive; no new fits or timings
were run. See [provenance](../03_results/ld_class_census_provenance.json).

To recount from your configured prepared LD and export a full table:

```sh
python supplementary/ld_class_census.py GWAS 0 1 build/gwas_census.jsonl
python supplementary/ld_class_census.py eQTL 0 1 build/eqtl_census.jsonl
python supplementary/export_ld_census.py \
  --primary-table 03_results/benchmark_primary.tsv \
  --output build/ld_class_census.tsv build/gwas_census.jsonl build/eqtl_census.jsonl
```

For sharded execution use shard indices 0 through NSHARDS−1 and a distinct output
file per shard. Outputs resume by case key; use fresh files for a new recount.
The distributed table already contains the completed census, so offline
reproduction (`make tables`, `make check`) needs no LD matrices or recount.
Derived comparisons retain the fixed 0.35 ≤ G/p < 0.45, 2,000 ≤ p ≤ 6,000
stratum and stable G/p tertiles in case-key order, with no variant-count weighting.
Equal G/p values may straddle tertile boundaries.

CS extraction retains the original stable order, long-double
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
