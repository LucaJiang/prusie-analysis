# prusie manuscript analyses

Citation metadata is provided in [CITATION.cff](CITATION.cff). See
[ARCHIVING.md](ARCHIVING.md) for archive structure, reproduction workflows, and
citation details (`analysis-v1.0.0`).

This repository contains the benchmark workflows, measurement records, and
figure- and table-generation code accompanying the manuscript comparing
[prusie](https://github.com/LucaJiang/prusie) 0.3.0 with `susieR` 0.16.6 across
200 FinnGen R13 type 2 diabetes (T2D) GWAS loci, 14,157 GTEx v11 Whole Blood
*cis*-eQTL eGenes, and six nested variant windows of a single locus.

All manuscript figures and summary tables can be regenerated directly from the
archived measurement records in `03_results/`. Re-running the fine-mapping fits
from scratch requires downloading the external source datasets (see
[Data access](#data-access)).

## Results at a glance

All timings were measured on exclusive AMD EPYC 7742 nodes using a single
numerical thread, with one warmup call and three retained timed calls per case
and backend.

| Metric | GWAS (200 loci) | eQTL (14,157 eGenes) |
| --- | ---: | ---: |
| Primary speed-up (summed `susieR` case medians / summed `prusie` case medians) | 42.86× | 26.60× |
| Cumulative speed-up across all retained calls | 49.20× | 26.58× |
| Median of per-case speed-up ratios | 47.74× | 27.53× |
| Median of per-case $\max |\Delta\alpha|$ | $8.7 \times 10^{-10}$ | $7.3 \times 10^{-10}$ |

- **Measurement boundaries:** Complete-call runtime measures a single
  `susie_rss` call on preloaded arrays (excluding module imports, disk I/O, and
  LD construction). Process peak memory records the high-water RSS of each
  fresh backend worker across input loading, one warmup call, and three timed
  calls.
- **Numerical agreement:** All 200 GWAS loci and 14,157 eQTL eGenes are
  included, including the 224 eQTL genes that reached the 100-iteration limit in
  both implementations (maximum $|\Delta\alpha| = 4.35 \times 10^{-5}$). One
  converged eQTL gene exhibits a PIP difference of $2.43 \times 10^{-4}$ due to
  a component near the purity/activity boundary ($\max |\Delta\alpha| =
  1.7 \times 10^{-9}$). Credible-set membership is identical between `prusie`
  and `susieR` across all 14,357 primary cases.
- **Sample sizes and LD column-class structure:** The GWAS analyses use
  $n = 486,019$ ($89,727$ cases and $396,292$ controls) with in-sample/reference
  LD estimated from 99 1000 Genomes FIN individuals; the eQTL analyses use
  $n = 800$ with LD from 525 EUR individuals. Grouping float64 LD columns that
  are identical up to sign yields $G$ exact equivalence classes per region. For
  eQTL, the median $G/p$ is $0.5934$ (IQR $0.5264\text{--}0.6586$), with 918
  genes in the comparison stratum $0.35 \le G/p < 0.45$ and
  $2,000 \le p \le 6,000$.

## Redraw figures and tables

Regenerating the summary tables and figures requires Python 3.12. Rendering the
software architecture diagram (`architecture.pdf`) also requires Cairo and
`librsvg` (`libcairo2` and `librsvg2-2` on Debian/Ubuntu).

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r environment/requirements-level-a.txt
make check     # recompute summaries, verify agreement with 03_results/, and check SHA256 digests
make all       # run make check, write tables to build/tables/, and render figures
```

`make check` verifies exact equality for non-numeric fields, relative agreement
within $10^{-12}$ for floating-point summaries across CPU/BLAS implementations,
and SHA256 checksums for the archived records and scripts.

| Manuscript item | Script in `04_figures_tables/` | Input in `03_results/` |
| --- | --- | --- |
| Fig. 1, runtime and memory | `make_figure1.py` | `benchmark_records.jsonl.gz` |
| Fig. 2, analysis report | `make_report_figure_compact.py` | `report_snippets/` |
| Software architecture | `make_architecture.py` | `architecture.svg` |
| Alpha/PIP agreement; fitting-stage time | `make_agreement_fitting_figures.py` | `benchmark_records.jsonl.gz` |
| Fig. S4, per-case speed-up distribution | `make_speedup_figure.py` | `benchmark_primary.tsv` |
| Fig. S6, peak memory | `make_memory_figure.py` | `benchmark_records.jsonl.gz` |
| Nested windows | `make_supplement_figures.py` | `benchmark_windows.tsv` |
| Full analysis report | `make_report_figure.py` | `report_snippets/` |
| Dataset, timing, agreement, and memory tables; manuscript statistics | `summarize_benchmark.py`, `manuscript_numbers.py` | `benchmark_records.jsonl.gz`, `primary_case_manifest.tsv`, `cs_membership.tsv` |
| Per-iteration, memory-by-size, LD column-class, and thread tables | `derived_analyses.py` | `benchmark_*.tsv`, `ld_class_census.tsv`, `thread_scaling.tsv` |
| Benchmark settings and software identities | — | `benchmark_manuscript_protocol.json` |

`manuscript_numbers.json` stores all summary statistics cited in the manuscript
at full precision. For the 10,000-variant nested window, `benchmark_windows.tsv`
uses the ten-call `susieR` replay in `window10000_retime.json` (Addressing a
single stalled original call), while `benchmark_records.jsonl.gz` retains the
initial three-call record.

## Rerun the analyses from source data

1. Obtain the source datasets listed below. Copy `config/example.yaml` to
   `config/local.yaml` and configure local paths.
2. Prepare the benchmark inputs using [`01_inputs/`](01_inputs/README.md) (locus
   and eGene selection, SNP quality control, allele alignment, and PLINK LD
   construction).
3. Execute the paired benchmark using [`02_benchmark/`](02_benchmark/README.md)
   with the pinned `prusie` 0.3.0 and `susieR` 0.16.6 environments described in
   [`environment/`](environment/README.md).
4. Run [`supplementary/`](supplementary/README.md) to generate the LD
   column-class census, credible-set membership comparisons, thread-scaling
   table, and primary case manifest.

`02_benchmark/archive/` preserves the exact execution scripts, timing
instrumentation patches, and SHA256 checksums from the manuscript benchmark
campaign.

## Data access

Raw association statistics, individual-level genotypes, and derived LD matrices
are not redistributed in this repository.

- **FinnGen R13** (type 2 diabetes summary statistics): Available via
  [FinnGen results access](https://www.finngen.fi/en/access_results) upon
  registration and acceptance of the data-access terms.
- **GTEx v11** (Whole Blood *cis*-eQTL summary statistics, GENCODE v47):
  Available from the [GTEx Portal](https://gtexportal.org/home/downloads/adult-gtex/qtl)
  under the [GTEx Portal Data License](https://gtexportal.org/home/license).
- **1000 Genomes** (GRCh38 high-coverage phased genotypes): Available from
  [IGSR](https://www.internationalgenome.org/data) under the
  [IGSR disclaimer](https://www.internationalgenome.org/IGSR_disclaimer/) and
  [EMBL-EBI terms of use](https://www.ebi.ac.uk/about/terms-of-use/).

The manuscript analyses used FinnGen R13 and 1000 Genomes data downloaded on
16 October 2026 and GTEx v11 data downloaded on 18 October 2026.

## License

The analysis code in this repository is licensed under GPL-3.0-or-later; see
[LICENSE](LICENSE). Archived `susieR` reference source files retain their
BSD-3-Clause license; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
