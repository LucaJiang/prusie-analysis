# prusie manuscript analyses

This repository holds the code and measurement records behind the prusie
manuscript. The manuscript compares [prusie](https://github.com/LucaJiang/prusie)
0.3.0 with susieR 0.16.6 on 200 FinnGen R13 type 2 diabetes loci, 14,157 GTEx v11
Whole Blood eGenes and six nested windows of one locus.

You can redraw every figure and table from the records here without R,
association data or new fits. Rerunning the fits themselves requires
downloading the source data (see [Data access](#data-access)).

## Results at a glance

All timings use exclusive AMD EPYC 7742 nodes and one numerical thread. Each
backend runs one warmup and three retained calls per case.

| | GWAS | eQTL |
| --- | ---: | ---: |
| Summed susieR time / summed prusie time, all retained calls | 49.20× | 26.58× |
| Same ratio using each case's median call | 42.86× | 26.60× |
| Median of per-case ratios | 47.74× | 27.53× |
| Median of per-case max \|Δα\| | 8.7 × 10⁻¹⁰ | 7.3 × 10⁻¹⁰ |

These ratios describe the evaluated cases. Complete-call time covers one
`susie_rss` call on loaded arrays; it excludes imports, disk reads and LD
construction. Peak memory is the whole backend process, including the
interpreter, loaded inputs and all calls, not just the fit.

Every case is kept, including the 224 eQTL genes where neither implementation
converged within 100 iterations. Among these, the largest alpha difference is
4.35 × 10⁻⁵. One converged eQTL gene has a PIP difference of 2.43 × 10⁻⁴,
caused by a component sitting at the activity threshold; its alpha difference
is 1.7 × 10⁻⁹. Credible-set membership is identical in every primary case.

### What these data can and cannot support

The GWAS fits pair an analysis sample of 486,019 (89,727 cases, 396,292
controls) with LD from only 99 FIN reference individuals, so the LD matrices
have rank at most 98. The eQTL fits use n = 800
with LD from 525 EUR individuals. The results here are valid comparisons of
runtime and numerical agreement on identical inputs. They are not biological
fine-mapping conclusions.

## Redraw figures and tables

You need Python 3.12. The architecture figure also needs Cairo and librsvg
(`libcairo2` and `librsvg2-2` on Debian/Ubuntu).

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r environment/requirements-level-a.txt
make check     # recompute summaries; values must equal those in 03_results/
make all       # also write tables to build/tables/ and redraw the figures
```

`make check` also verifies the checksums of the records and archived scripts.
CI runs the same command. Redrawn figures can differ from the manuscript PDFs in
fonts and pixels, but not in content.

| Manuscript item | Script in `04_figures_tables/` | Input in `03_results/` |
| --- | --- | --- |
| Fig. 1, runtime and memory | `make_figure1.py` | `benchmark_records.jsonl.gz` |
| Fig. 2, analysis report | `make_report_figure_compact.py` | `report_snippets/` |
| Software architecture | `make_architecture.py` | `architecture.svg` |
| Alpha/PIP agreement; fitting-stage time | `make_agreement_fitting_figures.py` | `benchmark_records.jsonl.gz` |
| Fig. S4, per-case speed-up distribution | `make_speedup_figure.py` | `benchmark_primary.tsv` |
| Nested windows | `make_supplement_figures.py` | `benchmark_windows.tsv` |
| Full analysis report | `make_report_figure.py` | `report_snippets/` |
| Dataset, timing, agreement and memory tables; numbers in the text | `summarize_benchmark.py`, `manuscript_numbers.py` | records, `primary_case_manifest.tsv`, `cs_membership.tsv` |
| Per-iteration, memory-by-size, LD column-class and thread tables | `derived_analyses.py` | `benchmark_*.tsv`, `ld_class_census.tsv`, `thread_scaling.tsv` |
| Fixed settings and software identities | — | `benchmark_manuscript_protocol.json` |

`manuscript_numbers.json` gives each number in the text at full precision.
The susieR result for the 10,000-variant window comes from a ten-call replay
(`window10000_retime.json`), because one original call stalled outside the
solver. The original calls remain in the records.

## Rerun the analyses

1. Download the data below. Copy `config/example.yaml` and set the paths.
2. Prepare the inputs with [`01_inputs/`](01_inputs/README.md). This covers locus
   and gene selection, strict SNP QC, allele alignment and LD reconstruction
   with PLINK.
3. Run the paired benchmark with [`02_benchmark/`](02_benchmark/README.md). It
   uses the prusie 0.3.0 wheel and the pinned susieR build described in
   [`environment/`](environment/README.md).
4. [`supplementary/`](supplementary/README.md) produces the LD column-class
   census, credible-set membership, thread-scaling table and case manifest.

`02_benchmark/archive/` holds the exact scripts and patches used for the
published measurements, with checksums. For new runs, use the configurable
scripts outside the archive.

**What has been re-executed from this repository.** Locus and gene selection,
reference sample lists and GENCODE TSS extraction have been rerun in full. The
first two GWAS loci and first two eQTL genes have been rerun end to end:
- The z-scores and variant order reproduced exactly.
- The PLINK LD matched the published hashes.
- The refitted alpha equalled the saved arrays.

Chromosome-wide reference QC, the full benchmark campaign, the timing
instrumentation build and the Slurm submission have not been re-executed from
this repository.

## Data access

No association statistics, genotypes or LD matrices are included.

- **FinnGen R13**, type 2 diabetes. Get it from
  [FinnGen results access](https://www.finngen.fi/en/access_results). You must
  register and accept the access conditions, including not attempting to
  identify participants.
- **GTEx v11**, Whole Blood cis-eQTL summary statistics (GENCODE 47). Get it from
  the [GTEx Portal](https://gtexportal.org/home/downloads/adult-gtex/qtl). Under
  the [GTEx Portal Data License](https://gtexportal.org/home/license), you must
  acknowledge the source and record the date you downloaded the data.
- **1000 Genomes**, GRCh38 high-coverage genotypes. Get it from
  [IGSR](https://www.internationalgenome.org/data). Use is subject to the
  [IGSR disclaimer](https://www.internationalgenome.org/IGSR_disclaimer/), the
  [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/) and any
  collection-specific conditions.

The published analyses used FinnGen and 1000 Genomes data obtained by
16 October 2026 and GTEx data obtained by 18 October 2026.

## License

GPL-3.0-or-later; see [LICENSE](LICENSE). The archived susieR source files keep
their BSD-3-Clause notice; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
