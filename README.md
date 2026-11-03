# prusie manuscript analyses

Code and published measurement records for the prusie manuscript: 200 FinnGen R13
type 2 diabetes loci, 14,157 GTEx v11 Whole Blood eGenes and six nested windows.
The [prusie package](https://github.com/LucaJiang/prusie) is maintained separately.

For the evaluated cases, the ratios of summed retained susieR/prusie complete-call
times are **49.20× for GWAS and 26.58× for eQTL**. The ratios of summed case medians
are **42.86× and 26.60×**; median paired case ratios are **47.74× and 27.53×**.
These definitions are distinct. Timings exclude imports, disk reads and LD
construction. Peak RSS includes the interpreter, libraries, loaded inputs, warmup
and retained fits; it is not fitting-only memory.

**Scientific scope.** The GWAS analysis n is **486,019** (89,727 cases and
396,292 controls), but its external LD uses **99 FIN reference individuals**.
A centered 99-sample genotype factor has rank at most 98; the first two checked
GWAS factors both have numerical rank 98. A full-panel rounded-LD rank census was
not performed. The eQTL analysis n is 800, with 525 EUR reference individuals.
These measurements support runtime and agreement comparisons on the evaluated
inputs. They are **not biological fine-mapping conclusions**.

Median maximum absolute alpha differences are 8.70 × 10⁻¹⁰ (GWAS) and
7.29 × 10⁻¹⁰ (eQTL). All 224 nonconverged eQTL cases remain in the summaries;
the largest alpha difference is 4.35 × 10⁻⁵. One converged eQTL case has a
2.43 × 10⁻⁴ PIP difference at the component activity threshold, despite an
alpha difference of 1.72 × 10⁻⁹. The saved credible-set membership comparison
reports identical membership across the primary panel. Alpha agreement does
not imply equality of every downstream quantity.

## Level A: tables and figures without fits

Use Python 3.12 and the pinned plotting environment:

```sh
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -r environment/requirements-level-a.txt
make check
make all
```

`make check` recomputes the primary/window TSVs, benchmark summary, derived
analyses and numerical claims, comparing decoded JSON and TSV values exactly.
It also authenticates the frozen inputs. `make tables` writes numerical outputs
to `build/tables/`. `make all` redraws eight manuscript
figures into `04_figures_tables/figures/`. No R, prusie, association data or fitting
is needed. CI runs Level A with Python 3.12.

| Manuscript output | Producing script in `04_figures_tables/` | Inputs in `03_results/` |
| --- | --- | --- |
| Figure 1: runtime and memory | `make_figure1.py` | `benchmark_records.jsonl.gz` |
| Figure 2: report | `make_report_figure_compact.py` | `report_snippets/` |
| Supplement: alpha/PIP vs variants and iterations; fitting time | `make_agreement_fitting_figures.py` | `benchmark_records.jsonl.gz` |
| Supplement: speed-up distribution | `make_speedup_figure.py` | `benchmark_primary.tsv` |
| Supplement: nested windows | `make_supplement_figures.py` | `benchmark_windows.tsv` |
| Supplement: full report figure | `make_report_figure.py` | `report_snippets/` |
| Dataset, timing, numerical consistency and memory tables; main-text numerical summaries | `summarize_benchmark.py`, `manuscript_numbers.py` | records, `primary_case_manifest.tsv`, `cs_membership.tsv` |
| Per-iteration, memory-by-size, LD-column classes and thread tables | `derived_analyses.py` | primary/window tables, `ld_class_census.tsv`, `thread_scaling.tsv` |
| Fixed settings and software identities | no calculation | `benchmark_manuscript_protocol.json`, `protocol_at_summary.json` |

`manuscript_numbers.json` exposes full precision values and both timing-ratio
definitions. The 10,000-variant R window uses the published ten-call replay in
`window10000_retime.json`; the original three calls remain in the record archive
and original benchmark summary. `protocol_at_summary.json` preserves the exact
protocol bytes identified by that summary. Figure checks compare panel content,
scales, labels and plotted patterns against the manuscript PDFs; pixel equality
is not required. A memory annotation is shortened to avoid clipping with fallback fonts. The report figures compose the supplied published screenshots.

## Level B: prepare inputs and rerun fits

1. Obtain the data below and set paths in a copy of `config/example.yaml`.
2. Follow [input preparation](01_inputs/README.md) for selection, strict SNP QC,
   alignment and PLINK LD reconstruction.
3. Follow the [paired benchmark](02_benchmark/README.md), using the frozen
   prusie 0.3.0 wheel and pinned susieR reference described in [environment](environment/README.md).
4. See [supplementary](supplementary/README.md) for the LD-class, credible-set,
   thread-scaling and nested-window producers.

No raw association statistics, genotypes or LD matrices are shipped. The
`02_benchmark/archive/` directory preserves exact execution evidence; use the
configurable scripts outside it for a new execution.

Local execution reproduced the complete locus/gene selection lists and the first
two GWAS loci and first two eQTL genes in manifest order. Their z-score and
variant-order hashes match exactly. Rebuilt PLINK LD matches the saved campaign
hashes for all four. Relative to the earlier normalized-genotype Gram matrices,
max |ΔR| is 5.88 × 10⁻¹⁵, 5.55 × 10⁻¹⁵, 1.54 × 10⁻¹⁴ and 1.62 × 10⁻¹⁴,
respectively. These are distinct floating-point LD constructions.

Reference sample selection and GENCODE TSS extraction were also rerun, with exact
frozen hashes. Full reference PVAR/PGEN validation and chromosome-wide reference
QC were **not re-executed**; the smoke reused prepared strict reference files.
Sanitized case JSON preserves scientific metadata and parameters but contains new
provenance paths and code identities, so its byte digest differs. The original
case-JSON digests were verified separately. Fitting-stage instrumentation, Slurm
execution and the full benchmark campaign were **not re-executed**. All eight final alpha arrays from the four local R/Python case pairs were
bitwise equal to their respective saved campaign arrays (max |Δalpha| = 0).
Between R and prusie, the four maximum differences were 1.60 × 10⁻⁹,
1.26 × 10⁻⁹, 6.58 × 10⁻⁹ and 7.43 × 10⁻¹³. Each backend ran one warmup
and three retained calls with one numerical thread and max_iter = 100.

## Data access

- **FinnGen R13**, type 2 diabetes: [FinnGen access](https://www.finngen.fi/en/access_results).
  Registration and acceptance of the access conditions are required, including
  the commitment not to identify participants. Statistics were available locally
  by 16 October 2026.
- **GTEx v11**, Whole Blood, GENCODE 47: [GTEx Portal](https://gtexportal.org/home/downloads/adult-gtex/qtl).
  Follow the [GTEx Portal Data License](https://gtexportal.org/home/license),
  acknowledging the source and acquisition date. These data were available locally
  by 18 October 2026; record your own acquisition date when downloading.
- **1000 Genomes**, GRCh38 high-coverage reference: [IGSR](https://www.internationalgenome.org/data).
  Follow the [IGSR disclaimer](https://www.internationalgenome.org/IGSR_disclaimer/),
  [EMBL–EBI terms](https://www.ebi.ac.uk/about/terms-of-use/) and collection-specific
  conditions. The reference was available locally by 16 October 2026.

The code is GPL-3.0-or-later. See [LICENSE](LICENSE) and
[third-party notices](THIRD_PARTY_NOTICES.md).
