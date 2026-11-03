"""Supplementary agreement and fitting-stage figures.

Writes:
    figures/supp_agreement.{pdf,svg,png}  — α-only numerical agreement
    figures/supp_agreement_niter.{pdf,svg,png} — the same against IBSS iterations
    figures/supp_fitting.{pdf,svg,png}    — fitting-stage runtime + aggregate speedup
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PAPER / ".cache/matplotlib"))
sys.path.insert(0, str(PAPER / "scripts"))

from _bench_plots import (
    fitting_by_variants_figure,
    numerical_agreement_figure,
)
from _bench_summaries import load_records

RECORDS = PAPER / "03_results/benchmark_records.jsonl.gz"
FIGURES = PAPER / "04_figures_tables/figures"


def main():
    records = load_records(RECORDS)
    numerical_agreement_figure(records, FIGURES, filename="supp_agreement")
    numerical_agreement_figure(records, FIGURES, filename="supp_agreement_niter", x="iterations")
    fitting_by_variants_figure(records, FIGURES, filename="supp_fitting")


if __name__ == "__main__":
    main()
