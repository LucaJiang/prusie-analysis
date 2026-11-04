"""Supplementary peak-memory figure.

Writes 04_figures_tables/figures/supp_memory.{pdf,svg,png}.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PAPER / ".cache/matplotlib"))
sys.path.insert(0, str(PAPER / "scripts"))

from _bench_plots import peak_memory_figure
from _bench_summaries import load_records

RECORDS = PAPER / "03_results/benchmark_records.jsonl.gz"
FIGURES = PAPER / "04_figures_tables/figures"


def main():
    peak_memory_figure(load_records(RECORDS), FIGURES)


if __name__ == "__main__":
    main()
