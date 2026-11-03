"""Supplementary figure: nested-window scaling.

Writes:
    figures/supp_scaling.{pdf,svg,png}   — nested-window runtime and iterations

The p = 10,000 susieR window entries come from the dedicated replay already
recorded in data/benchmark_windows.tsv.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PAPER / ".cache/matplotlib"))
sys.path.insert(0, str(PAPER / "scripts"))

import numpy as np

from _bench_style import (
    COLORS,
    DARK,
    DOUBLE_COLUMN_IN,
    HALO,
    HIGHLIGHT,
    INK,
    INK2,
    LABELS,
    MM,
    grid,
    note,
    panel_letter,
    plain_log_axis,
    save,
    style,
)
from matplotlib import pyplot as plt
from matplotlib import ticker
from matplotlib.lines import Line2D

DATA = PAPER / "03_results"
FIGURES = PAPER / "04_figures_tables/figures"


# ---------------------------------------------------------------------------
# supp_scaling — nested windows
# ---------------------------------------------------------------------------

WINDOW_P = [500, 1000, 2000, 4000, 10000, 19221]


def scaling():
    with (DATA / "benchmark_windows.tsv").open() as handle:
        rows = sorted(csv.DictReader(handle, delimiter="\t"), key=lambda row: int(row["p"]))
    assert [int(row["p"]) for row in rows] == WINDOW_P
    # The p = 10,000 susieR entries in benchmark_windows.tsv are the dedicated replay medians.
    replay = next(row for row in rows if int(row["p"]) == 10000)
    assert abs(float(replay["susieR_seconds"]) - 59.2686502905563) < 1e-9
    p = np.array(WINDOW_P, dtype=float)

    style("prusie-scaling")
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COLUMN_IN, 62 * MM))
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.18, top=0.8, wspace=0.34)
    ax_complete, ax_fitting, ax_iter = axes

    for ax, prefix, title in (
        (ax_complete, "", "Complete call"),
        (ax_fitting, "fitting_", "Fitting stage"),
    ):
        for backend in ("susieR", "prusie"):
            med = np.array([float(row[f"{backend}_{prefix}seconds"]) for row in rows])
            q1 = np.array([float(row[f"{backend}_{prefix}q1"]) for row in rows])
            q3 = np.array([float(row[f"{backend}_{prefix}q3"]) for row in rows])
            ax.fill_between(p, q1, q3, color=COLORS[backend], alpha=0.25, lw=0, zorder=2)
            ax.plot(p, med, color=DARK[backend], lw=1.3, marker="o", ms=3.4, mfc="white",
                    mec=DARK[backend], mew=1.0, path_effects=HALO, zorder=4)
        ax.set_xscale("log")
        ax.set_yscale("log")
        plain_log_axis(ax.yaxis, [0.01, 0.1, 1, 10, 100])
        ax.set_ylim(0.005, 1000)
        ax.set_ylabel("Time per call (s)")
        ax.set_title(title)
        grid(ax)

    iters = np.array([int(row["prusie_niter"]) for row in rows])
    assert np.array_equal(iters, [int(row["susieR_niter"]) for row in rows])
    ax_iter.plot(p, iters, color=INK2, lw=1.1, marker="o", ms=3.4, mfc="white", mec=INK2,
                 mew=1.0, zorder=4)
    for x, n in zip(p, iters, strict=True):
        ax_iter.annotate(str(n), (x, n), xytext=(0, 4), textcoords="offset points",
                         ha="center", va="bottom", fontsize=6, color=INK)
    ax_iter.set_xscale("log")
    ax_iter.set_ylim(0, 34)
    ax_iter.set_ylabel("Completed IBSS iterations")
    ax_iter.set_title("Iterations (identical for both)")
    grid(ax_iter)

    for ax, letter in zip(axes, "abc", strict=True):
        plain_log_axis(ax.xaxis, WINDOW_P)
        ax.xaxis.set_major_formatter(
            ticker.FixedFormatter(["500", "1k", "2k", "4k", "10k", "19.2k"])
        )
        ax.set_xlim(400, 25000)
        ax.set_xlabel("Leading variants, $p$")
        panel_letter(ax, letter, x=-0.24, y=1.025)
    handles = [
        Line2D([], [], color=DARK[b], marker="o", mfc="white", mec=DARK[b], mew=1.0, lw=1.3,
               ms=3.6, label=LABELS[b])
        for b in ("prusie", "susieR")
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.53, 1.0),
               handlelength=2.4, columnspacing=2.4)
    save(fig, FIGURES, "supp_scaling")


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    scaling()
    ledger_path = DATA / "figure_sources.json"
    ledger = json.loads(ledger_path.read_text())
    for name, source in (
        ("supp_scaling.pdf", "benchmark_windows.tsv"),
    ):
        ledger[name] = {
            "sha256": _hash(FIGURES / name),
            "source": source,
            "source_sha256": _hash(DATA / source),
            "generator": "scripts/make_supplement_figures.py",
            "hardware": "exclusive AMD EPYC 7742 compute nodes",
        }
    ledger["figure1.pdf"] = {
        "sha256": _hash(FIGURES / "figure1.pdf"),
        "source": "benchmark_records.jsonl.gz and figure1_traceability.json",
        "generator": "scripts/make_figure1.py",
        "hardware": "exclusive AMD EPYC 7742 compute nodes",
    }
    for supp in ("supp_agreement.pdf", "supp_fitting.pdf"):
        if (FIGURES / supp).exists():
            ledger[supp] = {
                "sha256": _hash(FIGURES / supp),
                "source": "benchmark_records.jsonl.gz",
                "source_sha256": _hash(DATA / "benchmark_records.jsonl.gz"),
                "generator": "scripts/make_agreement_fitting_figures.py",
                "hardware": "exclusive AMD EPYC 7742 compute nodes",
            }
    (FIGURES / "figure_sources.json").write_text(json.dumps(ledger, indent=2) + "\n")


if __name__ == "__main__":
    main()
