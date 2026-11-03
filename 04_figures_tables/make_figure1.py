"""Main-text Figure 1: per-case runtime (top) and backend peak memory (bottom).

Columns are the GWAS and eQTL panels. Every case is drawn; lines are medians
within log-spaced variant-count bins. Reads only data/benchmark_records.jsonl.gz.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from _bench_style import (
    COLORS,
    DARK,
    DOUBLE_COLUMN_IN,
    FIGURES,
    GIB,
    INK2,
    LABELS,
    MM,
    PAPER,
    RECORDS,
    binned_median,
    density_scatter,
    grid,
    load_records,
    median_line,
    note,
    panel_letter,
    plain_log_axis,
    primary,
    save,
    style,
)

PANELS = {
    "GWAS": {"title": "GWAS: 200 T2D loci (FinnGen R13)", "bins": 7, "min_bin": 5, "large_p": 10000},
    "eQTL": {"title": "eQTL: 14,157 genes (GTEx v11 whole blood)", "bins": 14, "min_bin": 20, "large_p": 5000},
}
BACKENDS = (("susieR", "R"), ("prusie", "prusie"))


def load_cases(records):
    cases = {}
    for panel in PANELS:
        rows = primary(records, panel)
        cases[panel] = {
            "p": np.array([r["p"] for r in rows], dtype=float),
            "time": {b: np.array([np.median(r[k]["times"]) for r in rows]) for b, k in BACKENDS},
            "time_sum": {b: float(sum(sum(r[k]["times"]) for r in rows)) for b, k in BACKENDS},
            "memory": {
                b: np.array([r[k]["peak_rss_diagnostic_bytes"] for r in rows]) / GIB
                for b, k in BACKENDS
            },
        }
    return cases


def draw(ax, panel, p, values, spec):
    order = np.argsort(p)
    for backend, _ in BACKENDS:
        y = values[backend]
        if panel == "GWAS":
            ax.scatter(p, y, s=6, color=COLORS[backend], alpha=0.45, linewidths=0, zorder=2)
        else:
            density_scatter(ax, p, y, COLORS[backend])
    for backend, _ in BACKENDS:
        x_bin, y_bin = binned_median(p[order], values[backend][order], spec["bins"], spec["min_bin"])
        median_line(ax, x_bin, y_bin, backend)
    ax.set_xscale("log")
    ax.set_yscale("log")
    grid(ax)


def main():
    style("prusie-figure1")
    cases = load_cases(load_records(RECORDS))
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COLUMN_IN, 120 * MM), sharex="col")
    fig.subplots_adjust(left=0.072, right=0.99, bottom=0.085, top=0.88, wspace=0.17, hspace=0.13)
    summary = {}
    for col, (panel, spec) in enumerate(PANELS.items()):
        c = cases[panel]
        p = c["p"]
        top, bottom = axes[0, col], axes[1, col]
        draw(top, panel, p, c["time"], spec)
        draw(bottom, panel, p, c["memory"], spec)

        line_p = np.geomspace(p.min(), p.max(), 200)
        bottom.plot(line_p, 8 * line_p**2 / GIB, color=INK2, lw=0.8, ls=(0, (4, 2)), zorder=3)

        m_r, m_p = c["memory"]["susieR"], c["memory"]["prusie"]
        ld = 8 * p**2
        large = p >= spec["large_p"]
        stats = {
            "cases": int(p.size),
            "cumulative_time_ratio": c["time_sum"]["susieR"] / c["time_sum"]["prusie"],
            "median_paired_memory_ratio": float(np.median(m_r / m_p)),
            "large_p": spec["large_p"],
            "large_p_cases": int(large.sum()),
            "large_p_median_peak_over_ld": {
                "prusie": float(np.median(m_p[large] * GIB / ld[large])),
                "susieR": float(np.median(m_r[large] * GIB / ld[large])),
            },
        }
        summary[panel] = stats
        note(top, f"Cumulative susieR/prusie runtime {stats['cumulative_time_ratio']:.1f}×")
        over = stats["large_p_median_peak_over_ld"]
        note(bottom, f"Median peak / LD ($p$ ≥ {spec['large_p']:,}): {over['prusie']:.2f} (prusie), "
                     f"{over['susieR']:.2f} (susieR)")

        top.set_title(spec["title"])
        bottom.set_xlabel("Variants in region, $p$")
        if panel == "GWAS":
            plain_log_axis(bottom.xaxis, [2000, 5000, 10000, 20000])
            bottom.set_xlim(1400, 22000)
            plain_log_axis(top.yaxis, [0.1, 1, 10, 100])
            top.set_ylim(0.02, 800)
            plain_log_axis(bottom.yaxis, [0.1, 0.2, 0.5, 1, 2, 5, 10, 20])
            bottom.set_ylim(0.1, 40)
            bottom.text(9000, 8 * 9000**2 / GIB * 0.8, "one LD matrix, $8p^2$ bytes", rotation=0,
                        ha="left", va="top", fontsize=6, color=INK2)
        else:
            plain_log_axis(bottom.xaxis, [10, 100, 1000])
            bottom.set_xlim(5, 8500)
            plain_log_axis(top.yaxis, [0.001, 0.01, 0.1, 1, 10, 100])
            top.set_ylim(5e-4, 400)
            plain_log_axis(bottom.yaxis, [0.05, 0.1, 0.2, 0.5, 1, 2, 5])
            bottom.set_ylim(0.05, 5)

    axes[0, 0].set_ylabel("Runtime per call (s)")
    axes[1, 0].set_ylabel("Peak memory (GiB)")
    for ax, letter in zip(axes.flat, "abcd", strict=True):
        panel_letter(ax, letter, x=-0.1 if ax in axes[:, 0] else -0.085, y=1.025)

    handles = [
        Line2D([], [], color=DARK["prusie"], marker="o", mfc="white", mec=DARK["prusie"], mew=1.0,
               lw=1.3, ms=3.6, label=LABELS["prusie"]),
        Line2D([], [], color=DARK["susieR"], marker="o", mfc="white", mec=DARK["susieR"], mew=1.0,
               lw=1.3, ms=3.6, label=LABELS["susieR"]),
        Line2D([], [], color=INK2, lw=0.8, ls=(0, (4, 2)), label="One dense float64 LD matrix"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.53, 1.0),
               handlelength=2.4, columnspacing=2.4)

    save(fig, FIGURES, "figure1")

    (FIGURES / "figure1_traceability.json").write_text(json.dumps({
        "source": "benchmark_records.jsonl.gz",
        "source_sha256": hashlib.sha256(RECORDS.read_bytes()).hexdigest(),
        "points": "One per primary case: median of three retained complete calls; backend-process peak RSS",
        "lines": "Median within log-spaced p bins with at least the stated case count",
        "runtime_annotation": "Sum of all retained susieR complete-call times / sum of all retained prusie times",
        "memory_annotation": "Median over cases with p >= large_p of backend-process peak RSS / 8p^2 bytes",
        "bins": {k: {"bins": v["bins"], "min_count": v["min_bin"]} for k, v in PANELS.items()},
        "summary": summary,
    }, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
