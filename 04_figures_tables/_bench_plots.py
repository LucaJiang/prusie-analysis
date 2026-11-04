"""Supplementary agreement and fitting-stage figures from benchmark_records.jsonl.gz."""

from __future__ import annotations

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
from matplotlib.lines import Line2D

from _bench_style import (
    COLORS,
    DARK,
    DOUBLE_COLUMN_IN,
    GIB,
    HALO,
    HIGHLIGHT,
    INK,
    INK2,
    LABELS,
    MM,
    MUTED,
    binned_median,
    density_scatter,
    grid,
    median_line,
    note,
    panel_letter,
    plain_log_axis,
    primary,
    save,
    style,
)

PANELS = ("GWAS", "eQTL")
TITLES = {"GWAS": "GWAS: 200 T2D loci", "eQTL": "eQTL: 14,157 genes"}
BINS = {"GWAS": (7, 5), "eQTL": (14, 20)}
FLOOR = 1e-17
ZERO_ROW = 10**-18.6
CASE_POINT = "#4d5d6c"


def sci(value):
    mantissa, exponent = f"{value:.2g}".split("e")
    return f"{mantissa}×10$^{{{int(exponent)}}}$"


def _error_axis(ax):
    ax.set_yscale("log")
    ax.set_ylim(10**-19.2, 3e-3)
    ticks = [ZERO_ROW, FLOOR, 1e-14, 1e-11, 1e-8, 1e-5]
    labels = ["0", "≤10$^{-17}$", "10$^{-14}$", "10$^{-11}$", "10$^{-8}$", "10$^{-5}$"]
    ax.yaxis.set_major_locator(ticker.FixedLocator(ticks))
    ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))
    ax.yaxis.set_minor_locator(ticker.NullLocator())
    ax.axhspan(10**-19.2, 10**-17.8, color="#f3f2ef", lw=0, zorder=0)
    for y in ticks[1:]:
        ax.axhline(y, color="#e6e5e1", lw=0.5, zorder=0)


def _iterations(rows):
    n = np.array([r["prusie"]["niter"][0] for r in rows], dtype=float)
    # Integer counts: spread coincident points horizontally within ±0.3 of the count.
    return n + np.random.default_rng(11).uniform(-0.3, 0.3, n.size)


def _agreement_panel(ax, rows, field, panel, p):
    e = np.array([r[field] for r in rows], dtype=float)
    limit = np.array([r["both_nonconverged_at_limit"] for r in rows])
    zero = e == 0
    shown = np.where(zero, ZERO_ROW, np.maximum(e, FLOOR))
    jitter = np.where(zero, 10 ** np.random.default_rng(7).uniform(-0.25, 0.25, e.size), 1.0)
    shown = shown * jitter
    regular = ~limit
    if panel == "GWAS":
        ax.scatter(p[regular], shown[regular], s=6, color=CASE_POINT, alpha=0.55, linewidths=0,
                   zorder=2)
    else:
        density_scatter(ax, p[regular], shown[regular], CASE_POINT, size=1.4, light=0.8)
    if limit.any():
        ax.scatter(p[limit], shown[limit], s=7, marker="D", facecolor="white",
                   edgecolor=HIGHLIGHT, linewidth=0.6, zorder=5)
    median = float(np.median(e))
    ax.axhline(median, color=INK2, lw=0.6, ls=(0, (1, 1.5)), zorder=3)
    _error_axis(ax)
    ax.set_xscale("log")
    return int(zero.sum()), int(((e > 0) & (e < FLOOR)).sum()), median


def numerical_agreement_figure(records, output, filename="supp_agreement", x="variants"):
    style("prusie-agreement")
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COLUMN_IN, 112 * MM), sharex="col")
    fig.subplots_adjust(left=0.1, right=0.99, bottom=0.09, top=0.9, wspace=0.16, hspace=0.16)
    for col, panel in enumerate(PANELS):
        rows = primary(records, panel)
        for row, (field, label) in enumerate(
            (("max_abs_alpha_error", "Max |Δα| (component probability)"),
             ("max_abs_pip_error", "Max |ΔPIP|"))
        ):
            ax = axes[row, col]
            p = (np.array([r["p"] for r in rows], dtype=float) if x == "variants"
                 else _iterations(rows))
            zeros, floored, median = _agreement_panel(ax, rows, field, panel, p)
            text = f"median {sci(median)}\nexact 0: {zeros:,}"
            if floored:
                text += f";  below 10$^{{-17}}$: {floored:,}"
            note(ax, text, x=0.03, y=0.27)
            if col == 0:
                ax.set_ylabel(label)
        axes[0, col].set_title(TITLES[panel])
        axes[1, col].set_xlabel("Variants in region, $p$" if x == "variants"
                                else "IBSS iterations")
    if x == "variants":
        plain_log_axis(axes[1, 0].xaxis, [2000, 5000, 10000, 20000])
        axes[1, 0].set_xlim(1400, 22000)
        plain_log_axis(axes[1, 1].xaxis, [10, 100, 1000])
        axes[1, 1].set_xlim(5, 8500)
    else:
        for ax in axes[1]:
            plain_log_axis(ax.xaxis, [2, 5, 10, 20, 50, 100])
            ax.set_xlim(1.5, 130)
    for ax, letter in zip(axes.flat, "abcd", strict=True):
        panel_letter(ax, letter, x=-0.13 if ax in axes[:, 0] else -0.09, y=1.025)
    handles = [
        Line2D([], [], ls="", marker="o", ms=3, mfc=CASE_POINT, mec="none",
               label="Converged"),
        Line2D([], [], ls="", marker="D", ms=3.5, mfc="white", mec=HIGHLIGHT, mew=0.9,
               label="Neither program converged (100 iterations)"),
        Line2D([], [], color=INK2, lw=0.6, ls=(0, (1, 1.5)), label="Panel median"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.55, 1.0),
               columnspacing=2.0, handletextpad=0.4)
    save(fig, output, filename)


def fitting_by_variants_figure(records, output, filename="supp_fitting"):
    style("prusie-fitting")
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COLUMN_IN, 66 * MM))
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.16, top=0.8, wspace=0.17)
    for ax, panel in zip(axes, PANELS, strict=True):
        rows = primary(records, panel)
        p = np.array([r["p"] for r in rows], dtype=float)
        order = np.argsort(p)
        values = {}
        for backend, key in (("susieR", "R"), ("prusie", "prusie")):
            y = np.array([np.median(r[key]["fitting_times"]) for r in rows])
            values[backend] = y
            if panel == "GWAS":
                ax.scatter(p, y, s=6, color=COLORS[backend], alpha=0.45, linewidths=0, zorder=2)
            else:
                density_scatter(ax, p, y, COLORS[backend])
        for backend, y in values.items():
            median_line(ax, *binned_median(p[order], y[order], *BINS[panel]), backend)
        ratio = sum(sum(r["R"]["fitting_times"]) for r in rows) / sum(
            sum(r["prusie"]["fitting_times"]) for r in rows
        )
        note(ax, f"Cumulative susieR/prusie fitting time {ratio:.1f}×")
        ax.set_xscale("log")
        ax.set_yscale("log")
        grid(ax)
        ax.set_title(TITLES[panel])
        ax.set_xlabel("Variants in region, $p$")
    plain_log_axis(axes[0].xaxis, [2000, 5000, 10000, 20000])
    axes[0].set_xlim(1400, 22000)
    plain_log_axis(axes[0].yaxis, [0.01, 0.1, 1, 10, 100])
    axes[0].set_ylim(0.01, 800)
    plain_log_axis(axes[1].xaxis, [10, 100, 1000])
    axes[1].set_xlim(5, 8500)
    plain_log_axis(axes[1].yaxis, [0.001, 0.01, 0.1, 1, 10, 100])
    axes[1].set_ylim(2e-4, 400)
    axes[0].set_ylabel("Fitting-stage time per call (s)")
    for ax, letter in zip(axes, "ab", strict=True):
        panel_letter(ax, letter, x=-0.1 if ax is axes[0] else -0.085, y=1.025)
    handles = [
        Line2D([], [], color=DARK[b], marker="o", mfc="white", mec=DARK[b], mew=1.0, lw=1.3,
               ms=3.6, label=LABELS[b])
        for b in ("prusie", "susieR")
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, bbox_to_anchor=(0.53, 1.0),
               handlelength=2.4, columnspacing=2.4)
    save(fig, output, filename)


def peak_memory_figure(records, output, filename="supp_memory"):
    style("prusie-memory")
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COLUMN_IN, 112 * MM), sharex="col")
    fig.subplots_adjust(left=0.085, right=0.99, bottom=0.09, top=0.9, wspace=0.16, hspace=0.16)
    for col, panel in enumerate(PANELS):
        rows = primary(records, panel)
        p = np.array([r["p"] for r in rows], dtype=float)
        order = np.argsort(p)
        peak = {b: np.array([r[key]["peak_rss_diagnostic_bytes"] for r in rows], dtype=float) / GIB
                for b, key in (("susieR", "R"), ("prusie", "prusie"))}
        ax = axes[0, col]
        for backend, y in peak.items():
            if panel == "GWAS":
                ax.scatter(p, y, s=6, color=COLORS[backend], alpha=0.45, linewidths=0, zorder=2)
            else:
                density_scatter(ax, p, y, COLORS[backend])
        for backend, y in peak.items():
            median_line(ax, *binned_median(p[order], y[order], *BINS[panel]), backend)
        line_p = np.geomspace(p.min() * 0.7, p.max() * 1.4, 200)
        ax.plot(line_p, 8 * line_p**2 / GIB, color=INK2, lw=0.8, ls=(0, (4, 2)), zorder=4)
        note(ax, f"median prusie {np.median(peak['prusie']) * 1024:.0f} MiB\n"
                 f"median susieR {np.median(peak['susieR']) * 1024:.0f} MiB")
        ax.set_yscale("log")
        plain_log_axis(ax.yaxis, [0.1, 0.2, 0.5, 1, 2, 5, 10, 20])
        ax.set_ylim(0.05, 30)
        grid(ax)
        ax.set_title(TITLES[panel])

        ax = axes[1, col]
        ratio = peak["susieR"] / peak["prusie"]
        if panel == "GWAS":
            ax.scatter(p, ratio, s=6, color=CASE_POINT, alpha=0.55, linewidths=0, zorder=2)
        else:
            density_scatter(ax, p, ratio, CASE_POINT, light=0.8)
        cx, cy = binned_median(p[order], ratio[order], *BINS[panel])
        ax.plot(cx, cy, color=INK, lw=1.3, marker="o", ms=3.4, mfc="white", mec=INK, mew=1.0,
                path_effects=HALO, zorder=5)
        median = float(np.median(ratio))
        ax.axhline(median, color=INK2, lw=0.6, ls=(0, (1, 1.5)), zorder=3)
        ax.axhline(1, color=MUTED, lw=0.6, zorder=1)
        ax.text(0.99, 1.12, "equal peak", transform=ax.get_yaxis_transform(), ha="right",
                va="bottom", fontsize=6, color=MUTED)
        q1, q3 = np.quantile(ratio, [0.25, 0.75])
        note(ax, f"median {median:.2f} [IQR {q1:.2f}–{q3:.2f}]\n"
                 f"range {ratio.min():.2f}–{ratio.max():.2f}")
        ax.set_ylim(0, 9)
        ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
        grid(ax)
        ax.set_xscale("log")
        ax.set_xlabel("Variants in region, $p$")
    axes[0, 0].set_ylabel("Process peak memory (GiB)")
    axes[1, 0].set_ylabel("Paired peak ratio, susieR / prusie")
    plain_log_axis(axes[1, 0].xaxis, [2000, 5000, 10000, 20000])
    axes[1, 0].set_xlim(1400, 22000)
    plain_log_axis(axes[1, 1].xaxis, [10, 100, 1000])
    axes[1, 1].set_xlim(5, 8500)
    for ax, letter in zip(axes.flat, "abcd", strict=True):
        panel_letter(ax, letter, x=-0.12 if ax in axes[:, 0] else -0.09, y=1.025)
    handles = [
        Line2D([], [], color=DARK[b], marker="o", mfc="white", mec=DARK[b], mew=1.0, lw=1.3,
               ms=3.6, label=LABELS[b])
        for b in ("prusie", "susieR")
    ] + [
        Line2D([], [], color=INK2, lw=0.8, ls=(0, (4, 2)), label="One dense LD matrix, $8p^2$ bytes"),
        Line2D([], [], color=INK, marker="o", mfc="white", mec=INK, mew=1.0, lw=1.3, ms=3.6,
               label="Binned median ratio"),
        Line2D([], [], color=INK2, lw=0.6, ls=(0, (1, 1.5)), label="Panel median ratio"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=5, bbox_to_anchor=(0.54, 1.0),
               handlelength=2.2, columnspacing=1.6, handletextpad=0.4)
    save(fig, output, filename)
