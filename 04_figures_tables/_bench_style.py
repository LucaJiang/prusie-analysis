"""Shared style and drawing helpers for the manuscript benchmark figures."""

from __future__ import annotations

import gzip
import json
import os
import subprocess
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PAPER / ".cache/matplotlib"))

import matplotlib

matplotlib.use("Agg")
import numpy as np
from matplotlib import colors as mcolors
from matplotlib import font_manager, ticker
from matplotlib import patheffects as pe
from matplotlib import pyplot as plt

# Okabe-Ito blue/vermillion; validated (light surface) for CVD and normal-vision separation.
PRUSIE = "#0072B2"
SUSIE = "#D55E00"
COLORS = {"prusie": PRUSIE, "susieR": SUSIE}
DARK = {"prusie": "#004a75", "susieR": "#8f3d00"}
LABELS = {"prusie": "prusie 0.3.0", "susieR": "susieR 0.16.6"}
HIGHLIGHT = "#CC79A7"

INK = "#1f1f1e"
INK2 = "#52514e"
MUTED = "#8a8984"
GRID = "#e6e5e1"
NEUTRAL = INK2

MM = 1.0 / 25.4
DOUBLE_COLUMN_IN = 180.0 * MM
SINGLE_COLUMN_IN = 86.0 * MM
GIB = 2.0**30

RECORDS = PAPER / "03_results/benchmark_records.jsonl.gz"
FIGURES = PAPER / "04_figures_tables/figures"


def _sans_font():
    try:
        path = subprocess.run(
            ["kpsewhich", "texgyreheros-regular.otf"], capture_output=True, text=True, check=False
        ).stdout.strip()
    except OSError:
        path = ""
    if path:
        for face in ("regular", "bold", "italic"):
            candidate = path.replace("regular", face)
            if Path(candidate).exists():
                font_manager.fontManager.addfont(candidate)
        return "TeX Gyre Heros"
    return "DejaVu Sans"


def style(salt="prusie-manuscript"):
    plt.rcParams.update(
        {
            "font.family": _sans_font(),
            "mathtext.fontset": "custom",
            "mathtext.it": "TeX Gyre Heros:italic",
            "mathtext.rm": "TeX Gyre Heros",
            "font.size": 7,
            "axes.labelsize": 7,
            "axes.titlesize": 7.5,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "axes.titlepad": 5,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 7,
            "text.color": INK,
            "axes.labelcolor": INK,
            "axes.edgecolor": INK2,
            "xtick.color": INK2,
            "ytick.color": INK2,
            "xtick.labelcolor": INK,
            "ytick.labelcolor": INK,
            "axes.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.minor.width": 0.4,
            "ytick.minor.width": 0.4,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "xtick.minor.size": 1.4,
            "ytick.minor.size": 1.4,
            "legend.frameon": False,
            "lines.solid_capstyle": "round",
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": salt,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def load_records(path=RECORDS):
    with gzip.open(Path(path), "rt") as stream:
        return [json.loads(line) for line in stream]


def primary(records, panel):
    return [r for r in records if r["panel"] == panel and r["kind"] == "primary"]


def plain_log_axis(axis, ticks):
    axis.set_major_locator(ticker.FixedLocator(ticks))
    axis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:,.10g}"))
    axis.set_minor_locator(ticker.LogLocator(base=10, subs=np.arange(2, 10)))
    axis.set_minor_formatter(ticker.NullFormatter())


def grid(ax, axis="y"):
    ax.grid(axis=axis, which="major", color=GRID, lw=0.5, zorder=0)
    ax.set_axisbelow(True)


def tint(colour, amount):
    """Mix `colour` towards white; amount 0 is the colour, 1 is white."""
    rgb = np.array(mcolors.to_rgb(colour))
    return tuple(rgb + (1 - rgb) * amount)


def density_scatter(ax, x, y, colour, size=1.4, light=0.82, bins=(70, 55)):
    """Scatter coloured by local log-space density: sparse points pale, dense points saturated."""
    lx, ly = np.log10(x), np.log10(y)
    h, xe, ye = np.histogram2d(lx, ly, bins=bins)
    ix = np.clip(np.digitize(lx, xe) - 1, 0, h.shape[0] - 1)
    iy = np.clip(np.digitize(ly, ye) - 1, 0, h.shape[1] - 1)
    d = np.log1p(h[ix, iy]) / np.log1p(h.max())
    order = np.argsort(d, kind="stable")
    base = np.array(mcolors.to_rgb(colour))
    rgb = base + (1 - base) * (light * (1 - d[order]))[:, None]
    ax.scatter(x[order], y[order], s=size, c=rgb, linewidths=0, rasterized=True, zorder=2)


def binned_median(p, y, n_bins, min_count):
    edges = np.geomspace(p.min(), p.max() * (1 + 1e-12), n_bins + 1)
    index = np.digitize(p, edges) - 1
    centres, medians = [], []
    for i in range(n_bins):
        mask = index == i
        if mask.sum() >= min_count:
            centres.append(np.median(p[mask]))
            medians.append(np.median(y[mask]))
    return np.array(centres), np.array(medians)


HALO = [pe.Stroke(linewidth=3.0, foreground="white"), pe.Normal()]


def median_line(ax, x, y, backend, **kwargs):
    ax.plot(
        x,
        y,
        color=DARK[backend],
        lw=1.3,
        marker="o",
        ms=3.4,
        mfc="white",
        mec=DARK[backend],
        mew=1.0,
        path_effects=HALO,
        zorder=5,
        **kwargs,
    )


def note(ax, text, x=0.03, y=0.97, **kwargs):
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.5,
        color=INK,
        path_effects=[pe.Stroke(linewidth=2.5, foreground="white"), pe.Normal()],
        zorder=7,
        **kwargs,
    )


def panel_letter(ax, letter, x=-0.14, y=1.04):
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=9,
        fontweight="bold",
        va="bottom",
        ha="left",
        color=INK,
    )


def save(fig, output, name, *, png_dpi=300, tight=False):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    extra = {"bbox_inches": "tight", "pad_inches": 0.02} if tight else {}
    fig.savefig(output / f"{name}.pdf", metadata={"CreationDate": None, "ModDate": None}, dpi=600, **extra)
    fig.savefig(output / f"{name}.svg", metadata={"Date": None}, dpi=600, **extra)
    fig.savefig(output / f"{name}.png", dpi=png_dpi, **extra)
    plt.close(fig)
