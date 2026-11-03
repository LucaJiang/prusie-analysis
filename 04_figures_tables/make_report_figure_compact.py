"""Compact variant of the offline-report figure for possible main-text use.

Writes figures/report_compact.{pdf,svg,png}. Uses the same screenshots as
make_report_figure.py; panel b omits the in-page navigation-link strip.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(PAPER / ".cache/matplotlib"))
sys.path.insert(0, str(PAPER / "scripts"))

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image

from _bench_style import DOUBLE_COLUMN_IN, INK, INK2, MM, PRUSIE, panel_letter, save, style

SNIPPETS = PAPER / "03_results/report_snippets"
FIGURES = PAPER / "04_figures_tables/figures"

SECTIONS = (
    ("What to check first", "b"),
    ("Main results and resolution", "b"),
    ("Analysis status", None),
    ("Figures", "c"),
    ("Sensitivity and LD diagnostics", None),
    ("Complete records and provenance", None),
    ("Machine-readable exports", None),
)
NAV_STRIP = (460, 800)


def schematic(ax):
    ax.set_axis_off()
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 4)
    for i, (title, link) in enumerate(SECTIONS):
        col, row = divmod(i, 4)
        x, y = col * 1.0, 3 - row
        ax.add_patch(FancyBboxPatch(
            (x + 0.01, y + 0.1), 0.96, 0.8, boxstyle="round,pad=0,rounding_size=0.03",
            facecolor="#eef4f9" if link else "#f5f5f3",
            edgecolor=PRUSIE if link else "#c9c8c3", lw=0.8 if link else 0.5,
            mutation_aspect=0.5,
        ))
        ax.text(x + 0.04, y + 0.5, f"{i + 1}  {title}", fontsize=5.8,
                fontweight="bold" if link else "normal", color=INK if link else INK2,
                va="center")
        if link:
            ax.text(x + 0.95, y + 0.5, link, fontsize=5.8, fontweight="bold",
                    color=PRUSIE, va="center", ha="right")


def snippet(ax, image):
    ax.imshow(image, interpolation="lanczos")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("#c9c8c3")
        spine.set_linewidth(0.6)


def checks_without_nav():
    image = np.asarray(Image.open(SNIPPETS / "checks_and_results.png").convert("RGB"))
    top, bottom = NAV_STRIP
    return np.concatenate([image[:top], image[bottom:]])


def main():
    style("prusie-report")
    b_img = checks_without_nav()
    c_img = np.asarray(Image.open(SNIPPETS / "locus_overview.png").convert("RGB"))
    width = DOUBLE_COLUMN_IN
    left, gap, right_w, top, bottom = 0.005, 0.035, 0.45, 0.95, 0.005
    mid_w = 1 - left - gap - right_w - 0.005
    c_in = right_w * width * c_img.shape[0] / c_img.shape[1]
    height = c_in / (top - bottom)
    fig = plt.figure(figsize=(width, height))
    b_h = mid_w * width * b_img.shape[0] / b_img.shape[1] / height
    ax_c = fig.add_axes((1 - 0.005 - right_w, bottom, right_w, top - bottom))
    ax_b = fig.add_axes((left, bottom, mid_w, b_h))
    a_bottom = bottom + b_h + 0.055
    ax_a = fig.add_axes((left, a_bottom, mid_w, top - a_bottom))

    schematic(ax_a)
    snippet(ax_b, b_img)
    snippet(ax_c, c_img)
    panel_letter(ax_a, "a", x=-0.012, y=1.02)
    panel_letter(ax_b, "b", x=-0.012, y=1.03)
    panel_letter(ax_c, "c", x=-0.04, y=1.0)
    save(fig, FIGURES, "report_compact")


if __name__ == "__main__":
    main()
