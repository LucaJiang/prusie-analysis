"""Distribution of per-case speed-ups s_i (total execution time), GWAS and eQTL.

Writes figures/supp_speedup_dist.{pdf,svg,png}.
"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figstyle as fs
import matplotlib.pyplot as plt

PRIMARY = fs.PAPER / '03_results/benchmark_primary.tsv'
XLIM = (8, 160)


def panel(ax, s, color, title):
    x = np.log10(s)
    bins = np.linspace(np.log10(XLIM[0]), np.log10(XLIM[1]), 41)
    ax.hist(x, bins=bins, density=True, color=color, alpha=0.35, edgecolor=fs.SURFACE, linewidth=0.4)
    grid = np.linspace(bins[0], bins[-1], 400)
    ax.plot(grid, gaussian_kde(x)(grid), color=color, lw=1.2)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.12)
    med = np.median(s)
    ax.axvline(np.log10(med), color=fs.INK2, lw=0.8, ls='--')
    ax.text(np.log10(med), ax.get_ylim()[1] * 0.97, f' median {med:.1f}×', color=fs.INK, va='top', ha='left')
    ticks = [10, 20, 50, 100]
    ax.set_xticks(np.log10(ticks), [f'{t}×' for t in ticks])
    ax.set_xlim(bins[0], bins[-1])
    ax.set_xlabel('Per-case speed-up $s_i$ (susieR / prusie), log scale')
    ax.set_ylabel('Density (per log$_{10}$ unit)')
    ax.set_title(f'{title} ({len(s):,} cases)')
    ax.grid(axis='x', visible=False)


def main():
    fs.apply()
    d = pd.read_csv(PRIMARY, sep='\t')
    fig, axes = plt.subplots(1, 2, figsize=(fs.DOUBLE_COLUMN, 62 * fs.MM), constrained_layout=True)
    for ax, letter, (name, label) in zip(axes, 'ab', [('GWAS', 'GWAS loci'), ('eQTL', 'eQTL genes')]):
        panel(ax, d.loc[d.panel == name, 'speedup'].to_numpy(), fs.PANEL[name], label)
        fs.panel_label(ax, letter)
    fs.save(fig, 'supp_speedup_dist')


if __name__ == '__main__':
    main()
