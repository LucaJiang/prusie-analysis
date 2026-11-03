"""Shared visual style for the manuscript and supplementary figures."""
from pathlib import Path
import os

PAPER = Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR', str(PAPER / '.cache/matplotlib'))
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Categorical slots validated together (all pairs, light surface): normal-vision and CVD separation pass.
PRUSIE = '#2a78d6'
SUSIER = '#eb6834'
GWAS = '#4a3aa7'
EQTL = '#199e70'
BACKEND = {'prusie': PRUSIE, 'susieR': SUSIER}
PANEL = {'GWAS': GWAS, 'eQTL': EQTL}

INK = '#1f1f1e'
INK2 = '#52514e'
MUTED = '#8a8984'
GRID = '#e6e5e1'
SURFACE = '#ffffff'

MM = 1 / 25.4
DOUBLE_COLUMN = 178 * MM
SINGLE_COLUMN = 86 * MM


def apply():
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 7,
        'axes.labelsize': 7, 'axes.titlesize': 7.5, 'axes.titleweight': 'bold',
        'axes.titlelocation': 'left', 'axes.titlepad': 6,
        'legend.fontsize': 6.5, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
        'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
        'text.color': INK, 'axes.linewidth': 0.6,
        'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
        'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.5, 'grid.linestyle': '-',
        'axes.axisbelow': True, 'axes.facecolor': SURFACE, 'figure.facecolor': SURFACE,
        'lines.linewidth': 1.2, 'lines.solid_capstyle': 'round', 'lines.solid_joinstyle': 'round',
        'legend.frameon': False, 'legend.handlelength': 1.2, 'legend.borderaxespad': 0.3,
        'pdf.fonttype': 42, 'svg.fonttype': 'none', 'svg.hashsalt': 'prusie-manuscript',
        'savefig.dpi': 600, 'figure.dpi': 150,
    })


def panel_label(ax, letter, x=-0.16, y=1.04):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom', ha='left')


def save(fig, stem, png_dpi=300):
    """Write deterministic PDF, SVG and PNG copies to figures/."""
    dest = PAPER / '04_figures_tables/figures'
    dest.mkdir(exist_ok=True)
    fig.savefig(dest / f'{stem}.pdf', metadata={'CreationDate': None, 'ModDate': None})
    fig.savefig(dest / f'{stem}.svg', metadata={'Date': None})
    fig.savefig(dest / f'{stem}.png', dpi=png_dpi)
    plt.close(fig)
