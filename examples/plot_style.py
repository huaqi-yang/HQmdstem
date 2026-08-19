#!/usr/bin/env python3
"""
Shared plotting style module — combines nature-skills palettes/helpers
with GPUMDkit formatting conventions (4-spine, inward ticks, compact sizing).

Usage:
    from plot_style import *
    apply_style()                # call once before any figure
    ...
    finalize_figure(fig, "figures/output", formats=["png","svg"], dpi=300)
"""

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# ═══════════════════════════════════════════════════════════════
# PALETTES (from nature-skills)
# ═══════════════════════════════════════════════════════════════

PALETTE = {
    "blue_main":      "#0F4D92",
    "blue_secondary": "#3775BA",
    "green_1": "#DDF3DE",
    "green_2": "#AADCA9",
    "green_3": "#8BCF8B",
    "red_1":   "#F6CFCB",
    "red_2":   "#E9A6A1",
    "red_strong": "#B64342",
    "neutral_light": "#CFCECE",
    "neutral_mid":   "#767676",
    "neutral_dark":  "#4D4D4D",
    "neutral_black": "#272727",
    "gold":   "#FFD700",
    "teal":   "#42949E",
    "violet": "#9A4D8E",
    "magenta":"#EA84DD",
}

DEFAULT_COLORS = [
    PALETTE["blue_main"],
    PALETTE["green_3"],
    PALETTE["red_strong"],
    PALETTE["teal"],
    PALETTE["violet"],
    PALETTE["neutral_light"],
]

PALETTE_NATURE_MATERIAL = {
    "aqua": "#77D7D1",
    "teal": "#33B5A5",
    "lilac": "#B9A7E8",
    "violet": "#7C6CCF",
    "callout_red": "#E53935",
    "neutral": "#D9D9D9",
}

# GPUMDkit scientific palette (PRL-style)
PRL_COLORS = ['#457B9D', '#D62828', '#2A9D8F', '#E9C46A']

# ═══════════════════════════════════════════════════════════════
# STYLE APPLICATION (nature-skills font + GPUMDkit formatting)
# ═══════════════════════════════════════════════════════════════

def apply_style(font_size=10, axes_linewidth=1.0, tick_length=3.5,
                tick_width=1.0, use_tex=False):
    """
    Apply combined nature-skills + GPUMDkit publication style.

    Font:  Arial (nature-skills + GPUMDkit agree)
    Spines: all 4 visible, linewidth ~1.0 (GPUMDkit convention)
    Ticks:  inward, top+right visible (GPUMDkit convention)
    Legend: frameon=False (both agree)
    SVG:    editable text (nature-skills)
    """
    # ── MANDATORY: editable SVG text + Arial font ──
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
    plt.rcParams['svg.fonttype'] = 'none'

    # ── GPUMDkit-compatible sizing ──
    plt.rcParams['font.size'] = font_size
    plt.rcParams['axes.labelsize'] = font_size + 1
    plt.rcParams['axes.titlesize'] = font_size + 2
    plt.rcParams['xtick.labelsize'] = font_size
    plt.rcParams['ytick.labelsize'] = font_size
    plt.rcParams['legend.fontsize'] = font_size - 1

    # ── Spine & axes (GPUMDkit: all 4 visible, inward ticks) ──
    plt.rcParams['axes.spines.right'] = True
    plt.rcParams['axes.spines.top'] = True
    plt.rcParams['axes.linewidth'] = axes_linewidth

    # ── Tick direction inward (GPUMDkit) ──
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['xtick.major.size'] = tick_length
    plt.rcParams['ytick.major.size'] = tick_length
    plt.rcParams['xtick.major.width'] = tick_width
    plt.rcParams['ytick.major.width'] = tick_width
    plt.rcParams['xtick.top'] = True
    plt.rcParams['ytick.right'] = True

    # ── Minor ticks ──
    plt.rcParams['xtick.minor.visible'] = True
    plt.rcParams['ytick.minor.visible'] = True
    plt.rcParams['xtick.minor.size'] = tick_length * 0.6
    plt.rcParams['ytick.minor.size'] = tick_length * 0.6

    # ── Legend ──
    plt.rcParams['legend.frameon'] = False

    # ── LaTeX ──
    if use_tex:
        plt.rcParams['text.usetex'] = True

    print(f"Style applied: Arial {font_size}pt, "
          f"axes_lw={axes_linewidth}, ticks=inward")


# ═══════════════════════════════════════════════════════════════
# HELPER FUNCTIONS (from nature-skills API)
# ═══════════════════════════════════════════════════════════════

def is_dark(hex_color, threshold=128):
    """Return True if hex color is dark (use white text on it)."""
    c = hex_color.lstrip('#')
    r, g, b = int(c[0:2], 16), int(c[2:4], 16), int(c[4:6], 16)
    return (0.299*r + 0.587*g + 0.114*b) < threshold


def add_panel_label(ax, label, x=-0.08, y=1.04, fontsize=12,
                    color='black', fontweight='bold'):
    """Place a Nature-style bold lowercase panel label (a, b, c)."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight=fontweight,
            color=color, ha='left', va='bottom')


def finalize_figure(fig, out_path, formats=None, dpi=300,
                    pad=1.5, bbox_inches=None, close=True):
    """
    Apply tight_layout and save figure. Creates parent directories.

    Parameters
    ----------
    out_path : str   — path without extension, or with extension
    formats  : list  — e.g. ['png', 'svg', 'pdf']. If None, uses extension of out_path.
    dpi      : int   — 300 standard for publication
    pad      : float — tight_layout pad
    """
    from pathlib import Path
    fig.tight_layout(pad=pad)
    base = Path(out_path)
    os.makedirs(base.parent, exist_ok=True)
    if formats is None:
        formats = [base.suffix.lstrip('.') or 'png']
        base = base.with_suffix('')
    saved = []
    for fmt in formats:
        p = str(base) + f'.{fmt}'
        kw = {}
        if bbox_inches is not None:
            kw['bbox_inches'] = bbox_inches
        fig.savefig(p, dpi=dpi, **kw)
        saved.append(p)
    if close:
        plt.close(fig)
    return saved


def set_axes_properties(ax_list, linewidth=1.5):
    """
    Apply GPUMDkit-style formatting to specific axes:
    inward ticks, all spines visible, top/right tick labels on.
    Call this per-axis when rcParams defaults need per-axis overrides.
    """
    for ax in (ax_list if isinstance(ax_list, (list, tuple))
               else [ax_list]):
        ax.tick_params(which='both', direction='in',
                       top=True, right=True)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(linewidth)


# ═══════════════════════════════════════════════════════════════
# CONVENIENCE: get axis with GPUMDkit defaults
# ═══════════════════════════════════════════════════════════════

def subplots(*args, **kwargs):
    """
    Drop-in replacement for plt.subplots with sensible GPUMDkit defaults.
    Passes through all args but sets dpi=100 if not specified.
    """
    if 'dpi' not in kwargs:
        kwargs['dpi'] = 100
    return plt.subplots(*args, **kwargs)


# ═══════════════════════════════════════════════════════════════
# TYPICAL FIGURE SIZE PRESETS (GPUMDkit conventions)
# ═══════════════════════════════════════════════════════════════

FIG_SINGLE      = (4.5, 3.5)    # standard single panel
FIG_SINGLE_WIDE = (6.0, 4.0)    # wider single panel
FIG_DUAL_ROW    = (14.0, 8.0)   # 2 rows stacked
FIG_DUAL_COL    = (10.0, 4.0)   # 2 columns side-by-side
FIG_TRIPLE_COL  = (14.0, 4.5)   # 3 columns side-by-side
FIG_GRID_2X2    = (12.0, 8.0)   # 2x2 grid
FIG_GRID_2X3    = (14.0, 8.0)   # 2x3 grid (thermo style)
FIG_GRID_3X2    = (14.0, 10.0)  # 3x2 grid
FIG_LARGE       = (10.0, 6.0)   # large single panel
