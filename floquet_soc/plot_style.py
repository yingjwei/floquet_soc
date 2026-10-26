# -*- coding: utf-8 -*-
"""PRL-style figure helpers for the command-line plots (port of phonkit style).

Style reference (not imported, kept as documentation):
    C:\\Users\\15756\\Desktop\\phonkit\\src\\phonkit\\band\\plot.py
    - _PRL_STYLE: Arial + stixsans, axes lw 1.2, tick 4 pt, labels 14 pt,
      ticks 11 pt, legend 11 pt, figure dpi 150;
    - Wong (2011) color-blind-safe palette;
    - single-column square figures (8 cm), 600 dpi, tight bounding box;
    - dashed Fermi line, dashed vertical high-symmetry lines, 12 pt 螕 labels;
    - band lines lw 0.6, alpha 0.85.

This module only affects presentation.  It contains no physics and is used by
``floquet_soc.cli`` for the ``bands`` / ``soc`` figures (v0.1.2-rc1+).
Spin-texture maps (bwr + arrows) are intentionally NOT restyled here.
"""
from __future__ import annotations

from contextlib import contextmanager

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
import numpy as np                                  # noqa: E402

# Global figure typography.  EVERY figure in the framework uses these sizes so
# the relative text scale is identical everywhere (band plots, spin textures,
# colorbars, difference plots, benchmarks, ...).
STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "mathtext.fontset": "stixsans",
    "axes.linewidth": 1.2,
    "xtick.major.size": 4,
    "xtick.major.width": 1.0,
    "ytick.major.size": 4,
    "ytick.major.width": 1.0,
    "axes.labelsize": 11,          # axis + colorbar labels
    "axes.titlesize": 10,          # kept for completeness; figures use NO titles
    "figure.titlesize": 10,
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.dpi": 150,
    "savefig.dpi": 600,
}

# One consistent type scale for every text element (fixes the old mix of
# 8/9/11/12/14/18 pt that looked uneven between different figure types).
FONTS = dict(label=11, tick=10, title=10, suptitle=10, legend=9, sym=10,
             cbar=11, cbar_tick=10,
             cbar_spin=11,     # spin-texture colorbar labels (unified scale)
             cbar_spin_tick=11)  # spin-texture colorbar ticks (~ tick size)

# Figure policy: no titles/suptitles anywhere in the framework output.

# Wong 2011 palette (same constants as phonkit + extended for multi-pair plots)
COLOR_SINGLE = "#333333"      # single band set (phonkit COLOR_FERMI)
COLOR_REF = "#0072B2"         # reference set (DFT / static)  [phonkit COLOR_DFT]
COLOR_CMP = "#D55E00"         # comparison set (Wannier / Floquet) [WANNIER]
COLOR_FERMI = "#333333"
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7",
           "#E69F00", "#56B4E9", "#F0E442", "#000000"]

FIG_SINGLE = (8.0 / 2.54, 8.0 / 2.54)        # PRL single column, square
FIG_2PANEL = (13.0 / 2.54, 5.2 / 2.54)       # 2 panels in a row
FIG_3PANEL = (16.0 / 2.54, 5.0 / 2.54)       # 3 panels in a row (wide)


@contextmanager
def apply_style():
    """Context manager applying the PRL rc defaults."""
    with plt.rc_context(STYLE):
        yield


def enforce_style():
    """Apply the global typography to the current matplotlib rcParams.

    Used by plotting modules that cannot easily wrap their body in a context
    manager (spin-texture maps, benchmark figures).  Keeps every figure's text
    scale identical.
    """
    plt.rcParams.update(STYLE)


def add_colorbar(fig, mappable, ax, label="", label_fontsize=None,
                tick_fontsize=None, **kw):
    """Create a colorbar with the global label/tick font sizes."""
    cbar = fig.colorbar(mappable, ax=ax, **kw)
    if label:
        cbar.set_label(label, fontsize=label_fontsize or FONTS["cbar"])
    cbar.ax.tick_params(labelsize=tick_fontsize or FONTS["cbar_tick"])
    return cbar


def kpath_distance(kpts, lattice=None):
    """Cumulative Cartesian k-distance in 1/Angstrom from fractional k-points.

    ``lattice`` rows are the real-space lattice vectors (Angstrom).  If it is
    None (analytic models without a lattice), a unit index axis is returned so
    the plot still works.
    """
    kpts = np.atleast_2d(np.asarray(kpts, float))
    if lattice is None:
        return np.arange(len(kpts), dtype=float)
    A = np.asarray(lattice, float)
    B = 2.0 * np.pi * np.linalg.inv(A).T          # rows b1,b2,b3
    kcart = kpts @ B
    d = np.linalg.norm(np.diff(kcart[:, :2], axis=0), axis=1)
    return np.concatenate([[0.0], np.cumsum(d)])


def sym_positions(kd, labels):
    """Positions of the high-symmetry labels along the (uniformly segmented) path."""
    n = len(kd)
    nl = len(labels) if labels else 0
    if nl <= 1:
        return np.array([kd[0]])
    nseg = (n - 1) // (nl - 1)
    pos = np.array([min(i * nseg, n - 1) for i in range(nl)], dtype=float)
    return kd[pos.astype(int)]


def add_sym_labels(ax, kd, labels, fontsize=None):
    """Dashed vertical high-symmetry lines + 螕 labels (phonkit _add_sym_labels)."""
    if fontsize is None:
        fontsize = FONTS["sym"]
    pos = sym_positions(kd, labels)
    ax.set_xticks(pos)
    ax.set_xticklabels([str(l).replace("GAMMA", "\u0393") for l in labels],
                       fontsize=fontsize)
    for p in pos:
        ax.axvline(p, color="black", ls="--", lw=0.4, alpha=0.4)
    ax.set_xlim(kd[0], kd[-1])


def add_reference_line(ax, y=0.0, lw=0.5, alpha=0.5):
    ax.axhline(y, color="gray", ls="--", lw=lw, alpha=alpha)


def save_figure(fig, path):
    """Save at 600 dpi with a tight box and close the figure."""
    fig.savefig(path, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)

