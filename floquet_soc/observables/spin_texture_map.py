# -*- coding: utf-8 -*-
"""Generic static spin-texture maps for any Wannier/analytic model.

This is the program behind ``floquet.py spin_texture``.  It is material
agnostic: it only needs a loaded model, a spin operator (``spin_at_k``) and a
band index.

Two grid modes
--------------
full   : a plane through the BZ (``plane_kpts``); plotted with the BZ polygon
         and clipped to it (the PbNV ``spin_static`` convention).
window : a square window around a chosen k point, in 1/Angstrom; plotted
         without a BZ polygon (valley / Gamma zooms).

Both modes share the same spin expectation:
    s_b(k) = ( <u_b|S_x(k)|u_b>, <u_b|S_y(k)|u_b>, <u_b|S_z(k)|u_b> )
with u_b the eigenvector of the (wsvec-aware) H(k).  Single-band spin is
gauge dependent when bands are degenerate; ``mask_split`` sets the in-plane
components to zero (drawn as dots) where the neighbouring gap is below that
threshold.
"""
from __future__ import annotations

import os

import numpy as np

from .spin_plot import plane_kpts, plot_spin_texture_2d


def resolve_band(cfg: dict, model, band=None) -> int:
    """Resolve ``band`` = None/'vbm'/'cbm'/int into a 0-based band index."""
    nw = int(model.nw)
    ne = cfg.get("num_electrons")
    if band is None:
        band = cfg.get("spin_plot_band")
        if band is None:
            band = (ne - 1) if ne else 0
    if isinstance(band, str):
        key = band.strip().lower()
        if key in ("vbm", "top", "valence"):
            if not ne:
                raise ValueError("vbm requires num_electrons in the config")
            band = ne - 1
        elif key in ("cbm", "bottom", "conduction"):
            if not ne:
                raise ValueError("cbm requires num_electrons in the config")
            band = ne
        else:
            band = int(key)
    band = int(band)
    if not (0 <= band < nw):
        raise ValueError(f"band index {band} outside [0, {nw - 1}]")
    return band


def full_kpts(cfg: dict, n: int = 41, kmax: float = 1.65) -> np.ndarray:
    """Fractional k points on the config plane covering the BZ."""
    return np.asarray(plane_kpts(cfg.get("kplane"), int(n), float(kmax)),
                      float)


def window_kpts(lattice, center_frac, half: float, n: int = 41):
    """Square window around ``center_frac``; returns (kfrac, kcart_xy).

    ``half`` is in 1/Angstrom, ``kcart_xy`` is the plotting coordinate.
    """
    lat = np.asarray(lattice, float)
    B = 2.0 * np.pi * np.linalg.inv(lat).T
    center = np.asarray(center_frac, float).reshape(3)
    ccart = center @ B
    q = np.linspace(-float(half), float(half), int(n))
    qx, qy = np.meshgrid(q, q)
    kcart = np.stack([ccart[0] + qx.ravel(),
                      ccart[1] + qy.ravel(),
                      np.zeros(qx.size)], axis=1)
    kfrac = kcart @ np.linalg.inv(B)
    return kfrac, kcart[:, :2]


def compute_map(model, spin_at_k, band: int, kpts, mask_split: float = 0.0,
                inplane_ratio: float = 0.0, chunk: int = 64):
    """Energies and band-resolved spin for every k point.

    Returns dict(kpts, energies, spins, spins_band, split, band).
    ``spins_band`` has shape (nk,3); ``spins`` all bands (nk,nw,3).
    """
    kpts = np.atleast_2d(np.asarray(kpts, float))
    nk = len(kpts)
    nw = int(model.nw)
    energies = np.empty((nk, nw))
    spins = np.empty((nk, nw, 3))
    for i0 in range(0, nk, max(1, int(chunk))):
        sl = slice(i0, min(i0 + chunk, nk))
        for j, k in enumerate(kpts[sl]):
            Hk = model.Hk(k)
            w, U = np.linalg.eigh(Hk)
            ops = np.asarray(spin_at_k(k))
            for b in range(nw):
                u = U[:, b]
                spins[i0 + j, b] = np.real(
                    [u.conj() @ ops[c] @ u for c in range(3)])
            energies[i0 + j] = w
    band_spins = spins[:, band, :].copy()
    # neighbour gap used only to suppress gauge-dependent arrows
    gap = np.full(nk, np.inf)
    for i in range(nk):
        e = energies[i]
        others = np.delete(e, band)
        gap[i] = float(np.min(np.abs(others - e[band])))
    if mask_split and mask_split > 0:
        bad = gap < float(mask_split)
        band_spins[bad, 0] = 0.0
        band_spins[bad, 1] = 0.0
    # user rule: |S_xy|/|S_z| below the cutoff -> treat as out-of-plane
    # (Sx=Sy=0, drawn as a dot); cutoff is a parameter, not a hard-coded
    # plotting convention.
    sxy = np.hypot(band_spins[:, 0], band_spins[:, 1])
    sz = np.abs(band_spins[:, 2])
    out_of_plane = sxy <= float(inplane_ratio) * np.maximum(sz, 1e-12)
    band_spins[out_of_plane, 0] = 0.0
    band_spins[out_of_plane, 1] = 0.0
    return dict(kpts=kpts, energies=energies, spins=spins,
                spins_band=band_spins, split=gap, band=int(band),
                inplane_ratio=float(inplane_ratio),
                n_out_of_plane=int(out_of_plane.sum()))


def plot_map(kpts, spins_band, path, *, lattice=None, kplane=None,
             colorbar_label="$S_z$", figsize=(3.2, 3.0), dot_ratio=None):
    """Plot one band's spin texture (BZ-clipped when lattice is given).

    ``dot_ratio`` is the arrow-vs-dot threshold; when omitted the plot
    default (0.25) is used.  CLI callers pass the same ``inplane_ratio`` so
    the user cutoff is what decides "out of plane".
    """
    kw = {} if dot_ratio is None else dict(dot_ratio=float(dot_ratio))
    plot_spin_texture_2d(kpts, np.asarray(spins_band, float), path,
                         lattice=lattice, kplane=kplane,
                         colorbar_label=colorbar_label, figsize=figsize, **kw)
    print(f"figure: {path}")
    return path
