# -*- coding: utf-8 -*-
"""Electron spectral function A(k,E) for Floquet states (stage B).

The Floquet spectral function is built by unfolding the photon sectors onto
the original energy axis:

    A(k,E) = Σ_{α,m} w_{αm} · G_σ(E - (eps_α - m ħω))

with w_{αm} the photon-sector weight of Floquet mode α (sum_m w = 1) and G_σ
a Gaussian of width σ.  This is a *weight density* (integral over E gives the
total weight in the window); the energy convention is the same as in
``solvers.evolution`` (E = eps - m ħω).
"""
from __future__ import annotations

import numpy as np


def _accumulate(E0, weights, E_grid, sigma):
    """Gaussian-broadened accumulation of delta peaks onto E_grid."""
    E = np.asarray(E_grid, float)
    nE = len(E)
    dE = float(E[1] - E[0])
    half = int(max(4 * sigma / dE, 3))
    kern = np.exp(-0.5 * (np.arange(-half, half + 1) * dE / sigma) ** 2)
    kern /= kern.sum() * dE
    out = np.zeros(nE)
    for e, w in zip(np.atleast_1d(E0), np.atleast_1d(weights)):
        if w == 0.0:
            continue
        i0 = int(np.round((e - E[0]) / dE))
        lo, hi = i0 - half, i0 + half + 1
        if hi < 0 or lo >= nE:
            continue
        klo, khi = max(0, -lo), min(2 * half + 1, nE - lo)
        lo, hi = max(0, lo), min(nE, hi)
        out[lo:hi] += w * kern[klo:khi]
    return out


def spectral_function(eps_u, ms, wm, E_grid, sigma, hw):
    """A(k,E) from the photon decomposition of a single k point.

    Parameters
    ----------
    eps_u : (nw,) unsorted quasienergies
    ms : (nm,) photon indices
    wm : (nw,nm) photon weights
    E_grid : (nE,) energy grid
    sigma : Gaussian broadening (eV)
    hw : photon energy (eV)

    Returns (nE,) spectral weight density [1/eV].
    """
    ms = np.asarray(ms)
    wm = np.asarray(wm)
    E_grid = np.asarray(E_grid, float)
    A = np.zeros(len(E_grid))
    for im, m in enumerate(ms):
        A += _accumulate(np.asarray(eps_u) - m * hw, wm[:, im], E_grid, sigma)
    return A


def spectral_map(kpts, floquet_data, E_grid, sigma, hw):
    """A(k,E) on a k path: returns (nk,nE)."""
    return np.array([spectral_function(f["eps_u"], f["wspec"]["ms"],
                                       f["wspec"]["wm"], E_grid, sigma, hw)
                     for f in floquet_data])
