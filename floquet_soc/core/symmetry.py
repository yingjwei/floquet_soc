# -*- coding: utf-8 -*-
"""Symmetry helpers available in stage A.

Only time-reversal k pairing is implemented (needed to validate spin matrices:
for a TRS-preserving model eig[S(-k)] = -eig[S(k)] in the band basis).
Point-group / magnetic-symmetry constraints on effective SOC terms belong to
stage C (status: not implemented).
"""
from __future__ import annotations

import numpy as np


def find_kpairs(kpts, tol=1e-8):
    """Pair each k with its -k partner (mod reciprocal lattice vector).

    Returns dict index -> partner index (or -1 if no partner in the list).
    """
    K = np.mod(np.atleast_2d(np.asarray(kpts, float)), 1.0)
    nk = len(K)
    partner = -np.ones(nk, dtype=int)
    for i in range(nk):
        target = np.mod(-K[i], 1.0)
        d = np.minimum(np.abs(K - target), 1.0 - np.abs(K - target))
        dist = d.sum(axis=1)
        j = int(np.argmin(dist))
        if dist[j] < tol:
            partner[i] = j
    return partner


def time_reversal_spectra_antisymmetry(eigsets, partner):
    """Check eig(S(-k)) = -eig(S(k)) for spectra given per k.

    eigsets : (nk, ncomp, nband) sorted eigenvalues.
    Returns max deviation, only over k pairs (i<j).
    """
    E = np.asarray(eigsets, float)
    worst = 0.0
    for i, j in enumerate(partner):
        if j <= i or j < 0:
            continue
        worst = max(worst, float(np.max(np.abs(E[j] + E[i]))))
    return worst
