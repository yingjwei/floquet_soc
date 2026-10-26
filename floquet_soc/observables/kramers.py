# -*- coding: utf-8 -*-
"""Kramers-degeneracy diagnostics (PT-symmetric non-magnetic systems).

Band index convention: 0-based numpy columns.  Kramers pairs are
(0,1), (2,3), ... i.e. 1-based (1,2), (3,4), ... (E_{2n-1}=E_{2n}).
"""
from __future__ import annotations

import numpy as np


def kramers_delta(energies):
    """delta(k) = |E_{2n}(k)-E_{2n+1}(k)| (0-based); shape (nk, nb//2)."""
    E = np.atleast_2d(np.asarray(energies, float))
    nb = E.shape[1] - E.shape[1] % 2
    return np.abs(E[:, 0:nb:2] - E[:, 1:nb:2])


def inter_doublet_gaps(energies):
    """Gap between consecutive Kramers doublets: E_{2n+2}-E_{2n+1} (0-based)."""
    E = np.atleast_2d(np.asarray(energies, float))
    nb = E.shape[1] - E.shape[1] % 2
    if nb < 4:
        return np.empty((E.shape[0], 0))
    return E[:, 2:nb:2] - E[:, 1:nb - 1:2]


def kramers_summary(energies, kpts=None):
    """Statistics of delta(k) and the worst k point (fractional)."""
    d = kramers_delta(energies)
    k = None if kpts is None else np.atleast_2d(np.asarray(kpts, float))
    worst_idx = int(np.argmax(d.max(axis=1)))
    return dict(
        nk=int(d.shape[0]), n_pairs=int(d.shape[1]),
        delta_max=float(d.max()),
        delta_median=float(np.median(d)),
        delta_p95=float(np.percentile(d, 95)),
        worst_k=(None if k is None else k[worst_idx].tolist()))
