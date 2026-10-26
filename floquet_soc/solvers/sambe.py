# -*- coding: utf-8 -*-
"""Sambe Floquet Hamiltonian (stage B) — independent cross-check of the
direct time-evolution solver.

Builds the extended-space Hamiltonian

    H_Sambe[m,n] = H_{m-n} + m hbar omega delta_{mn}

with the Fourier components H_p of H(k,t) obtained from
``evolution.hkt_components`` (Wannier: Bessel factorization; analytic:
numerical Fourier transform), then diagonalizes H_Sambe.  Its eigenvalues
folded into the first Floquet zone are the quasienergies (each physical value
appears once per photon replica).

This implementation deliberately shares only the H(k,t) definition with
``evolution.py`` and uses a completely different algorithm (extended-space
diagonalization instead of time-ordered evolution), so agreement between the
two is a meaningful numerical check.
"""
from __future__ import annotations

import numpy as np

from .evolution import hkt_components, floquet_modes, fold_quasienergy


def fourier_components(model, drive, k, n_harm=3, n_step=400):
    """Fourier components H_p, p = -n_harm..n_harm, for the Sambe matrix.

    ``n_step`` is accepted for API compatibility; the components are computed
    exactly (Wannier Bessel) or on a dense internal grid (analytic).
    """
    comps = hkt_components(model, drive, k, pmax=n_harm)
    Hp = {p: comps[p] for p in range(-n_harm, n_harm + 1) if p in comps}
    herm = max(float(np.max(np.abs(h - h.conj().T))) for h in Hp.values())
    Hp["herm_max"] = herm
    return Hp


def sambe_matrix(Hp, hw, n_harm):
    """Extended-space Sambe Hamiltonian, dimension nw*(2N+1)."""
    nw = Hp[0].shape[0]
    nh = 2 * n_harm + 1
    dim = nw * nh
    M = np.zeros((dim, dim), dtype=complex)
    for a in range(nh):
        m = a - n_harm
        for b in range(nh):
            n = b - n_harm
            p = m - n
            if p in Hp:
                M[a * nw:(a + 1) * nw, b * nw:(b + 1) * nw] = Hp[p]
        M[a * nw:(a + 1) * nw, a * nw:(a + 1) * nw] += m * hw * np.eye(nw)
    return M


def sambe_modes(model, drive, k, n_harm=3, n_step=400):
    """Diagonalize the Sambe Hamiltonian.

    Returns eigenvalues (unfolded), quasienergies folded to [-hw/2,hw/2),
    and the weight of each eigenvector in the m=0 photon sector.
    """
    Hp = fourier_components(model, drive, k, n_harm=n_harm, n_step=n_step)
    M = sambe_matrix(Hp, drive.hw, n_harm)
    ev, vec = np.linalg.eigh(0.5 * (M + M.conj().T))
    nw = model.nw
    folded = fold_quasienergy(ev, drive.hw)
    m0 = np.zeros(len(ev))
    nh = 2 * n_harm + 1
    c = n_harm * nw                       # m = 0 block offset
    for j in range(len(ev)):
        m0[j] = float(np.sum(np.abs(vec[c:c + nw, j]) ** 2))
    return dict(eigenvalues=ev, quasienergies=folded, m0_weights=m0,
                matrix=M, dim=len(ev), n_harm=n_harm,
                herm_max=Hp.get("herm_max", 0.0))


def sambe_vs_evolution(model, drive, k, n_harm=3, n_step_sambe=400,
                       n_step_evo=300, pmax_evo=None):
    """Compare Sambe quasienergies with direct-evolution quasienergies.

    For every evolution quasienergy eps_j, find the nearest Sambe eigenvalue
    modulo hbar*omega; return the maximum mismatch (eV).
    """
    evo = np.sort(floquet_modes(model, drive, k, n_step=n_step_evo,
                                pmax=pmax_evo)["eps"])
    sam = sambe_modes(model, drive, k, n_harm=n_harm, n_step=n_step_sambe)
    lam = sam["eigenvalues"]
    worst = 0.0
    detail = []
    for e in evo:
        d = np.abs(e - (lam + np.round((e - lam) / drive.hw) * drive.hw))
        j = int(np.argmin(d))
        worst = max(worst, float(d[j]))
        detail.append((float(e), float(lam[j]), float(d[j])))
    return dict(max_mismatch=worst, detail=detail, sambe=sam,
                evolution_eps=evo)
