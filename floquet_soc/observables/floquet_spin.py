# -*- coding: utf-8 -*-
"""Floquet spin observables (stage B).

Periodic-averaged spin texture (micromotion included!):

    S_bar_alpha(k) = (1/T) ∫_0^T <u_alpha k(t)| S |u_alpha k(t)> dt

where ``u_alpha k(t) = e^{i eps_alpha t} U(t)|psi_alpha>`` is the periodic
Floquet mode.  Using only the t=0 eigenvector of U(T) is generally wrong; the
difference between ``S_bar`` and the t=0 spin ``S_t0`` is reported so the
micromotion contribution is explicit.

Spin-resolved spectral density (signed and non-negative versions kept apart):

    A^{S_i}(k,E)   = Σ_{α,m} w_{αm} * s_{αm,i}         * δ(E - eps_α + mħω)
    A^{|S_i|}(k,E) = Σ_{α,m} w_{αm} * |s_{αm,i}|       * δ(E - eps_α + mħω)

where ``s_{αm,i} = <u^{(m)}_α|S_i|u^{(m)}_α> / w_{αm}`` is the spin
polarization of photon component m and ``w_{αm}=||u^{(m)}||²``.  The signed
quantity is the physical spin spectral weight (may be negative); the absolute
one is a non-negative intensity.  They must not be mixed.
"""
from __future__ import annotations

import numpy as np


def _spin_expectation(u, S):
    """<u|S|u> for a matrix of columns u (nw, nstate) -> (nstate,) real."""
    return np.real(np.einsum("ia,ij,ja->a", u.conj(), S, u))


def floquet_spin_observables(model, drive, k, spin_ops, n_step=150, m_sp=16):
    """Stage-B spin observables at one k point.

    Returns dict with:
        eps, eps_u, Ebar, S_bar (nw,3), S_t0 (nw,3),
        ms, wm (nw,nm), spin_weights (nm,nw,3) signed,
        micromotion_norm (nw,), unitarity_error, hw.
    """
    from ..solvers.evolution import floquet_modes
    fl = floquet_modes(model, drive, k, n_step=n_step, want_wspec=True,
                       m_sp=m_sp)
    mic = fl["micromotion"]
    Ut, t, T, dt = mic["Ut"], mic["t"], mic["T"], mic["dt"]
    eps_u = fl["eps_u"]
    NS = Ut.shape[0]
    nw = model.nw

    # u(t) = e^{i eps t} U(t)|psi>
    ph = np.exp(1j * np.outer(t, eps_u))              # (NS,nw)
    u = Ut * ph[:, None, :]                           # (NS,nw,nstate)
    S_bar = np.zeros((nw, 3))
    for c in range(3):
        acc = np.zeros(nw)
        for n in range(NS):
            acc += _spin_expectation(u[n], spin_ops[c])
        S_bar[:, c] = acc * dt / T
    norms = np.real(np.einsum("nia,nia->na", u.conj(), u))
    micromotion_norm = norms.max(axis=0) - norms.min(axis=0)

    # t=0 spin (t=0 Floquet mode = eigenvector of U(T))
    evec = fl["evec"]
    S_t0 = np.zeros((nw, 3))
    for c in range(3):
        S_t0[:, c] = _spin_expectation(evec, spin_ops[c])

    # photon components u^{(m)} and their spin polarization
    ms = fl["wspec"]["ms"]
    nm = len(ms)
    spin_weights = np.zeros((nm, nw, 3))
    wm = fl["wspec"]["wm"]
    for im, m in enumerate(ms):
        phm = np.exp(1j * np.outer(t, eps_u - m * drive.hw))
        um = np.einsum("na,nia->ia", phm, Ut) / NS
        for c in range(3):
            spin_weights[im, :, c] = _spin_expectation(um, spin_ops[c])

    # The Schur decomposition returns eigenvalues in an arbitrary order; sort
    # every per-state quantity by quasienergy so that state labels are stable
    # and comparable with the static band ordering (sorted eigenvalues).
    order = np.argsort(eps_u)
    return dict(eps=fl["eps"], eps_u=eps_u[order], Ebar=fl["Ebar"][order],
                S_bar=S_bar[order], S_t0=S_t0[order],
                ms=ms, wm=wm[order], w_sum=fl["wspec"]["w_sum"][order],
                spin_weights=spin_weights[:, order, :],
                micromotion_norm=micromotion_norm[order],
                evec=fl["evec"][:, order],
                unitarity_error=fl["unitarity_error"], hw=drive.hw,
                w_main=fl["wspec"]["w_main"][order])


def spin_spectral_density(obs, E_grid, sigma):
    """Signed / non-negative spin-resolved spectral densities (see docstring).

    Energy convention: E_{αm} = eps_α - m ħω.
    """
    E = np.asarray(E_grid, float)
    nE = len(E)
    dE = float(E[1] - E[0])
    half = int(max(4 * sigma / dE, 3))
    kern = np.exp(-0.5 * (np.arange(-half, half + 1) * dE / sigma) ** 2)
    kern /= kern.sum() * dE
    A_s = np.zeros((3, nE))
    A_a = np.zeros((3, nE))
    eps, hw = obs["eps_u"], obs["hw"]
    for im, m in enumerate(obs["ms"]):
        E0 = eps - m * hw
        for a in range(len(eps)):
            i0 = int(np.round((E0[a] - E[0]) / dE))
            lo, hi = i0 - half, i0 + half + 1
            if hi < 0 or lo >= nE:
                continue
            klo, khi = max(0, -lo), min(2 * half + 1, nE - lo)
            lo, hi = max(0, lo), min(nE, hi)
            for c in range(3):
                w = obs["spin_weights"][im, a, c]
                A_s[c, lo:hi] += w * kern[klo:khi]
                A_a[c, lo:hi] += abs(w) * kern[klo:khi]
    return dict(A_signed=A_s, A_abs=A_a, E_grid=E, sigma=sigma)


def periodic_spin_texture(model, drive, spin_ops, kpts, n_step=150,
                          m_sp=16):
    """S_bar_alpha(k) along a list of k points (+ static comparison).

    Returns dict:
        kpts, eps (nk,nw), S_bar (nk,nw,3), S_t0 (nk,nw,3),
        S_static (nk,nw,3)  [spin texture of the undriven H(k)],
        S_drive (delta drive strength used for the static reference).
    """
    kpts = np.atleast_2d(np.asarray(kpts, float))
    out = None
    from ..observables.spin_texture import spin_texture
    for ik, k in enumerate(kpts):
        obs = floquet_spin_observables(model, drive, k, spin_ops,
                                       n_step=n_step, m_sp=m_sp)
        stat = spin_texture(model.Hk(k), spin_ops)
        if out is None:
            nw = model.nw
            out = dict(kpts=kpts, eps=np.empty((len(kpts), nw)),
                       S_bar=np.empty((len(kpts), nw, 3)),
                       S_t0=np.empty((len(kpts), nw, 3)),
                       S_static=np.empty((len(kpts), nw, 3)))
        out["eps"][ik] = np.sort(obs["eps"])
        out["S_bar"][ik] = obs["S_bar"]
        out["S_t0"][ik] = obs["S_t0"]
        out["S_static"][ik] = stat["spins"]
    return out
