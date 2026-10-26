# -*- coding: utf-8 -*-
"""Direct time-evolution Floquet solver (stage B).

Ported from the numerically validated PbNV engine
(``PbNV/floquet_dressed.py``), generalized to any model exposing a time
dependent Hamiltonian through ``drives.peierls.hkt``:

    H(k,t)  ->  U(k,T)  ->  Floquet modes  u_alpha(t) = e^{i eps_alpha t} U(t)|psi>

Conventions (do not mix):
    U(T)|psi> = e^{-i eps T}|psi>
    u^{(m)} = (1/T) int e^{-i m omega t} u(t) dt  ->  u(t)=sum_m e^{i m omega t} u^{(m)}
    physical state: psi(t) = e^{-i eps t} u(t) = sum_m e^{-i(eps - m omega)t} u^{(m)}
    => photon-sector energy  E_{alpha m} = eps_alpha - m hbar omega
    => mean energy            Ebar_alpha = eps_alpha - omega sum_m m w_{alpha m}

Numerical scheme (unchanged from the validated engine): midpoint
exponentialization, midpoint quadrature for the time-averaged energy, Schur
decomposition of U(T).

H(k,t) is represented by its first Fourier components H_p:
  * Wannier models: Bessel factorization of the Peierls phase (same as the
    validated PbNV engine, generalized to arbitrary polarization);
  * analytic continuum models: numerical Fourier transform of H(k+A(t))
    (minimal coupling) with a dense midpoint grid;
  * ``pmax`` truncates |p| and is recorded in the diagnostics.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import schur

from ..drives.peierls import hkt, peierls_fourier_components


def hkt_components(model, drive, k, pmax=None):
    """Fourier components H_p of H(k,t), p = -pmax..pmax."""
    if pmax is None:
        pmax = 3 if hasattr(model, "HR") else 8
    if hasattr(model, "HR"):
        return peierls_fourier_components(model, drive, k, pmax=pmax)
    # analytic continuum model: H(k + A(t)) with minimal coupling
    n_samp = max(4 * (2 * pmax + 1), 400)
    T = drive.period
    t = (np.arange(n_samp) + 0.5) * T / n_samp
    H = np.empty((n_samp, model.nw, model.nw), dtype=complex)
    for n, tn in enumerate(t):
        Ht = hkt(model, drive, k, float(tn))
        H[n] = 0.5 * (Ht + Ht.conj().T)
    Hp = {}
    for p in range(-pmax, pmax + 1):
        Hp[p] = np.tensordot(np.exp(-1j * p * drive.omega * t), H,
                             axes=(0, 0)) / n_samp
    return Hp


def sample_hkt(model, drive, k, n_step, pmax=None):
    """Hamiltonian at midpoint times t_n=(n+1/2)T/N from its Fourier components.

    Returns (H_samples (N,nw,nw), t_mid (N,), herm_max, T, dt, pmax).
    """
    comps = hkt_components(model, drive, k, pmax=pmax)
    pmax = max(abs(p) for p in comps)
    T = drive.period
    dt = T / n_step
    t = (np.arange(n_step) + 0.5) * dt
    H = np.zeros((n_step, model.nw, model.nw), dtype=complex)
    for p, Hp in comps.items():
        H += Hp[None, :, :] * np.exp(1j * p * drive.omega * t)[:, None, None]
    herm = float(np.max(np.abs(H - H.conj().transpose(0, 2, 1))))
    return H, t, herm, T, dt, pmax


def _eig_unitary(U):
    """Eigen-decompose a unitary matrix via Schur (stable); fallback to eig."""
    Tw, Z = schur(U, output="complex")
    off = float(np.max(np.abs(Tw - np.diag(np.diag(Tw)))))
    if not np.isfinite(off) or off > 1e-8:
        ev, evec = np.linalg.eig(U)
        return ev, evec, off
    return np.diag(Tw), Z, off


def fold_quasienergy(eps, hw):
    """Fold quasienergies into the first Floquet zone ``[-hw/2, hw/2)``.

    Quasienergies are defined modulo ``hw``; this is the single source of truth
    for that fold, shared by ``evolution`` and ``sambe`` so the two solvers use
    the same interval (previously one used ``(-hw/2, hw/2]`` and the other
    ``[-hw/2, hw/2)``).
    """
    return (np.asarray(eps) + hw / 2) % hw - hw / 2


def floquet_modes(model, drive, k, n_step=150, want_wspec=False, m_sp=16,
                  pmax=None):
    """Floquet modes at one k point (conventions in the module docstring).

    Returns dict with eps (sorted/folded), eps_u (unsorted), Ebar (unsorted),
    evec, U, unitarity_error, herm_max, schur_offdiag, and -- if
    ``want_wspec`` -- wspec (ms, wm, w_sum, w0, w_main, m_dom) and micromotion
    (Ut, t, T, dt) so that spin averages can include micromotion.
    """
    H, t, herm, T, dt, pmax = sample_hkt(model, drive, k, n_step, pmax=pmax)
    nw = model.nw
    U = np.eye(nw, dtype=complex)
    V = np.eye(nw, dtype=complex)
    A = np.zeros((nw, nw), dtype=complex)
    Ut_hist = np.empty((n_step, nw, nw), dtype=complex) if want_wspec else None

    for n in range(n_step):
        Ht = H[n]
        w, Q = np.linalg.eigh(Ht)
        Qh = Q.conj().T
        Uhalf = (Q * np.exp(-0.5j * w * dt)) @ Qh
        Vmid = Uhalf @ V
        A += (Vmid.conj().T @ (Ht @ Vmid)) * dt
        Ustep = (Q * np.exp(-1j * w * dt)) @ Qh
        U = Ustep @ U
        V = Ustep @ V
        if want_wspec:
            Ut_hist[n] = Vmid

    ev, evec, schur_off = _eig_unitary(U)
    eps_u = fold_quasienergy((-np.angle(ev) / T).real, drive.hw)
    Ebar = np.real(np.einsum("ia,ij,ja->a", evec.conj(), A / T, evec))
    eps = np.sort(eps_u)

    out = dict(eps=eps, eps_u=eps_u, Ebar=Ebar, evec=evec, U=U,
               unitarity_error=float(np.max(np.abs(U.conj().T @ U
                                                 - np.eye(nw)))),
               herm_max=herm, schur_offdiag=schur_off,
               k=np.asarray(k, float), hw=drive.hw,
               pmax=int(pmax),
               model_kind=getattr(model, "kind", "unknown"),
               drive=drive.as_dict())

    if want_wspec:
        phi = 2 * np.pi * (np.arange(n_step) + 0.5) / n_step
        Ut = np.einsum("nij,ja->nia", Ut_hist, evec)       # (N,nw,nw)
        tmid = phi * (T / (2 * np.pi))
        ph0 = np.exp(1j * np.outer(tmid, eps_u))
        ms = np.arange(-m_sp, m_sp + 1)
        wm = np.empty((nw, len(ms)))
        for im, m in enumerate(ms):
            ph = ph0 * np.exp(-1j * m * phi)[:, None]
            um = np.einsum("na,nia->ia", ph, Ut) / n_step
            wm[:, im] = np.real(np.einsum("ia,ia->a", um.conj(), um))
        i0 = int(np.where(ms == 0)[0][0])
        out["wspec"] = dict(ms=ms, wm=wm, w_sum=wm.sum(axis=1),
                            w0=wm[:, i0], w_main=wm.max(axis=1),
                            m_dom=ms[np.argmax(wm, axis=1)])
        out["micromotion"] = dict(Ut=Ut, t=tmid, T=T, dt=dt)
    return out


def quasienergies(model, drive, k, n_step=150, pmax=None):
    return floquet_modes(model, drive, k, n_step=n_step, pmax=pmax)["eps"]


def mean_energies(model, drive, k, n_step=150, pmax=None):
    return floquet_modes(model, drive, k, n_step=n_step, pmax=pmax)["Ebar"]
