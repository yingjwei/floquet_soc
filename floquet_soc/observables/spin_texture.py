# -*- coding: utf-8 -*-
"""Static / driven spin textures.

Conventions
-----------
* Spin operators are given in **Pauli units** (eigenvalue bound 1), i.e. the
  analytic models use ``S_i = sigma_i / 2`` for hbar units; when reading a
  ``.spn`` file the convention is the VASP one (also Pauli-like, bound 1).
  ``units`` is recorded in the returned dict and must be quoted with results.
* A *band-basis* spin texture is the diagonal of the spin matrix in the
  eigenbasis of H(k):  S_n(k) = <psi_nk| S |psi_nk>.
"""
from __future__ import annotations

import numpy as np


def band_eigenbasis(Hk: np.ndarray):
    ev, evec = np.linalg.eigh(np.asarray(Hk))
    return ev, evec


def spin_texture(Hk: np.ndarray, ops: np.ndarray) -> dict:
    """Spin texture of Hamiltonian eigenstates.

    Parameters
    ----------
    Hk : (nb,nb) complex Hermitian Hamiltonian at one k.
    ops : (3,nb,nb) spin operators in the same basis.

    Returns dict with energies (nb,), spins (nb,3), eigenvectors (nb,nb).
    """
    Hk = np.asarray(Hk)
    ops = np.asarray(ops)
    nb = Hk.shape[0]
    if ops.shape != (3, nb, nb):
        raise ValueError(f"spin operators shape {ops.shape} incompatible "
                         f"with H(k) {Hk.shape}")
    ev, evec = np.linalg.eigh(Hk)
    spins = np.empty((nb, 3), dtype=float)
    for c in range(3):
        # <psi_n| S_c |psi_n> = sum_ij evec[i,n]^* ops[c,i,j] evec[j,n]
        spins[:, c] = np.real(np.einsum("in,ij,jn->n", evec.conj(), ops[c],
                                        evec))
    return dict(energies=ev, spins=spins, evec=evec)


def spin_texture_path(model, ops, kpts) -> dict:
    """Spin texture along a list of k points.

    ``model`` must provide ``Hk(k)``; ``ops`` are the spin operators in the
    same basis (3,nb,nb).
    """
    kpts = np.atleast_2d(np.asarray(kpts, float))
    nk = len(kpts)
    out = None
    for ik, k in enumerate(kpts):
        r = spin_texture(model.Hk(k), ops)
        if out is None:
            nb = len(r["energies"])
            out = dict(energies=np.empty((nk, nb)),
                       spins=np.empty((nk, nb, 3)))
        out["energies"][ik] = r["energies"]
        out["spins"][ik] = r["spins"]
    out["kpts"] = kpts
    out["units"] = "Pauli (bound 1); divide by 2 for hbar units"
    return out


def spin_polarization_magnitude(spins: np.ndarray) -> np.ndarray:
    """|S_n(k)| for an array (nk,nb,3) or (nb,3)."""
    return np.linalg.norm(np.asarray(spins), axis=-1)


def spin_momentum_locking_angle(spins, kpts, band_indices=None):
    """Angle (rad) between spin and the in-plane momentum for each state.

    Only meaningful for 2D models; the angle is defined as the angle between
    the in-plane spin and k, so 90 degrees means S _|_ k (Rashba-like).
    """
    spins = np.asarray(spins)
    kpts = np.atleast_2d(np.asarray(kpts, float))
    if band_indices is None:
        band_indices = range(spins.shape[1])
    out = []
    for ik, k in enumerate(kpts):
        kn = np.array([k[0], k[1], 0.0])
        nk = np.linalg.norm(kn)
        if nk < 1e-12:
            continue
        for ib in band_indices:
            s = np.array([spins[ik, ib, 0], spins[ik, ib, 1], 0.0])
            sn = np.linalg.norm(s)
            if sn < 1e-12:
                continue
            out.append(float(np.arccos(np.clip(np.dot(s, kn) / (sn * nk),
                                              -1, 1))))
    return np.array(out)


def spin_texture_wannier_batch(model, spin_op, kpts, chunk=128):
    """Batched spin texture for a WannierSpinOperator (interpolated).

    S(k) = sum_R S(R) exp(i k.R) is built with a single BLAS matmul per chunk;
    H(k) is diagonalised in batches and the diagonal spin expectations
    <sigma_i>_b = <b|S_i(k)|b> are evaluated with matmuls.

    Returns dict(energies=(nk,nw), spins=(nk,nw,3), chunk=chunk).
    """
    kpts = np.atleast_2d(np.asarray(kpts, float))
    nk = len(kpts)
    nw = model.nw
    nR = len(spin_op.Rvecs)
    SRf = spin_op.S_R.reshape(nR, 3 * nw * nw)
    energies = np.empty((nk, nw))
    spins = np.empty((nk, nw, 3))
    for start in range(0, nk, chunk):
        stop = min(start + chunk, nk)
        Kb = kpts[start:stop]
        ph = np.exp(2j * np.pi * (Kb @ spin_op.Rvecs.T))       # (m,nR)
        ops = (ph @ SRf).reshape(len(Kb), 3, nw, nw)
        ops = 0.5 * (ops + ops.conj().transpose(0, 1, 3, 2))
        Hk = model.Hk_batch(Kb)
        w, V = np.linalg.eigh(Hk)
        energies[start:stop] = w
        for c in range(3):
            MV = np.matmul(ops[:, c], V)                       # (m,nw,nw)
            spins[start:stop, :, c] = np.real(
                np.einsum("kib,kib->kb", V.conj(), MV))
    return dict(energies=energies, spins=spins, chunk=chunk)
