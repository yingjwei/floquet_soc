# -*- coding: utf-8 -*-
"""Persistent spin texture (PST) metrics.

Definition
----------
For a set of target states (momentum region, band set) with unit spin
directions ``s_hat_i(k)``:

    M_ij = (1/N_k) sum_k s_hat_i(k) s_hat_j(k)

* perfect uniaxial collinearity: eigenvalues of M are (1, 0, 0); the
  eigenvector of the largest eigenvalue is the PST axis;
* isotropic (fully non-collinear): (1/3, 1/3, 1/3);
* in-plane Rashba: (1/2, 1/2, 0).

The metric measures **collinearity only**.  It is not a spin lifetime and
cannot be used to infer one without a scattering/dynamics model.
"""
from __future__ import annotations

import numpy as np


def pst_collinearity(spins, weights=None, min_polarization=1e-6):
    """Collinearity matrix M and its eigen-decomposition.

    Parameters
    ----------
    spins : (n,3) array of spin expectation values (any unit, e.g. Pauli).
    weights : optional (n,) weights (e.g. spectral weight).
    min_polarization : ignore states with |S| below this value.

    Returns dict with M, eigenvalues (descending), axes, collinearity
    (largest eigenvalue), pst_axis, n_states, mean_polarization.
    """
    S = np.atleast_2d(np.asarray(spins, float))
    norm = np.linalg.norm(S, axis=1)
    keep = norm > min_polarization
    if weights is None:
        w = np.ones(keep.sum())
    else:
        w = np.asarray(weights, float)[keep]
        w = w / w.sum()
    S = S[keep] / norm[keep, None]
    if len(S) < 3:
        raise ValueError("need at least 3 polarised states for PST metric")
    if weights is None:
        M = S.T @ S / len(S)
        mean_pol = float(np.mean(norm[keep]))
    else:
        M = (S * w[:, None]).T @ S
        mean_pol = float(np.sum(norm[keep] * w))
    ev, axes = np.linalg.eigh(M)                 # ascending
    order = np.argsort(ev)[::-1]
    ev = ev[order]
    axes = axes[:, order]
    return dict(M=M, eigenvalues=ev, axes=axes,
                collinearity=float(ev[0]), pst_axis=axes[:, 0],
                n_states=int(keep.sum()), mean_polarization=mean_pol)


def axis_angle(axis1, axis2):
    """Angle (rad) between two axes (sign-insensitive)."""
    a = np.asarray(axis1, float)
    b = np.asarray(axis2, float)
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    return float(np.arccos(np.clip(abs(np.dot(a, b)), -1.0, 1.0)))


def pst_change(pst_before, pst_after):
    """Compare two PST analyses: preservation / rotation / degradation."""
    rot = axis_angle(pst_before["pst_axis"], pst_after["pst_axis"])
    lam_b = pst_before["eigenvalues"]
    lam_a = pst_after["eigenvalues"]
    return dict(axis_rotation_deg=float(np.degrees(rot)),
                collinearity_before=float(lam_b[0]),
                collinearity_after=float(lam_a[0]),
                delta_collinearity=float(lam_a[0] - lam_b[0]),
                preserved=bool(rot < np.radians(5)
                                and abs(lam_a[0] - lam_b[0]) < 0.05))


def analyse_texture(energies, spins, bands=None, emin=None, emax=None,
                    min_polarization=0.1, min_splitting=0.0):
    """PST analysis on a selected (band, energy) region of a spin texture.

    Parameters
    ----------
    energies : (nk,nb) energies
    spins    : (nk,nb,3) spin expectation values (Pauli units)
    bands    : optional list of band indices; default all
    emin, emax : optional energy window
    min_polarization : ignore states below this |S| (default 0.1), because
        PST collinearity of weakly polarized states is meaningless.

    Returns dict of ``pst_collinearity`` plus selection bookkeeping:
        n_states, mean_polarization, bands, emin, emax, min_polarization
    """
    E = np.atleast_2d(np.asarray(energies, float))
    S = np.asarray(spins, float)
    nk, nb = E.shape
    if bands is None:
        bands = list(range(nb))
    bands = [b for b in bands if 0 <= b < nb]
    emin = -np.inf if emin is None else emin
    emax = np.inf if emax is None else emax
    vecs = []
    for ik in range(nk):
        for b in bands:
            if not (emin <= E[ik, b] <= emax):
                continue
            # exclude near-degenerate states: there the spin direction is
            # ill-defined (eigenvector mixing inside the degenerate subspace)
            others = np.delete(E[ik], b)
            split = float(np.min(np.abs(others - E[ik, b])))
            if split < min_splitting:
                continue
            vecs.append(S[ik, b])
    vecs = np.array(vecs)
    res = pst_collinearity(vecs, min_polarization=min_polarization)
    res.update(dict(bands=bands, emin=float(emin), emax=float(emax),
                    min_polarization=float(min_polarization),
                    min_splitting=float(min_splitting)))
    return res


def gauge_invariance_test(model, spin_at_k, kpts, bands, seeds=(1, 2, 3, 4, 5),
                          deg_tol=5e-3):
    """Gauge-invariance test for single-band spin / PST metrics.

    Within near-degenerate subspaces (|E_i-E_j| < deg_tol eV) the eigenvectors
    are arbitrary; random U(2) rotations must not change any physical
    observable.  Returns the spread of the PST collinearity / axis and the
    gauge-invariant subspace spin vectors  s_p = 0.5*Tr(P_p sigma).
    """
    import numpy as np
    kpts = np.atleast_2d(np.asarray(kpts, float))
    nb = int(model.nw)
    E = np.empty((len(kpts), nb))
    S0 = np.empty((len(kpts), nb, 3))
    V = np.empty((len(kpts), nb, nb), complex)
    OPS = np.empty((len(kpts), 3, nb, nb), complex)
    for ik, k in enumerate(kpts):
        ev, vec = np.linalg.eigh(model.Hk(k))
        ops_k = spin_at_k(k)
        E[ik] = ev
        V[ik] = vec
        OPS[ik] = ops_k
        for c in range(3):
            M = np.conj(vec.T) @ ops_k[c] @ vec
            S0[ik, :, c] = np.real(np.diag(M))

    def metric(S):
        r = analyse_texture(E, S, bands=bands)
        return dict(M=r["eigenvalues"].tolist(),
                    coll=float(r["collinearity"]),
                    axis=np.asarray(r["pst_axis"], float).tolist())

    base = metric(S0)
    # near-degenerate pairs (sorted eigenvector columns (0,1),(2,3),...)
    pair_list = []
    for ik in range(len(kpts)):
        ev = E[ik]
        for a in range(0, nb - 1, 2):
            if abs(ev[a] - ev[a + 1]) < deg_tol:
                pair_list.append((ik, a, a + 1))
    # gauge-invariant subspace spin
    s_inv = []
    for (ik, a, b) in pair_list:
        vec = V[ik][:, [a, b]]
        P = vec @ np.conj(vec.T)
        s_inv.append([0.5 * np.real(np.trace(P @ OPS[ik, c]))
                      for c in range(3)])
    s_inv = np.asarray(s_inv, float) if s_inv else np.zeros((0, 3))
    seeds_out = []
    rng = np.random.default_rng
    for seed in seeds:
        rr = rng(seed)
        S = S0.copy()
        for ik in range(len(kpts)):
            vec = V[ik].copy()
            for (ik2, a, b) in pair_list:
                if ik2 != ik:
                    continue
                W = rr.normal(size=(2, 2)) + 1j * rr.normal(size=(2, 2))
                W, _ = np.linalg.qr(W)
                vec[:, [a, b]] = V[ik][:, [a, b]] @ W
            for c in range(3):
                M = np.conj(vec.T) @ OPS[ik, c] @ vec
                S[ik, :, c] = np.real(np.diag(M))
        m = metric(S)
        m["seed"] = int(seed)
        seeds_out.append(m)
    coll = [m["coll"] for m in seeds_out] or [base["coll"]]
    ax0 = np.asarray(base["axis"], float)
    ax0 = ax0 / max(np.linalg.norm(ax0), 1e-12)
    dots = []
    for m in seeds_out:
        ax = np.asarray(m["axis"], float)
        dots.append(abs(float(np.dot(ax / max(np.linalg.norm(ax), 1e-12),
                                     ax0))))
    return dict(base=base, seeds=seeds_out,
                n_pairs=len(pair_list),
                coll_spread=float(max(coll) - min(coll)),
                coll_range=[float(min(coll)), float(max(coll))],
                axis_dot_min=(float(min(dots)) if dots else 1.0),
                invariant_spin=s_inv,
                invariant_mean_norm=float(np.linalg.norm(s_inv, axis=1).mean())
                if len(s_inv) else 0.0,
                invariant_max_norm=float(np.linalg.norm(s_inv, axis=1).max())
                if len(s_inv) else 0.0)
