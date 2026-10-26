# -*- coding: utf-8 -*-
"""SOC-related band splitting and effective-SOC fitting (stage C).

Definitions
-----------
Static splitting:
    Delta^0_{n}(k) = E_{n+1}(k) - E_n(k)          (adjacent bands)
or, for a user-selected spin-split pair, the splitting of the two bands that
are closest in energy *and* have opposite spin expectation (pairing rule is
explicit and reported).

Driven case:
    Delta^F_{n}(k) = quasienergy splitting of the same pair
                     (after overlap matching to static bands)
    delta Delta(k) = Delta^F - Delta^0

Effective SOC fitting
---------------------
A *generic* linear model is fitted in a user-chosen low-energy window:

    H(k) ~= sum_j c_j f_j(k) M_j

with basis functions f_j(k) and operator matrices M_j supplied by the user
(e.g. {1, kx, ky, kx*ky, ...} x {I, sigma_x, sigma_y, sigma_z} for a two-band
model).  The fit is a least-squares problem over the selected k region; the
residual is reported so the model choice can be judged.  No universal
Rashba/Dresselhaus formula is imposed on arbitrary materials.

SOC on/off
----------
Two Hamiltonians may be subtracted only if they share the Wannier basis and
gauge.  ``soc_off_difference`` refuses direct subtraction otherwise and
returns NEEDS_INPUT with the missing files listed.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares


# ------------------------------------------------------------------ splitting
def adjacent_splitting(energies):
    """Delta_n(k) = E_{n+1}(k) - E_n(k); energies (nk, nb)."""
    E = np.atleast_2d(np.asarray(energies, float))
    return np.diff(E, axis=1)


def spin_paired_splitting(energies, spins, emin=None, emax=None,
                          max_energy_gap=0.5):
    """Pair neighbouring bands by (i) small energy gap and (ii) opposite spin.

    Parameters
    ----------
    energies : (nk,nb) band energies
    spins    : (nk,nb,3) spin expectation (Pauli units) of the same bands

    Returns dict with arrays (nk, npair):
        i, j, gap, spin_antiparallel, energy_mid
    Only pairs with both energies inside [emin, emax] and gap <= max_energy_gap
    are kept; the pairing is done per k (greedy nearest-energy neighbour with
    anti-parallel spin preference), and the rule is recorded.
    """
    E = np.atleast_2d(np.asarray(energies, float))
    S = np.asarray(spins, float)
    nk, nb = E.shape
    if S.shape[:2] != (nk, nb):
        raise ValueError("energies/spins shape mismatch")
    emin = -np.inf if emin is None else emin
    emax = np.inf if emax is None else emax
    pairs = []
    for ik in range(nk):
        order = np.argsort(E[ik])
        used = np.zeros(nb, bool)
        for a in range(nb - 1):
            n = order[a]
            if used[n]:
                continue
            if not (emin <= E[ik, n] <= emax):
                continue
            # candidates: nearest energies above, prefer antiparallel spin
            cand = [m for m in order[a + 1:] if not used[m]
                    and E[ik, m] - E[ik, n] <= max_energy_gap]
            if not cand:
                continue
            sn = S[ik, n]
            scores = []
            for m in cand:
                sm = S[ik, m]
                denom = max(np.linalg.norm(sn) * np.linalg.norm(sm), 1e-12)
                cos = float(np.dot(sn, sm) / denom)
                gap = float(E[ik, m] - E[ik, n])
                scores.append((gap / max(max_energy_gap, 1e-12)
                               + 0.5 * (cos + 1.0), m, gap, cos))
            score, m, gap, cos = min(scores)
            used[n] = used[m] = True
            pairs.append(dict(k=ik, i=int(n), j=int(m), gap=gap,
                              spin_antiparallel=(-cos),
                              energy_mid=float(0.5 * (E[ik, n] + E[ik, m]))))
    return pairs


# ------------------------------------------------------------------ fitting
def fit_linear_model(Hk_func, kpts, basis_funcs, operators,
                     weights=None, k_weight_func=None):
    """Least-squares fit  H(k) ~= sum_j c_j f_j(k) M_j.

    Parameters
    ----------
    Hk_func : callable k -> (n,n) matrix (Hermitian)
    kpts : (nk,3) points used for the fit
    basis_funcs : list of callables f_j(k) -> float
    operators : list of (n,n) matrices M_j (same length)
    weights : optional array (nk,) of weights
    k_weight_func : optional callable k -> float extra weight

    Returns dict: coefficients c, residual (Frobenius), rms, condition number,
    n_terms, kpts.
    """
    kpts = np.atleast_2d(np.asarray(kpts, float))
    Hk = [np.asarray(Hk_func(k), complex) for k in kpts]
    nb = Hk[0].shape[0]
    nk = len(kpts)
    nterm = len(basis_funcs)
    if len(operators) != nterm:
        raise ValueError("basis_funcs and operators must have the same length")
    for op in operators:
        if np.asarray(op).shape != (nb, nb):
            raise ValueError("operator dimension mismatch")
    M = np.stack([np.asarray(op, complex) for op in operators])
    H = np.stack(Hk)                                   # (nk,nb,nb)
    w = np.ones(nk)
    if weights is not None:
        w *= np.asarray(weights, float)
    if k_weight_func is not None:
        w *= np.array([float(k_weight_func(k)) for k in kpts])
    # element-wise least squares:  H(k)[a,b] ~= sum_j c_j f_j(k) M_j[a,b]
    X = np.zeros((nk * nb * nb, nterm), complex)
    for j in range(nterm):
        f = np.array([float(basis_funcs[j](k)) for k in kpts])
        X[:, j] = (w[:, None, None] * f[:, None, None]
                   * M[j][None, :, :]).reshape(-1)
    Y = (w[:, None, None] * H).reshape(-1)
    c, *_ = np.linalg.lstsq(X, Y, rcond=None)
    c = np.real_if_close(c, tol=1e4).real
    model = np.einsum("j,jab->ab", c, M)
    resid = np.array([np.linalg.norm(w[ik] * (H[ik] - sum(
        c[j] * float(basis_funcs[j](kpts[ik])) * M[j]
        for j in range(nterm)))) for ik in range(nk)])
    rms = float(np.sqrt(np.mean(resid ** 2)))
    return dict(coefficients=c, residual=resid, rms=rms,
                condition=float(np.linalg.cond(X)), n_terms=nterm,
                kpts=kpts)


def soc_off_difference(scalar_config, soc_config):
    """Refuse to subtract Hamiltonians from different basis/gauge runs.

    Returns a dict with status 'NEEDS_INPUT' unless both configurations point
    to the *same* Wannier basis (same win/hr/centres hashes) or a matched
    linear interpolation H(lambda)=H_scalar+lambda*H_SOC is provided.
    """
    import hashlib
    import os

    def digest(cfg):
        files = cfg.get("files", {})
        h = hashlib.sha1()
        for key in ("win", "hr", "centres"):
            p = files.get(key)
            if not (p and os.path.exists(p)):
                return None
            f = open(p, "rb")
            h.update(f.read(1 << 16))
            f.close()
        return h.hexdigest()

    if scalar_config is None or soc_config is None:
        return dict(status="NEEDS_INPUT",
                    reason="provide SOC-on and SOC-off Wannier runs "
                           "(and confirm they share the same Wannier basis)")
    d1, d2 = digest(scalar_config), digest(soc_config)
    if d1 is None or d2 is None:
        return dict(status="NEEDS_INPUT",
                    reason="Wannier files (win/hr/centres) are missing; cannot "
                           "verify a common basis")
    if d1 != d2:
        return dict(status="NEEDS_INPUT",
                    reason="SOC-on and SOC-off Wannier bases differ "
                           "(win/hr/centres hashes do not match); direct "
                           "subtraction would be meaningless")
    return dict(status="OK",
                reason="same Wannier basis: H(lambda)=H_scalar+lambda*H_SOC "
                       "may be constructed")


# ------------------------------------------------------- gap-edge pairs
def gap_edge_pairs(num_electrons, num_wann):
    """0-based band pairs that straddle the band gap.

    For ``num_electrons`` filled spinor bands the physically relevant SOC
    splittings near the gap are only:

    * valence pair  ``(VBM-1, VBM)``   = ``(ne-2, ne-1)``;
    * conduction pair ``(CBM, CBM+1)`` = ``(ne, ne+1)``.

    Returns the list of pairs that exist in the model (0 <= index < num_wann).
    """
    ne = int(num_electrons)
    nw = int(num_wann)
    pairs = []
    if 0 <= ne - 2 and ne - 1 < nw:
        pairs.append((ne - 2, ne - 1))
    if 0 <= ne and ne + 1 < nw:
        pairs.append((ne, ne + 1))
    return pairs


def edge_pair_splitting(energies, pairs, spins=None):
    """Per-k spectral splitting of explicitly chosen band pairs.

    Unlike :func:`spin_paired_splitting` (which scans/greedily pairs adjacent
    bands), this evaluates exactly the requested ``(i, j)`` pairs, e.g. the two
    gap-edge pairs ``(VBM-1, VBM)`` and ``(CBM, CBM+1)``.

    Parameters
    ----------
    energies : (nk, nb) band energies
    pairs : iterable of (i, j) 0-based band indices
    spins : optional (nk, nb, 3) spin expectation values (Pauli units)

    Returns
    -------
    list of dict, one per pair:
        i, j, gap (nk,), max, min, median, k_argmax;
        if spins is given: spin_cos (nk,), spin_antiparallel (nk,),
        spin_antiparallel_median.
    """
    E = np.atleast_2d(np.asarray(energies, float))
    S = None if spins is None else np.asarray(spins, float)
    out = []
    for (i, j) in pairs:
        i, j = int(i), int(j)
        if not (0 <= i < E.shape[1] and 0 <= j < E.shape[1]):
            continue
        gap = E[:, j] - E[:, i]
        rec = dict(i=i, j=j, gap=gap, max=float(gap.max()),
                   min=float(gap.min()), median=float(np.median(gap)),
                   k_argmax=int(np.argmax(gap)))
        if S is not None:
            na = np.linalg.norm(S[:, i], axis=1)
            nb = np.linalg.norm(S[:, j], axis=1)
            cos = np.einsum("ia,ia->i", S[:, i], S[:, j]) / np.maximum(
                na * nb, 1e-12)
            rec["spin_cos"] = cos
            rec["spin_antiparallel"] = -cos
            rec["spin_antiparallel_median"] = float(np.median(-cos))
        out.append(rec)
    return out
