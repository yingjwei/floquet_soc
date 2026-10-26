# -*- coding: utf-8 -*-
"""Linear Rashba fitting for gap-edge band pairs.

Ported from the user's standalone script
``C:\\Users\\15756\\Desktop\\_桌面归档\\Python脚本\\rashba_refined_fitting.py``
(originally written for VASP ``BAND.dat`` + ``KLABELS``), adapted to the
Wannier model used by this framework.

Method
------
For a Kramers pair that is degenerate at a high-symmetry anchor (e.g. Γ) the
two branches disperse as

    E_a(k) = a·k² + E0 - α·|k|,      E_b(k) = a·k² + E0 + α·|k|,

so a parabolic background plus a **linear** term in |k| is fitted to each
branch over |k| <= kmax (k measured from the anchor, units 1/Å).  The fit
window is grown from kmin while both branches keep R² >= r2_thr AND the two
fitted |α| agree within alpha_ratio_thr; the largest accepted window is used
and α_R = (|α_a| + |α_b|)/2 in eV·Å.

If no window satisfies the criteria the result is reported as NOT_LINEAR
(the caller decides whether to emit WARN / NOT_APPLICABLE).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import curve_fit


def rashba_parabolic_model(k, a, e0, alpha):
    """E(k) = a·k² + E0 + α·|k| (k measured from the Rashba-degenerate point)."""
    k = np.asarray(k, float)
    return a * k * k + e0 + alpha * np.abs(k)


def fit_branch(k, e, kmax):
    """Least-squares fit of one branch over |k| <= kmax; returns (popt, r2)."""
    k = np.asarray(k, float)
    e = np.asarray(e, float)
    mask = np.abs(k) <= kmax
    kk, ee = k[mask], e[mask]
    if kk.size < 4:
        raise ValueError(f"only {kk.size} points inside |k| <= {kmax}")
    p0 = (0.5, float(ee.mean()), 0.0)
    popt, _ = curve_fit(rashba_parabolic_model, kk, ee, p0=p0, maxfev=20000)
    resid = ee - rashba_parabolic_model(kk, *popt)
    ss = float(np.sum(resid ** 2))
    tot = float(np.sum((ee - ee.mean()) ** 2))
    r2 = 1.0 - ss / tot if tot > 0 else 0.0
    return popt, r2


def splitting_linearity(k, e_lower, e_upper, kmax):
    """Fit ΔE = 2·α_R·|k| (through the origin) and return (alpha_R, R²).

    This is the direct "is the pair splitting linear?" test; it is used as a
    hard validity criterion in addition to the branch fits, because two nearly
    degenerate branches can be fitted with arbitrary same-sign α at high R².
    """
    k = np.asarray(k, float)
    delta = np.asarray(e_upper, float) - np.asarray(e_lower, float)
    mask = np.abs(k) <= kmax
    kk, dd = np.abs(k[mask]), delta[mask]
    if kk.size < 4 or np.sum(kk * kk) <= 0:
        return None, -np.inf
    s = float(np.sum(kk * dd) / np.sum(kk * kk))     # ΔE ≈ s·|k|, α_R = s/2
    resid = dd - s * kk
    ss = float(np.sum(resid ** 2))
    tot = float(np.sum((dd - dd.mean()) ** 2))
    r2 = 1.0 - ss / tot if tot > 0 else 0.0
    return 0.5 * s, r2


def auto_fit_window(k, e_lower, e_upper, kmin=0.03, kmax_max=0.15,
                    step=0.005, r2_thr=0.999, alpha_ratio_thr=0.15,
                    split_r2_thr=0.99):
    """Largest |k| window where the pair has a genuine linear Rashba splitting.

    Requirements for a window to be accepted:
      1. both branch fits have R² >= r2_thr;
      2. the two fitted α have OPPOSITE signs (E_± = a k² + E0 ± α|k|);
      3. the two |α| agree within alpha_ratio_thr;
      4. the direct splitting ΔE = 2·α_R·|k| fit has R² >= split_r2_thr.

    Returns a dict on success, or None if even the first window fails.
    """
    k = np.asarray(k, float)
    lo = np.asarray(e_lower, float)
    up = np.asarray(e_upper, float)
    best = None
    for kmax in np.arange(kmin, kmax_max + 1e-12, step):
        if np.sum(np.abs(k) <= kmax) < 4:
            continue
        try:
            p_lo, r2_lo = fit_branch(k, lo, kmax)
            p_up, r2_up = fit_branch(k, up, kmax)
        except (RuntimeError, ValueError):
            break
        if r2_lo < r2_thr or r2_up < r2_thr:
            break
        a_lo, a_up = abs(p_lo[2]), abs(p_up[2])
        if max(a_lo, a_up) <= 0.0:
            break
        if p_lo[2] * p_up[2] >= 0.0:        # no Rashba-like ±α splitting
            break
        ratio = abs(a_lo - a_up) / max(a_lo, a_up)
        if ratio > alpha_ratio_thr:
            break
        alpha_split, r2_split = splitting_linearity(k, lo, up, kmax)
        if alpha_split is None or r2_split < split_r2_thr:
            break
        alpha_r = 0.5 * (a_lo + a_up)
        if alpha_split > 0 and abs(alpha_split - alpha_r) > \
                alpha_ratio_thr * max(alpha_split, 1e-12):
            break
        best = dict(kmax=float(kmax),
                    alpha_r=float(alpha_r),
                    alpha_from_split=float(alpha_split),
                    r2_split=float(r2_split),
                    alpha_lower=float(p_lo[2]),
                    alpha_upper=float(p_up[2]),
                    r2_lower=float(r2_lo),
                    r2_upper=float(r2_up),
                    popt_lower=p_lo,
                    popt_upper=p_up)
    return best


def fit_gap_edge_pair(k, e_pair, kmin=0.03, kmax_max=0.15, step=0.005,
                      r2_thr=0.999, alpha_ratio_thr=0.15,
                      split_r2_thr=0.99, fallback_kmax=0.08):
    """Fit both branches of one gap-edge pair; returns a result dict.

    ``e_pair`` is (2, nk): the two branch energies on the same k grid.
    The lower/upper branch is ordered by energy at the smallest |k| point.

    Status:
      * ``LINEAR``   - a window satisfying all linearity criteria was found;
      * ``FALLBACK`` - no valid window; the original script would still fit at
        ``fallback_kmax`` (with a printed warning), so we do the same and keep
        the fit curves/alpha_R but expose the poor R² and status;
      * ``NOT_LINEAR`` - fallback fit itself failed (e.g. too few points).
    """
    k = np.asarray(k, float)
    e = np.atleast_2d(np.asarray(e_pair, float))
    order = np.argsort(np.abs(k))
    small = e[:, order[:2]].mean(axis=1)
    idx = np.argsort(small)
    lower, upper = e[idx[0]], e[idx[1]]
    res = auto_fit_window(k, lower, upper, kmin=kmin, kmax_max=kmax_max,
                          step=step, r2_thr=r2_thr,
                          alpha_ratio_thr=alpha_ratio_thr,
                          split_r2_thr=split_r2_thr)
    if res is not None:
        out = dict(status="LINEAR", k=k, lower=lower, upper=upper)
        out.update(res)
        return out
    # --- fallback identical to the original script's behaviour -------------
    try:
        p_lo, r2_lo = fit_branch(k, lower, fallback_kmax)
        p_up, r2_up = fit_branch(k, upper, fallback_kmax)
        alpha_split, r2_split = splitting_linearity(k, lower, upper,
                                                   fallback_kmax)
        a_lo, a_up = abs(p_lo[2]), abs(p_up[2])
        return dict(status="FALLBACK", k=k, lower=lower, upper=upper,
                    kmax=float(fallback_kmax),
                    alpha_r=float(0.5 * (a_lo + a_up)),
                    alpha_from_split=float(alpha_split)
                    if alpha_split is not None else float("nan"),
                    r2_split=float(r2_split),
                    alpha_lower=float(p_lo[2]),
                    alpha_upper=float(p_up[2]),
                    r2_lower=float(r2_lo),
                    r2_upper=float(r2_up),
                    popt_lower=p_lo, popt_upper=p_up)
    except (RuntimeError, ValueError):
        return dict(status="NOT_LINEAR", k=k, lower=lower, upper=upper,
                    alpha_r=None, kmax=None)
