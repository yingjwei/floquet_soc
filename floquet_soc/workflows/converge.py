# -*- coding: utf-8 -*-
"""Stage D: convergence study for the Floquet solvers.

For one model / drive / k point, vary the time-step number ``n_step`` and the
photon-window ``m_sp`` and compare against the most accurate run.  Convergence
is judged on quantities that are *not* affected by quasienergy folding:
mean energies ``Ebar``, periodic-average spin ``Sbar`` and photon weight sum.
A PASS requires monotone-ish decrease and final deviations below tolerance.
"""
from __future__ import annotations

import numpy as np

from ..observables.floquet_spin import floquet_spin_observables


def _invariants(obs):
    """Basis-mixing invariant sums (trace of A/T and of the spin operators)."""
    return dict(
        sum_ebar=float(np.sum(np.asarray(obs["Ebar"], float))),
        sum_sbar=np.asarray(np.sum(obs["S_bar"], axis=0), float),
        sum_w=float(np.sum(np.asarray(obs["w_sum"], float))),
    )


def _matched_differences(ref, obs, window, min_purity=0.999,
                         gap_tol=1e-3):
    """Match reference states to ``obs`` by overlap, excluding near-degenerate
    states whose eigenvectors are mixed by the Schur decomposition.

    Returns (dev_ebar, dev_spin, n_matched, min_purity_used).
    """
    E_ref = np.asarray(ref["Ebar"], float)
    S_ref = np.asarray(ref["S_bar"], float)
    ev_ref = np.asarray(ref["evec"], complex)
    eps_ref = np.sort(np.asarray(ref["eps"], float))
    E = np.asarray(obs["Ebar"], float)
    S = np.asarray(obs["S_bar"], float)
    ev = np.asarray(obs["evec"], complex)
    ov = np.abs(ev_ref.conj().T @ ev) ** 2          # (ref, run)
    best = np.argmax(ov, axis=1)
    pur = ov[np.arange(len(best)), best]
    # nearest quasienergy neighbour gap of the reference state
    gaps = np.full(len(eps_ref), np.inf)
    gaps[1:] = np.minimum(gaps[1:], np.diff(eps_ref))
    gaps[:-1] = np.minimum(gaps[:-1], np.diff(eps_ref))
    sel = ((E_ref > window[0]) & (E_ref < window[1])
           & (pur > min_purity) & (gaps > gap_tol))
    if not np.any(sel):
        return np.nan, np.nan, np.nan, 0, 0.0
    dev_e = float(np.max(np.abs(E[best[sel]] - E_ref[sel])))
    dev_s = float(np.max(np.abs(S[best[sel]] - S_ref[sel])))
    # photon-window-sensitive observable: sum of photon weights per state
    w_ref = np.asarray(ref["w_sum"], float)
    w_obs = np.asarray(obs["w_sum"], float)
    dev_w = float(np.max(np.abs(w_obs[best[sel]] - w_ref[sel])))
    return dev_e, dev_s, dev_w, int(sel.sum()), float(np.min(pur[sel]))


def run_convergence(model, drive, k, spin_ops, n_step_values, m_sp_values,
                    tol_ebar=1e-5, tol_spin=1e-4, window=(-2.0, 2.0),
                    min_purity=0.5, tol_wsum=1e-6):
    """Run n_step and m_sp convergence checks at one k point.

    Returns dict with:
        n_step_rows, m_sp_rows (lists of dict), reference (n_step, m_sp),
        pass (bool), max_dev_ebar, max_dev_spin, wsum_err
    """
    n_step_values = sorted(int(x) for x in n_step_values)
    m_sp_values = sorted(int(x) for x in m_sp_values)
    ns_ref, ms_ref = n_step_values[-1], m_sp_values[-1]
    ref = floquet_spin_observables(model, drive, k, spin_ops,
                                   n_step=ns_ref, m_sp=ms_ref)

    def _wsum(obs):
        E = np.asarray(obs["Ebar"], float)
        w = np.asarray(obs["w_sum"], float)
        sel = (E > window[0]) & (E < window[1])
        return float(np.max(np.abs(w[sel] - 1.0))) if np.any(sel) else np.nan

    inv_ref = _invariants(ref)

    n_rows = []
    for ns in n_step_values:
        obs = floquet_spin_observables(model, drive, k, spin_ops,
                                       n_step=ns, m_sp=ms_ref)
        de, ds, dw, nm, pur = _matched_differences(ref, obs, window,
                                                    min_purity)
        n_rows.append(dict(n_step=ns, dev_ebar=de, dev_spin=ds, dev_wsum=dw,
                           n_matched=nm, min_purity=pur,
                           wsum_err=_wsum(obs), _invariants=_invariants(obs),
                           unitarity_error=float(obs["unitarity_error"])))
    m_rows = []
    for ms in m_sp_values:
        obs = floquet_spin_observables(model, drive, k, spin_ops,
                                       n_step=ns_ref, m_sp=ms)
        de, ds, dw, nm, pur = _matched_differences(ref, obs, window,
                                                    min_purity)
        # per-state photon-weight deficit in the target energy window
        Ebar = np.asarray(obs["Ebar"], float)
        wsum = np.asarray(obs["w_sum"], float)
        sel = (Ebar > window[0]) & (Ebar < window[1])
        deficit = 1.0 - wsum[sel]
        order = np.argsort(-np.abs(deficit))
        idx = np.where(sel)[0][order]
        worst = [dict(index=int(i), Ebar=float(Ebar[i]),
                      deficit=float(1.0 - wsum[i]))
                 for i in idx[:5]]
        m_rows.append(dict(m_sp=ms, dev_ebar=de, dev_spin=ds, dev_wsum=dw,
                           n_matched=nm, min_purity=pur,
                           wsum_err=_wsum(obs),
                           w_deficit_max=float(np.max(np.abs(deficit)))
                           if deficit.size else 0.0,
                           w_deficit_median=float(np.median(np.abs(deficit)))
                           if deficit.size else 0.0,
                           n_deficit_gt_1e3=int(np.sum(np.abs(deficit) > 1e-3)),
                           n_states_window=int(deficit.size),
                           worst_states=worst,
                           _invariants=_invariants(obs)))
    for row in n_rows + m_rows:
        inv = row.pop("_invariants")
        row["dev_sum_ebar"] = abs(inv["sum_ebar"] - inv_ref["sum_ebar"])
        row["dev_sum_sbar"] = float(np.max(np.abs(inv["sum_sbar"]
                                                  - inv_ref["sum_sbar"])))
        row["dev_sum_w"] = abs(inv["sum_w"] - inv_ref["sum_w"])
    used = [r for r in n_rows + m_rows if r["n_matched"] > 0]
    dev_ebar = max([r["dev_ebar"] for r in used]) if used else np.nan
    dev_spin = max([r["dev_spin"] for r in used]) if used else np.nan
    dev_inv = max([max(r["dev_sum_ebar"], r["dev_sum_sbar"])
                   for r in n_rows + m_rows])
    dev_wsum = max([r["dev_wsum"] for r in n_rows + m_rows
                    if np.isfinite(r["dev_wsum"])] or [np.nan])
    # "production" setting = the largest setting below the reference run
    prod = []
    if len(n_rows) > 1:
        prod.append(n_rows[-2])
    if len(m_rows) > 1:
        prod.append(m_rows[-2])
    prod_ok = all((r["n_matched"] > 0 and r["dev_ebar"] < tol_ebar
                   and r["dev_spin"] < tol_spin
                   and r["dev_wsum"] < tol_wsum) for r in prod)
    ok_m = [r["m_sp"] for r in m_rows
            if r.get("w_deficit_max", 0.0) < 1e-3]
    accepted_m_sp = (min(ok_m) if ok_m else None)
    passed = bool(prod_ok and (dev_inv < 1e-6) and dev_wsum < tol_wsum)
    return dict(n_step_rows=n_rows, m_sp_rows=m_rows,
                reference=dict(n_step=ns_ref, m_sp=ms_ref),
                invariants=inv_ref, pass_=passed, max_dev_ebar=dev_ebar,
                max_dev_spin=dev_spin, max_dev_invariant=dev_inv,
                max_dev_wsum=dev_wsum,
                accepted_m_sp=accepted_m_sp,
                production=[r.get("n_step", r.get("m_sp")) for r in prod],
                wsum_err=max(r["wsum_err"] for r in n_rows + m_rows))
