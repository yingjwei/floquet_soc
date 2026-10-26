# -*- coding: utf-8 -*-
"""Stage D: generic parameter scans over drive parameters.

The scan runs a *cheap, well-defined* observable at one or a few k points:
sorted quasienergies (spread), periodic-average spin polarization |S_bar|, and
the photon-weight normalisation error. It is deliberately not a "bands vs hw"
plot: for each (hw, gamma) the Floquet spectrum at fixed k is compared with
itself and with the static spin texture, so no band tracking is involved.
"""
from __future__ import annotations

import numpy as np

from ..drives.peierls import Drive
from ..observables.floquet_spin import floquet_spin_observables


def run_scan(model, drive_cfg, hw_values, gamma_values, kpts, spin_ops,
             n_step=150, m_sp=16, window=(-2.0, 2.0)):
    """Scan drive frequency and amplitude at fixed k points.

    Parameters
    ----------
    drive_cfg : dict passed to ``Drive.from_config`` minus hw/gamma
    hw_values, gamma_values : lists of floats (eV, A0*a_ref)
    kpts : (nk,3) fractional k points
    spin_ops : (3,nw,nw) spin operators (Pauli units) or a callable k->ops
    """
    kpts = np.atleast_2d(np.asarray(kpts, float))
    rows = []
    for hw in hw_values:
        for gamma in gamma_values:
            cfg = dict(drive_cfg)
            cfg.update(dict(hw=float(hw), gamma=float(gamma)))
            drive = Drive.from_config(cfg, a_ref=getattr(model, "a_ref", 1.0))
            for ik, k in enumerate(kpts):
                ops_k = spin_ops(k) if callable(spin_ops) else spin_ops
                obs = floquet_spin_observables(
                    model, drive, k, ops_k, n_step=n_step, m_sp=m_sp)
                eps = np.sort(np.asarray(obs["eps"], float))
                Ebar = np.asarray(obs["Ebar"], float)
                S_bar = np.asarray(obs["S_bar"], float)
                wsum = np.asarray(obs["w_sum"], float)
                sel = (Ebar > window[0]) & (Ebar < window[1])
                if np.any(sel):
                    pol = np.linalg.norm(S_bar[sel], axis=1)
                    wsum_err = float(np.max(np.abs(wsum[sel] - 1.0)))
                    mean_sz = float(np.mean(S_bar[sel, 2]))
                    n_win = int(sel.sum())
                else:
                    pol = np.array([np.nan])
                    wsum_err = np.nan
                    mean_sz = np.nan
                    n_win = 0
                rows.append(dict(
                    hw=float(hw), gamma=float(gamma), ik=int(ik),
                    kpt=[float(x) for x in k],
                    eps_spread=float(eps.max() - eps.min()),
                    eps_min=float(eps.min()), eps_max=float(eps.max()),
                    mean_pol=float(np.mean(pol)),
                    max_pol=float(np.max(pol)),
                    mean_sz=float(mean_sz),
                    wsum_err=float(wsum_err),
                    n_states_in_window=n_win,
                    unitarity_error=float(obs["unitarity_error"]),
                ))
    return rows


def rows_to_arrays(rows):
    keys = ("hw", "gamma", "eps_spread", "mean_pol", "max_pol", "mean_sz",
            "wsum_err", "unitarity_error")
    return {k: np.array([r[k] for r in rows], float) for k in keys}


def scan_grid(rows):
    hws = np.unique([r["hw"] for r in rows])
    gammas = np.unique([r["gamma"] for r in rows])
    return hws, gammas


def observable_map(rows, quantity="mean_pol", ik=0):
    """(ngamma, nhw) map of one observable for plotting (gamma rows, hw cols)."""
    hws, gammas = scan_grid(rows)
    Z = np.full((len(gammas), len(hws)), np.nan)
    for r in rows:
        if r["ik"] != ik:
            continue
        i = int(np.argmin(np.abs(gammas - r["gamma"])))
        j = int(np.argmin(np.abs(hws - r["hw"])))
        Z[i, j] = r[quantity]
    return hws, gammas, Z
