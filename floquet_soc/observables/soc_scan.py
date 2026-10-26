# -*- coding: utf-8 -*-
"""Standard Floquet-SOC scan: alpha_R(gamma), Delta_spin(gamma), PST(gamma).

Definitions (kept explicit to avoid over-claiming):

* ``alpha_R``: effective valence-doublet slope along Gamma->K, fitted from
  period-averaged mean energies (same as ``floquet_points --rashba``).
* ``Delta_spin``: spin-related splitting of the valence pair and conduction
  pair at K/K', i.e. |Ebar_upper - Ebar_lower|, and its change w.r.t. static.
  It is *not* a microscopic SOC strength.
* ``PST``: collinearity matrix of unit spins over a small k region for the
  VBM band; largest eigenvalue (collinearity), principal axis, mean
  polarization, and axis rotation w.r.t. the static texture.

All three are computed for the same (gamma, polarization) grid and written to
one JSON/CSV plus figures, so the scan is directly comparable across materials.
"""
from __future__ import annotations

import csv
import json
import os

import numpy as np

from ..solvers.evolution import floquet_modes
from ..observables.spin_plot import plane_kpts
from ..observables.pst import pst_collinearity, axis_angle
from .floquet_points import (resolve_points, run_points, make_drive,
                             rashba_effective)
from .floquet_spin import floquet_spin_observables


def _edge_gap(E, lo, hi):
    """E[hi] - E[lo], NaN when the indices fall outside the model."""
    if not (0 <= lo < len(E) and 0 <= hi < len(E)):
        return float("nan")
    return float(E[hi] - E[lo])


def _static_edges(model, k, ne):
    E, U = np.linalg.eigh(model.Hk(k))
    return E, U


def _pst_region(model, spin_at, cfg, ne, hw, gamma, pol, n=7, kmax=0.35,
                n_step=120, m_sp=16):
    kpts = plane_kpts(cfg.get("kplane"), n, kmax)
    band = int(ne) - 1
    S = np.empty((len(kpts), 3))
    for i, k in enumerate(kpts):
        E, U = np.linalg.eigh(model.Hk(k))
        ops = np.asarray(spin_at(k))
        if gamma == 0:
            u = U[:, band]
            S[i] = np.real([u.conj() @ ops[c] @ u for c in range(3)])
        else:
            obs = floquet_spin_observables(
                model, make_drive(hw, gamma, pol, getattr(model, "a_ref", 1.0)),
                k, ops, n_step=n_step, m_sp=m_sp)
            j = int(np.argmax(np.abs(U[:, band].conj() @ obs["evec"])))
            S[i] = obs["S_bar"][j]
    pst = pst_collinearity(S)
    return dict(collinearity=float(pst["eigenvalues"][0]),
                axis=np.asarray(pst["pst_axis"]).tolist(),
                mean_polarization=float(pst["mean_polarization"]),
                n_states=int(pst["n_states"]),
                max_polarization=float(np.linalg.norm(S, axis=1).max()))


def run_scan(model, spin_at, cfg, hw, gammas, pols, points=None,
             n_step=120, m_sp=16, pst_n=7, pst_kmax=0.35):
    ne = cfg.get("num_electrons")
    if not ne:
        raise ValueError("soc_scan needs num_electrons")
    pts = resolve_points(cfg, points or "K,Kp,Gamma")
    # static K/K' splittings
    static = {}
    for lab, k in pts:
        E, U = _static_edges(model, k, ne)
        static[lab] = dict(
            k=np.asarray(k).tolist(),
            val=_edge_gap(E, int(ne) - 2, int(ne) - 1),
            cond=_edge_gap(E, int(ne), int(ne) + 1),
            E=E.tolist(), U=None)
    # Floquet targets at the same points
    fl = run_points(model, spin_at, cfg, pts, hw, gammas, pols,
                    n_step=n_step, m_sp=m_sp)
    # Delta_spin is defined on the gap-edge quadruple.  Small models (fewer
    # than 4 bands) may only provide part of it: report NaN there instead of
    # raising KeyError on a hard-coded name.
    def _split(t, hi, lo):
        if hi not in t or lo not in t:
            return float("nan")
        return abs(t[hi]["Ebar"] - t[lo]["Ebar"])

    rows_delta = []
    for lab in fl["points"]:
        for key, case in fl["points"][lab]["cases"].items():
            t = case["targets"]
            valF = _split(t, "VBM", "VBM-1")
            condF = _split(t, "CBM+1", "CBM")
            rows_delta.append(dict(
                point=lab, gamma=case["gamma"], pol=case["pol"],
                val_static=static[lab]["val"], val_F=valF,
                delta_val=valF - static[lab]["val"],
                cond_static=static[lab]["cond"], cond_F=condF,
                delta_cond=condF - static[lab]["cond"]))
    # alpha_R(gamma) from the Gamma->K effective doublet slope
    kp = [k for l, k in pts if l == "K"]
    direction = kp[0] if kp else np.array([1 / 3, 1 / 3, 0.0])
    alpha = {}
    for pol in pols:
        alpha[pol] = rashba_effective(
            model, spin_at, cfg, hw, gammas, direction_frac=direction,
            n_step=n_step, m_sp=m_sp, pol=pol)
    # PST baseline and driven
    base = _pst_region(model, spin_at, cfg, ne, hw, 0.0, "circ+",
                       n=pst_n, kmax=pst_kmax, n_step=n_step, m_sp=m_sp)
    pst_rows = [dict(region="VBM", gamma=0.0, pol="static",
                     collinearity=base["collinearity"],
                     axis=base["axis"],
                     mean_polarization=base["mean_polarization"],
                     n_states=base["n_states"],
                     axis_rotation_deg=0.0)]
    for gamma in gammas:
        if gamma == 0:
            continue
        for pol in pols:
            p = _pst_region(model, spin_at, cfg, ne, hw, gamma, pol,
                            n=pst_n, kmax=pst_kmax, n_step=n_step,
                            m_sp=m_sp)
            pst_rows.append(dict(
                region="VBM", gamma=float(gamma), pol=pol,
                collinearity=p["collinearity"], axis=p["axis"],
                mean_polarization=p["mean_polarization"],
                n_states=p["n_states"],
                axis_rotation_deg=float(axis_angle(
                    np.asarray(base["axis"]), np.asarray(p["axis"])))))
    return dict(hw=float(hw), gammas=[float(g) for g in gammas],
                pols=list(pols), points=[l for l, _ in pts],
                static=static, alpha=alpha,
                delta_spin=rows_delta, pst=pst_rows)


def save(res, outdir):
    os.makedirs(outdir, exist_ok=True)
    jp = os.path.join(outdir, "soc_scan.json")
    with open(jp, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2, ensure_ascii=False,
                  default=lambda o: o.tolist() if hasattr(o, "tolist")
                  else str(o))
    cp = os.path.join(outdir, "soc_scan.csv")
    with open(cp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["kind", "point/region", "gamma", "pol", "value1",
                    "value2", "value3"])
        for pol, rows in res["alpha"].items():
            for g, row in rows.items():
                w.writerow(["alpha_R", "Gamma-K", g, pol,
                            row["alpha_eff_ev_a"], row["r2"],
                            row["min_overlap"]])
        for r in res["delta_spin"]:
            w.writerow(["Delta_spin", r["point"], r["gamma"], r["pol"],
                        r["val_F"], r["delta_val"], r["cond_F"]])
        for r in res["pst"]:
            w.writerow(["PST", r["region"], r["gamma"], r["pol"],
                        r["collinearity"], r["axis_rotation_deg"],
                        r["mean_polarization"]])
    make_figures(res, outdir)
    return dict(json=jp, csv=cp,
                figure=os.path.join(outdir, "soc_scan.png"))


def make_figures(res, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .. import plot_style as ps
    ps.enforce_style()
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 2.9),
                             constrained_layout=True)
    # alpha_R
    ax = axes[0]
    for pol, rows in res["alpha"].items():
        g = [float(k) for k in rows]
        a = [rows[k]["alpha_eff_ev_a"] * 1000 for k in rows]
        ax.plot(g, a, "o-", label=pol)
    ax.set_xlabel("gamma")
    ax.set_ylabel("alpha_eff (meV A)")
    ax.legend(frameon=False, fontsize=7)
    # Delta_spin: K and K' valence pair, one line style per polarization
    ax = axes[1]
    styles = ["-", "--", ":", "-."]
    for ip, pol in enumerate(res["pols"]):
        for pt in ("K", "K′"):
            rows = [r for r in res["delta_spin"]
                    if r["point"] == pt and r["pol"] == pol]
            if not rows:
                continue
            rows = sorted(rows, key=lambda r: r["gamma"])
            g = [r["gamma"] for r in rows]
            dv = [r["delta_val"] * 1000 for r in rows]
            ax.plot(g, dv, marker="o", ls=styles[ip % len(styles)], lw=0.9,
                    ms=2.5, label=f"{pt} {pol}")
    ax.set_xlabel("gamma")
    ax.set_ylabel("delta Delta_spin (meV)")
    ax.legend(frameon=False, fontsize=5)
    # PST: include the static baseline at gamma=0
    ax = axes[2]
    base = next((r for r in res["pst"] if r["pol"] == "static"), None)
    for pol in res["pols"]:
        rows = [r for r in res["pst"] if r["pol"] == pol]
        if not rows:
            continue
        rows = sorted(rows, key=lambda r: r["gamma"])
        g = [r["gamma"] for r in rows]
        c = [r["collinearity"] for r in rows]
        if base is not None:
            g = [0.0] + g
            c = [base["collinearity"]] + c
        ax.plot(g, c, "o-", ms=2.5, label=pol)
    ax.set_xlabel("gamma")
    ax.set_ylabel("PST collinearity")
    if base is not None:
        vals = [base["collinearity"]] + [r["collinearity"] for r in res["pst"]
                                         if r["pol"] != "static"]
        lo, hi = min(vals), max(vals)
        pad = max(1e-4, 0.1 * (hi - lo))
        ax.set_ylim(lo - pad, hi + pad)
    ax.legend(frameon=False, fontsize=7)
    ps.save_figure(fig, os.path.join(outdir, "soc_scan.png"))
