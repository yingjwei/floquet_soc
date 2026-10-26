# -*- coding: utf-8 -*-
"""Band-resolved Floquet spin scan at selected k points.

Generic program behind ``floquet.py floquet_points``.  For every
(point, gamma, polarization) it computes the driven Floquet modes and matches
the four gap-edge bands (VBM-1, VBM, CBM, CBM+1) by overlap with the undriven
eigenvectors.  Reported per band:

    Ebar (period-averaged energy), S_bar, dS_bar vs static,
    dominant photon weight w_main, w_sum, micromotion, overlap.

Also aggregates valley asymmetry between K and K' when both are present.
Nothing here assumes a material; only cfg['num_electrons'], cfg['kpath'],
cfg['drive'] and the spin operator are used.
"""
from __future__ import annotations

import json
import os

import numpy as np

from ..drives.peierls import Drive
from .floquet_spin import floquet_spin_observables


def _is_float(x):
    try:
        float(x)
        return True
    except ValueError:
        return False


def _kpath_map(cfg):
    out = {}
    for item in cfg.get("kpath") or []:
        try:
            lab, k = item
        except Exception:                                     # noqa: BLE001
            continue
        out[str(lab)] = np.asarray(k, float).reshape(3)
    return out


def resolve_points(cfg, spec=None):
    """Resolve point labels or explicit triples into [(label, k_frac)]."""
    kmap = _kpath_map(cfg)
    alias = {"gamma": "Γ", "g": "Γ", "Γ": "Γ", "k": "K", "kp": "K′",
             "k'": "K′", "kprime": "K′", "m": "M"}
    if not spec:
        spec = []
        for lab in ("Γ", "K", "K′"):
            if lab in kmap:
                spec.append(lab)
    if isinstance(spec, str):
        parts = [x.strip() for x in spec.split(",")]
        if len(parts) == 3 and all(_is_float(x) for x in parts):
            spec = [spec]
        elif ";" in spec:
            spec = [s.strip() for s in spec.split(";") if s.strip()]
        else:
            spec = [s.strip() for s in spec.split(",") if s.strip()]
    points = []
    for item in spec:
        if isinstance(item, (list, tuple, np.ndarray)) and len(item) == 3:
            points.append((str(len(points)), np.asarray(item, float)))
            continue
        s = str(item)
        key = alias.get(s.lower(), s)
        if key == "K′" and key not in kmap and "K" in kmap:
            # TRS partner of K: -K modulo a reciprocal lattice vector
            points.append((key, (-kmap["K"]) % 1.0))
        elif key in kmap:
            points.append((key, kmap[key]))
        else:
            parts = s.replace(";", ",").split(",")
            if len(parts) == 3:
                points.append((s, np.asarray([float(x) for x in parts])))
            else:
                raise ValueError(f"cannot resolve k point {item!r}; "
                                 f"known labels: {sorted(kmap)}")
    return points


def make_drive(hw, gamma, pol, a_ref):
    p = pol.strip().lower()
    if p in ("circ+", "sigma+", "s+", "circular+"):
        return Drive(hw=hw, gamma=gamma, a_ref=a_ref,
                     polarization="circular", sign=+1)
    if p in ("circ-", "sigma-", "s-", "circular-"):
        return Drive(hw=hw, gamma=gamma, a_ref=a_ref,
                     polarization="circular", sign=-1)
    if p in ("linear", "lin"):
        return Drive(hw=hw, gamma=gamma, a_ref=a_ref, polarization="linear")
    raise ValueError(f"unknown polarization {pol!r}")


def target_bands(cfg, model):
    """Gap-edge bands **and their names**, from a single source of truth.

    Returns ``(bands, names)`` with ``names`` a dict ``band -> label``.

    * ``num_electrons`` present: the canonical VBM-1/VBM/CBM/CBM+1 quadruple,
      clipped to the available band range (so a 2-band model yields just the
      valence pair and *keeps the VBM-1/VBM names*, it never renames band 0).
    * ``num_electrons`` absent: the lowest four bands, honestly labelled
      ``band0..band3`` instead of pretending they are gap edges.
    """
    ne = cfg.get("num_electrons")
    if not ne:
        bands = list(range(min(4, model.nw)))
        return bands, {b: f"band{b}" for b in bands}
    edges = [int(ne) - 2, int(ne) - 1, int(ne), int(ne) + 1]
    names = {b: nm for b, nm in zip(edges, ["VBM-1", "VBM", "CBM", "CBM+1"]) if 0 <= b < model.nw}
    return [b for b in edges if 0 <= b < model.nw], names


def run_points(model, spin_at_k, cfg, points, hw, gammas, pols,
               n_step=150, m_sp=16, subspace=False, subspace_tol=5e-3):
    ne = cfg.get("num_electrons")
    bands, names = target_bands(cfg, model)
    results = {}
    match_tol = float(cfg.get("floquet_match_tol", 5e-3))
    for label, k in points:
        E, U = np.linalg.eigh(model.Hk(k))
        gaps = {}
        for b in bands:
            others = np.delete(E, b)
            gaps[b] = float(np.min(np.abs(others - E[b])))
        S = np.asarray(spin_at_k(k))
        static_spin = {}
        for b in bands:
            u = U[:, b]
            static_spin[b] = np.real([u.conj() @ S[c] @ u for c in range(3)])
        point = {}
        for gamma in gammas:
            for pol in pols:
                drive = make_drive(hw, gamma, pol,
                                   getattr(model, "a_ref", 1.0))
                obs = floquet_spin_observables(
                    model, drive, k, spin_at_k(k),
                    n_step=n_step, m_sp=m_sp)
                tab = {}
                for b in bands:
                    match = U[:, b].conj() @ obs["evec"]
                    j = int(np.argmax(np.abs(match)))
                    ov = float(np.abs(match[j]))
                    Sb = np.asarray(obs["S_bar"][j])
                    dS = Sb - static_spin[b]
                    tab[names[b]] = dict(
                        band=int(b), overlap=ov,
                        static_gap=gaps[b],
                        reliable=bool(gaps[b] >= match_tol),
                        Ebar=float(obs["Ebar"][j]),
                        eps=float(obs["eps_u"][j]),
                        w_main=float(np.max(obs["wm"][j])),
                        w_sum=float(obs["w_sum"][j]),
                        S_bar=Sb.tolist(),
                        dS_bar=dS.tolist(),
                        dS_norm=float(np.linalg.norm(dS)),
                        micromotion=float(obs["micromotion_norm"][j]))
                if subspace:
                    for nm, group in (("valence_pair", (int(ne) - 2,
                                                        int(ne) - 1)),
                                      ("conduction_pair", (int(ne),
                                                          int(ne) + 1))):
                        group = tuple(b for b in group if 0 <= b < model.nw)
                        if len(group) < 2:
                            continue
                        P = U[:, list(group)] @ U[:, list(group)].conj().T
                        w = np.array([float(np.linalg.norm(
                            P @ obs["evec"][:, j])) ** 2
                            for j in range(obs["evec"].shape[1])])
                        idx = np.argsort(w)[::-1][:len(group)]
                        ws = float(w[idx].sum())
                        if ws <= 1e-14:
                            continue
                        Sb = np.einsum("i,ic->c", w[idx] / ws,
                                       obs["S_bar"][idx])
                        tab[nm] = dict(
                            bands=list(group), n_states=int(len(group)),
                            overlap=float(np.sqrt(ws / len(group))),
                            Ebar=float(np.sum(w[idx] * obs["Ebar"][idx]) / ws),
                            S_bar=Sb.tolist(),
                            dS_norm=float(np.linalg.norm(
                                Sb - np.mean([static_spin[b] for b in group],
                                             axis=0))),
                            w_main=float(np.sum(
                                w[idx] * obs["w_main"][idx]) / ws),
                            w_sum=float(np.sum(w[idx] * obs["w_sum"][idx])
                                        / ws),
                            reliable=True, subspace=True,
                            static_gap=float(min(np.min(np.abs(
                                np.delete(E, b) - E[b])) for b in group)))
                point[f"g{gamma:g}_{pol}"] = dict(
                    targets=tab, unitarity=float(obs["unitarity_error"]),
                    hw=float(hw), gamma=float(gamma), pol=pol)
        # doublet mean-energy splitting at this point
        for key, entry in point.items():
            obs_eps = [entry["targets"][names[b]]["Ebar"] for b in bands]
            entry["doublet"] = dict(
                valence_dEbar=float(abs(obs_eps[1] - obs_eps[0]))
                if len(bands) >= 2 else None,
                conduction_dEbar=float(abs(obs_eps[3] - obs_eps[2]))
                if len(bands) >= 4 else None)
        results[label] = dict(k=k.tolist(), cases=point)
    return dict(hw=float(hw), bands=bands, band_names=names,
                points=results, n_step=int(n_step), m_sp=int(m_sp))


def _doublet_energy(obs, U, bands):
    """Mean-energy splitting of the two Floquet states matching a subspace."""
    P = U[:, list(bands)] @ U[:, list(bands)].conj().T
    w = np.array([float(np.linalg.norm(P @ obs["evec"][:, j]))
                  for j in range(obs["evec"].shape[1])])
    idx = np.argsort(w)[::-1][:len(bands)]
    idx = idx[np.argsort(obs["Ebar"][idx])]
    return float(obs["Ebar"][idx[-1]] - obs["Ebar"][idx[0]]), float(w[idx].min())


def rashba_effective(model, spin_at_k, cfg, hw, gammas,
                     direction_frac=None, tmin=5e-3, tmax=3e-2, npts=7,
                     n_step=150, m_sp=16, pol="circ+"):
    """Effective valence-doublet splitting along Gamma -> K vs gamma.

    Returns alpha_eff = slope/2 of the mean-energy splitting (not a
    quasi-energy difference), plus R^2 and the matched-subspace overlap.
    """
    ne = cfg.get("num_electrons")
    if not ne:
        raise ValueError("rashba scan needs num_electrons")
    bands = (int(ne) - 2, int(ne) - 1)
    lat = np.asarray(model.lattice, float)
    B = 2.0 * np.pi * np.linalg.inv(lat).T
    if direction_frac is None:
        direction_frac = np.array([1 / 3, 1 / 3, 0.0])
    u = (np.asarray(direction_frac, float) @ B)
    u = u / np.linalg.norm(u)
    ts = np.linspace(float(tmin), float(tmax), int(npts))
    kf = np.array([(t * u) @ np.linalg.inv(B) for t in ts])
    out = {}
    for gamma in gammas:
        drive = make_drive(hw, gamma, pol, getattr(model, "a_ref", 1.0))
        dd, ovs = [], []
        for k in kf:
            E, U = np.linalg.eigh(model.Hk(k))
            obs = floquet_spin_observables(model, drive, k, spin_at_k(k),
                                           n_step=n_step, m_sp=m_sp)
            d, ov = _doublet_energy(obs, U, bands)
            dd.append(d)
            ovs.append(ov)
        dd = np.asarray(dd)
        A = 2.0 * ts
        alpha = float(np.sum(A * dd) / np.sum(A * A))
        pred = alpha * A
        r2 = 1 - np.sum((dd - pred) ** 2) / np.sum((dd - dd.mean()) ** 2)
        out[f"{gamma:g}"] = dict(alpha_eff_ev_a=alpha, r2=float(r2),
                                 dE_meV=(dd * 1000).tolist(),
                                 t=list(map(float, ts)),
                                 min_overlap=float(min(ovs)))
    return out


def local_grid_at_hw(model, spin_at_k, cfg, center_frac, half, n, hw,
                     gamma, pol, n_step=150, m_sp=16):
    """Small 2D grid around a valley: Ebar, S_z and |S_xy| of the VBM band."""
    ne = cfg.get("num_electrons")
    band = (int(ne) - 1) if ne else (model.nw - 1)
    lat = np.asarray(model.lattice, float)
    B = 2.0 * np.pi * np.linalg.inv(lat).T
    ccart = np.asarray(center_frac, float).reshape(3) @ B
    q = np.linspace(-float(half), float(half), int(n))
    qx, qy = np.meshgrid(q, q)
    drive = make_drive(hw, gamma, pol, getattr(model, "a_ref", 1.0))
    rows = []
    for dx, dy in zip(qx.ravel(), qy.ravel()):
        kcart = np.array([ccart[0] + dx, ccart[1] + dy, 0.0])
        k = kcart @ np.linalg.inv(B)
        E, U = np.linalg.eigh(model.Hk(k))
        obs = floquet_spin_observables(model, drive, k, spin_at_k(k),
                                       n_step=n_step, m_sp=m_sp)
        match = U[:, band].conj() @ obs["evec"]
        j = int(np.argmax(np.abs(match)))
        ov = float(np.abs(match[j]))
        Sb = np.asarray(obs["S_bar"][j])
        rows.append(dict(q=[float(dx), float(dy)],
                         Ebar=float(obs["Ebar"][j]),
                         Sz=float(Sb[2]),
                         Sxy=float(np.hypot(Sb[0], Sb[1])),
                         overlap=ov))
    return rows


def valley_asymmetry(res):
    labels = {lab: res["points"][lab]["k"] for lab in res["points"]}
    out = {}
    for lab in labels:
        base = lab.lower().replace("′", "p")
        partner = None
        for other in labels:
            ob = other.lower().replace("′", "p")
            if (ob == "kp" and base == "k") or (ob == "k" and base == "kp"):
                partner = other
        if partner is None:
            continue
        for key, entry in res["points"][lab]["cases"].items():
            other = res["points"][partner]["cases"].get(key)
            if other is None:
                continue
            rows = {}
            for nm, t in entry["targets"].items():
                t2 = other["targets"][nm]
                rows[nm] = dict(
                    dEbar=float(abs(t["Ebar"] - t2["Ebar"])),
                    spin_contrast=float(np.linalg.norm(
                        np.asarray(t["S_bar"]) - np.asarray(t2["S_bar"]))),
                    dS_first=t["dS_norm"], dS_second=t2["dS_norm"])
            out.setdefault(key, {})[f"{lab}|{partner}"] = rows
    return out


def save_and_plot(res, asym, outdir):
    os.makedirs(outdir, exist_ok=True)
    json_path = os.path.join(outdir, "floquet_points.json")
    payload = dict(res)
    payload["valley_asymmetry"] = asym
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False,
                  default=lambda o: o.tolist() if hasattr(o, "tolist")
                  else str(o))
    # figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = list(res["points"])
    names = [res["band_names"][b] for b in res["bands"]]
    cases = list(next(iter(res["points"].values()))["cases"])
    gammas = [res["points"][labels[0]]["cases"][c]["gamma"] for c in cases]
    pols = [res["points"][labels[0]]["cases"][c]["pol"] for c in cases]
    # 1) band-edge spin change
    fig, axes = plt.subplots(1, len(labels), figsize=(3.2 * len(labels), 2.7),
                             constrained_layout=True)
    if len(labels) == 1:
        axes = [axes]
    for ax, lab in zip(axes, labels):
        for nm in names:
            vals = [res["points"][lab]["cases"][c]["targets"][nm]["dS_norm"]
                    for c in cases]
            ax.plot(range(len(cases)), vals, "o-", label=nm)
        ax.set_xticks(range(len(cases)))
        ax.set_xticklabels([f"{g:g}\n{p}" for g, p in zip(gammas, pols)],
                           fontsize=6)
    ax.set_ylabel("|dS|")
    axes[0].legend(frameon=False, fontsize=7)
    fig.savefig(os.path.join(outdir, "floquet_points_spin.png"),
                dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    # 2) valley asymmetry (only when both K and K' are present)
    out = dict(json=json_path,
               spin=os.path.join(outdir, "floquet_points_spin.png"))
    if asym:
        fig, ax = plt.subplots(figsize=(5.2, 2.8),
                               constrained_layout=True)
        for key, rows in asym.items():
            for nm in names:
                vals = [rows[k][nm]["dEbar"] * 1000 for k in rows]
                ax.plot(range(len(vals)), vals, "o-", label=f"{key}:{nm}")
        ax.set_ylabel("valley |dE| (meV)")
        ax.set_xticks(range(len(asym)))
        ax.set_xticklabels(list(asym), rotation=30, fontsize=7)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(frameon=False, fontsize=6)
        p = os.path.join(outdir, "floquet_points_valley.png")
        fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        out["valley"] = p
    # NB: `out` is built once above; every optional figure must be *added* to
    # it.  Re-binding `out = dict(...)` here used to drop the `valley` key that
    # the `if asym:` block had just set and to advertise a path that may never
    # have been written.
    extra = res.get("extra") or {}
    if extra.get("rashba"):
        fig, ax = plt.subplots(figsize=(3.6, 2.8),
                               constrained_layout=True)
        gs = [float(g) for g in extra["rashba"]]
        al = [extra["rashba"][g]["alpha_eff_ev_a"] * 1000 for g in
              extra["rashba"]]
        ax.plot(gs, al, "o-", color="#0072B2")
        ax.set_xlabel("gamma")
        ax.set_ylabel("alpha_eff (meV A)")
        p = os.path.join(outdir, "floquet_points_rashba.png")
        fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        out["rashba"] = p
    if extra.get("grid"):
        labels = list(extra["grid"])
        fig, axes = plt.subplots(1, len(labels), figsize=(3.0 * len(labels), 2.7),
                                 constrained_layout=True)
        if len(labels) == 1:
            axes = [axes]
        for ax, lab in zip(axes, labels):
            rows = extra["grid"][lab]
            qx = np.array([r["q"][0] for r in rows])
            qy = np.array([r["q"][1] for r in rows])
            sz = np.array([r["Sz"] for r in rows])
            cf = ax.tricontourf(qx, qy, sz,
                                levels=np.linspace(-1, 1, 41),
                                cmap="bwr", vmin=-1, vmax=1)
            ax.set_aspect("equal")
            ax.set_xticks([]); ax.set_yticks([])
            ax.text(0.03, 0.97, lab, transform=ax.transAxes,
                    fontsize=9, va="top")
        cbar = fig.colorbar(cf, ax=list(axes), fraction=0.03, pad=0.02,
                            label="S_z")
        cbar.set_ticks(np.linspace(-1, 1, 5))
        p = os.path.join(outdir, "floquet_points_grid.png")
        fig.savefig(p, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        out["grid"] = p
    return out
