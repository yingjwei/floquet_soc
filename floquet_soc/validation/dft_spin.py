# -*- coding: utf-8 -*-
"""Compare a model spin texture with a DFT spin-texture table.

Generic for any material: the DFT file only needs the documented columns

    #kx ky Energy Sx Sy Sz
    #Grid-Size: nx ny
    #Band-Index: n      (1-based)

The model side uses the same (wsvec-aware) H(k) and spin operator as
``floquet.py spin_texture``.  Reported metrics:

* energy offset / residual after removing a constant reference shift;
* raw spin-component RMS/max and correlations;
* global scale fit (DFT files often carry a non-Pauli normalization);
* normalized-direction RMS (gauge/scale invariant), overall and for
  non-degenerate points;
* one representative valley (K) direction comparison.

Single-band spin is gauge dependent at degeneracies; points with a neighbour
gap below ``mask_split`` are excluded from the direction claim.
"""
from __future__ import annotations

import json
import os

import numpy as np

from ..observables.spin_plot import (plot_spin_texture_2d,
                                     plot_spin_texture_panels)


def _corr_or_nan(a, b) -> float:
    """Pearson correlation; nan when either series is constant (0/0)."""
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def read_spin_texture_dat(path: str) -> dict:
    rows = []
    band = None
    grid = None
    with open(path, encoding="utf-8", errors="ignore") as fh:
        lines = fh.readlines()
    for line in lines:
        if line.strip().startswith("#"):
            low = line.lower()
            if "band-index" in low:
                band = int(line.split(":")[1])
            if "grid-size" in low:
                grid = [int(x) for x in line.split(":")[1].split()]
            continue
        q = line.split()
        if len(q) >= 6:
            rows.append([float(x) for x in q[:6]])
    if not rows:
        raise ValueError(f"{path}: no data rows")
    data = np.asarray(rows, float)
    return dict(kx=data[:, 0], ky=data[:, 1], energy=data[:, 2],
                spin=data[:, 3:6], band_index=band, grid=grid,
                n_points=len(data), path=path)


def compare(model, spin_at_k, dft: dict, band: int | None = None,
            mask_split: float = 0.05, chunk: int = 64) -> dict:
    kxy = np.column_stack([dft["kx"], dft["ky"]])
    lat = np.asarray(model.lattice, float)
    B = 2.0 * np.pi * np.linalg.inv(lat).T
    kfrac = np.column_stack([kxy, np.zeros(len(kxy))]) @ np.linalg.inv(B)
    nw = int(model.nw)
    E = np.empty((len(kfrac), nw))
    S = np.empty((len(kfrac), nw, 3))
    for i0 in range(0, len(kfrac), max(1, chunk)):
        for j, k in enumerate(kfrac[i0:i0 + chunk]):
            w, U = np.linalg.eigh(model.Hk(k))
            ops = np.asarray(spin_at_k(k))
            for b in range(nw):
                u = U[:, b]
                S[i0 + j, b] = np.real(
                    [u.conj() @ ops[c] @ u for c in range(3)])
            E[i0 + j] = w
    if band is None:
        band = (dft["band_index"] - 1) if dft["band_index"] else (nw - 1)
    band = int(band)
    Ew, Sw = E[:, band], S[:, band]
    Sd, Ed = dft["spin"], dft["energy"]
    # gap to the *nearest* other band (Kramers partner for the valence edge)
    others = np.delete(E, band, axis=1)
    split = np.min(np.abs(others - E[:, band][:, None]), axis=1)
    off = float(np.mean(Ed - Ew))
    resid = Ed - Ew - off
    # per-component scale: undefined (nan) when the DFT table has no weight in
    # that component -- guard the 0/0 instead of letting numpy print a warning
    scale = []
    for c in range(3):
        denom = float(np.sum(Sd[:, c] ** 2))
        num = float(np.sum(Sw[:, c] * Sd[:, c]))
        scale.append(num / denom if denom > 0 else float("nan"))
    denom_all = float(np.sum(Sd * Sd))
    scale_all = float(np.sum(Sw * Sd) / denom_all) if denom_all > 0 else \
        float("nan")
    nwv = Sw / np.linalg.norm(Sw, axis=1)[:, None]
    ndv = Sd / np.linalg.norm(Sd, axis=1)[:, None]
    sel = np.abs(split) > float(mask_split)
    dir_all = float(np.sqrt(np.mean(np.sum((nwv - ndv) ** 2, axis=1))))
    dir_sel = float(np.sqrt(np.mean(
        np.sum((nwv[sel] - ndv[sel]) ** 2, axis=1)))) if sel.any() else None
    out = dict(
        dft_file=dft["path"], band_wannier=band,
        band_dft=dft["band_index"], grid=dft["grid"],
        n_points=int(len(kxy)),
        energy=dict(offset_ev=off,
                    rms_resid_ev=float(np.sqrt(np.mean(resid ** 2))),
                    max_abs_resid_ev=float(np.max(np.abs(resid)))),
        spin_raw=dict(rms=[float(x) for x in np.sqrt(np.mean(
                        (Sw - Sd) ** 2, axis=0))],
                      max=[float(x) for x in np.max(np.abs(Sw - Sd),
                                                    axis=0)],
                      corr=[_corr_or_nan(Sw[:, c], Sd[:, c])
                            for c in range(3)]),
        scale=dict(component=scale, global_fit=scale_all,
                   norm_w=float(np.linalg.norm(Sw, axis=1).mean()),
                   norm_dft=float(np.linalg.norm(Sd, axis=1).mean()),
                   norm_w_max=float(np.linalg.norm(Sw, axis=1).max()),
                   norm_dft_max=float(np.linalg.norm(Sd, axis=1).max())),
        direction=dict(normalized_rms_all=dir_all,
                       normalized_rms_non_degenerate=dir_sel,
                       mask_split_ev=float(mask_split),
                       n_non_degenerate=int(sel.sum()),
                       n_dft_Sz_exact_zero=int(np.sum(
                           np.abs(Sd[:, 2]) < 1e-9))),
        kxy=kxy, E=E, S=S, split=split,
        dft_spin=Sd, dft_energy=Ed,
    )
    return out


def save_and_plot(res: dict, model, spin_at_k, outdir: str,
                  prefix: str = "dft_spin_compare",
                  plots: bool = True) -> dict:
    os.makedirs(outdir, exist_ok=True)
    band = res["band_wannier"]
    kxy = res["kxy"]
    S = res["S"]
    Sd = res["dft_spin"]
    Ed = res["dft_energy"]
    Ew = res["E"][:, band]
    npz = os.path.join(outdir, prefix + ".npz")
    np.savez(npz, kxy=kxy, band=band, E=res["E"], S=S,
             energies_dft=Ed, spins_dft=Sd)
    out = dict(json=os.path.join(outdir, prefix + ".json"), npz=npz)
    if plots:
        tex = os.path.join(outdir, prefix + "_texture.png")
        plot_spin_texture_panels(
            [dict(kpts=kxy, spins=Sd), dict(kpts=kxy, spins=S[:, band, :])],
            tex)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(6.4, 2.9),
                                 constrained_layout=True)
        sel = np.abs(res["split"]) > res["direction"]["mask_split_ev"]
        ax = axes[0]
        ax.scatter(Sd[sel, 2], S[:, band, :][sel, 2], s=4,
                   c=np.linalg.norm(kxy[sel], axis=1), cmap="viridis",
                   rasterized=True)
        xx = np.linspace(Sd[sel, 2].min(), Sd[sel, 2].max(), 10)
        for c, col in ((0, "#0072B2"), (1, "#D55E00"), (2, "#333333")):
            ax.plot(xx, res["scale"]["component"][c] * xx, color=col, lw=1.0,
                    label="$S_%s$: x%.3f" % ("xyz"[c],
                                             res["scale"]["component"][c]))
        ax.set_xlabel("DFT $S_z$")
        ax.set_ylabel("Wannier $S_z$")
        ax.legend(frameon=False, fontsize=7)
        ax = axes[1]
        resid = (Ed - Ew - res["energy"]["offset_ev"]) * 1000
        ax.hist(resid, bins=40, color="#0072B2", alpha=0.8)
        ax.set_xlabel("$E_{\\rm DFT}-E_{\\rm Wannier}-\\langle\\Delta E\\rangle$ (meV)")
        ax.set_ylabel("count")
        scatter = os.path.join(outdir, prefix + "_scatter.png")
        fig.savefig(scatter, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        out["texture"] = tex
        out["scatter"] = scatter
    meta = {k: v for k, v in res.items()
            if k not in ("kxy", "E", "S", "split", "dft_spin",
                         "dft_energy")}
    meta["outputs"] = out
    with open(out["json"], "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False,
                  default=lambda o: o.tolist() if hasattr(o, "tolist")
                  else str(o))
    return out
