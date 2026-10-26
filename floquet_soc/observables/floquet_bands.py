# -*- coding: utf-8 -*-
"""Floquet band comparison along a k path.

The physically unambiguous "after" object is the spectral function A(k,E):
quasienergies are defined modulo hbar*omega and cannot be drawn as ordinary
bands without folding/continuation.  We therefore show

  (a) static bands E(k)                      (before)
  (b) Floquet spectral weight A(k,E) with the static bands overlaid
      (after; photon replicas appear at E = eps_F + m*hbar*omega)

Quasienergies/mean energies are still stored in the NPZ for analysis, but are
not used as band labels.
"""
from __future__ import annotations

import json
import os

import numpy as np

from ..solvers.evolution import floquet_modes
from .spectral import spectral_function


def _match_columns(prev_evec, cur_evec):
    """Greedy max-overlap bijection: carried state -> current column index.

    ``floquet_modes`` returns the eigenvector columns in an order that is
    arbitrary per k (Schur).  Sorting each k by quasienergy destroys state
    identity across a band crossing; instead we follow each state by the
    overlap ``|<prev_state|cur_state>|``.  Greedy + used-column bookkeeping is
    O(n^2) and exact away from exact degeneracies (where the order is anyway
    ambiguous); ``scipy.optimize.linear_sum_assignment`` is the rigorous
    fallback for large-degeneracy systems.
    """
    ov = np.abs(prev_evec.conj().T @ cur_evec)
    nw = ov.shape[0]
    perm = np.empty(nw, dtype=int)
    taken = np.zeros(nw, dtype=bool)
    for a in range(nw):
        score = ov[a].copy()
        score[taken] = -1.0
        b = int(np.argmax(score))
        perm[a] = b
        taken[b] = True
    return perm


def compute_path(model, drive, kpts, n_step=100, m_sp=16, sigma=0.01,
                 ne=None):
    kpts = np.atleast_2d(np.asarray(kpts, float))
    nk, nw = len(kpts), int(model.nw)
    E = np.empty((nk, nw))
    for i, k in enumerate(kpts):
        E[i] = np.linalg.eigvalsh(model.Hk(k))
    ref = float(E[:, ne - 1].max()) if ne else float(E.max())
    E_grid = np.linspace(ref - 1.2, ref + 1.8, 700)
def _unwrap_along_path(eps_tracked, hw):
    """Remove +/- hw jumps, independently per column (one tracked state each).

    ``eps_tracked`` is folded ([-hw/2, hw/2)) and state-tracked across k; a
    state that drifts past the zone edge folds back, producing a +/- hw step.
    This returns the same states unfolded to a continuous value, which is the
    correct object for plotting ordinary bands or taking k differences.
    """
    eps_tracked = np.asarray(eps_tracked, float)
    out = np.empty_like(eps_tracked)
    out[0] = eps_tracked[0]
    offset = np.zeros(eps_tracked.shape[1])
    for i in range(1, len(eps_tracked)):
        offset -= np.rint((eps_tracked[i] - out[i - 1]) / hw)
        out[i] = eps_tracked[i] + offset * hw
    return out


def compute_path(model, drive, kpts, n_step=100, m_sp=16, sigma=0.01,
                 ne=None):
    kpts = np.atleast_2d(np.asarray(kpts, float))
    nk, nw = len(kpts), int(model.nw)
    E = np.empty((nk, nw))
    for i, k in enumerate(kpts):
        E[i] = np.linalg.eigvalsh(model.Hk(k))
    ref = float(E[:, ne - 1].max()) if ne else float(E.max())
    E_grid = np.linspace(ref - 1.2, ref + 1.8, 700)
    eps = np.empty((nk, nw))
    Ebar = np.empty((nk, nw))
    wsum = np.empty((nk, nw))
    A = np.empty((nk, len(E_grid)))
    A0 = np.empty((nk, len(E_grid)))
    prev_evec = None
    for i, k in enumerate(kpts):
        fl = floquet_modes(model, drive, k, n_step=n_step, want_wspec=True,
                           m_sp=m_sp)
        evec = np.asarray(fl["evec"])
        eps_u = np.asarray(fl["eps_u"])
        if i == 0:
            perm = np.argsort(eps_u)
        else:
            perm = _match_columns(prev_evec, evec)
        eps[i] = eps_u[perm]                      # folded + state-tracked
        Ebar[i] = np.asarray(fl["Ebar"])[perm]
        wsum[i] = np.asarray(fl["wspec"]["w_sum"])[perm]
        prev_evec = evec[:, perm]
        A[i] = spectral_function(fl["eps_u"], fl["wspec"]["ms"],
                                 fl["wspec"]["wm"], E_grid, sigma,
                                 drive.hw)
        ms = np.asarray(fl["wspec"]["ms"])
        im0 = int(np.where(ms == 0)[0][0])
        A0[i] = spectral_function(fl["eps_u"], [0],
                                  fl["wspec"]["wm"][:, im0:im0 + 1],
                                  E_grid, sigma, drive.hw)
    eps_unwrapped = _unwrap_along_path(eps, drive.hw)
    return dict(kpts=kpts, E=E, eps=eps, eps_unwrapped=eps_unwrapped,
                Ebar=Ebar, w_sum=wsum, A=A, A0=A0, E_grid=E_grid, sigma=sigma,
                hw=float(drive.hw), gamma=float(drive.gamma),
                polarization=drive.polarization)


def save_and_plot(res, labels, outdir, prefix="floquet_bands"):
    os.makedirs(outdir, exist_ok=True)
    npz = os.path.join(outdir, prefix + ".npz")
    np.savez(npz, **res)
    with open(os.path.join(outdir, prefix + ".json"), "w",
              encoding="utf-8") as f:
        json.dump(dict(hw=res["hw"], gamma=res["gamma"],
                       polarization=res["polarization"], labels=labels,
                       sigma=res["sigma"], n_k=len(res["kpts"])), f,
                  indent=2, ensure_ascii=False)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .. import plot_style as ps
    ps.enforce_style()
    E, A, Eg = res["E"], res["A0"], res["E_grid"]
    ne = res.get("ne")
    ref = float(E[:, ne - 1].max()) if ne else float(E.max())
    nk, nw = E.shape
    x = np.arange(nk)
    ylo, yhi = -2.0, 3.0
    sel = (Eg - ref >= ylo) & (Eg - ref <= yhi)
    Egw = Eg[sel] - ref
    Aw = A[:, sel]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), sharex=True,
                             sharey=True, constrained_layout=True)
    ax0, ax1 = axes
    # (a) all static bands (clipped by the common y window)
    for b in range(nw):
        ax0.plot(x, E[:, b] - ref, color="#777777", lw=0.8)
    ax0.set_ylabel("$E-E_{\\rm VBM}$ (eV)")
    ax0.set_xlabel("k path")
    # (b) spectral function + all static bands
    vmax = float(np.percentile(Aw, 99.5)) if Aw.max() > 0 else 1.0
    cf = ax1.pcolormesh(x, Egw, Aw.T, shading="auto", cmap="inferno",
                        vmin=0.0, vmax=vmax)
    for b in range(nw):
        ax1.plot(x, E[:, b] - ref, color="#9ecae1", lw=0.6, alpha=0.7)
    ax1.set_xlabel("k path")
    cb = fig.colorbar(cf, ax=ax1, fraction=0.05, pad=0.02)
    cb.set_label("A0(k,E)  (m=0)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    # high-symmetry positions: exactly the segment boundaries, no margins
    pos = []
    if labels:
        nseg = max(1, (nk - 1) // max(1, len(labels) - 1))
        pos = [min(i * nseg, nk - 1) for i in range(len(labels))]
    else:
        pos = [0, nk - 1]
    for ax in axes:
        ax.set_xlim(0, nk - 1)
        ax.set_ylim(ylo, yhi)
        for p in pos[1:-1]:
            ax.axvline(p, color="#999999", lw=0.6, ls="--", zorder=0)
        ax.set_xticks(pos)
        ax.set_xticklabels(labels if labels else [""] * len(pos))
    fig.savefig(os.path.join(outdir, prefix + ".png"), dpi=300,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return dict(npz=npz, figure=os.path.join(outdir, prefix + ".png"))
