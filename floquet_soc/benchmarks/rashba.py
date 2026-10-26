# -*- coding: utf-8 -*-
"""Benchmark 2: 2-band Rashba model.

    H(k) = (hbar^2 k^2 / 2m) I + alpha_R (k_y sigma_x - k_x sigma_y)

Analytic expectations:
  * band splitting  Delta E(k) = 2 alpha_R |k|;
  * spin-momentum locking: <S> = +- (1/2) (k_y, -k_x, 0)/|k|, so S_z = 0,
    S . k = 0, |S| = 1/2 (hbar units);
  * PST metric of one band over a full circle: eigenvalues (1/2, 1/2, 0).
"""
from __future__ import annotations

import numpy as np

from .. import plot_style as ps

from ..core.spin_operator import SPIN_HALF, SIGMA
from ..core.validation import (ValidationReport, PASS, FAIL, WARN)
from ..observables.spin_texture import spin_texture
from ..observables.pst import pst_collinearity


class RashbaModel:
    kind = "analytic"
    name = "Rashba"

    def __init__(self, hbar2_over_2m=0.5, alpha=0.1):
        self.c = float(hbar2_over_2m)
        self.alpha = float(alpha)
        self.nw = 2

    def Hk(self, k):
        k = np.asarray(k, float).reshape(3)
        kx, ky = k[0], k[1]
        H = self.c * (kx * kx + ky * ky) * np.eye(2, dtype=complex)
        H += self.alpha * (ky * SIGMA[0] - kx * SIGMA[1])
        return H

    def Hk_batch(self, kpts):
        kpts = np.atleast_2d(np.asarray(kpts, float))
        return np.array([self.Hk(k) for k in kpts])

    def static_bands(self, kpts):
        return np.linalg.eigvalsh(self.Hk_batch(kpts))

    def spin_ops(self):
        return SPIN_HALF

    def splitting(self, k):
        ev = np.linalg.eigvalsh(self.Hk(k))
        return float(ev[1] - ev[0])

    def spin_at(self, k):
        return spin_texture(self.Hk(k), SPIN_HALF)


def run(outdir=None):
    model = RashbaModel()
    rep = ValidationReport("Benchmark 2: Rashba SOC model")
    ks = np.linspace(0.02, 0.30, 57)
    alpha, c = model.alpha, model.c

    # 1) splitting = 2 alpha k
    err = 0.0
    for k in ks:
        err = max(err, abs(model.splitting([k, 0, 0]) - 2 * alpha * k))
    rep.add("rashba_splitting_analytic", PASS if err < 1e-10 else FAIL,
            f"max |DeltaE - 2 alpha k| = {err:.2e} eV (tol 1e-10)")

    # 2) spin-momentum locking, |S|=1/2, Sz=0, S.k=0
    dev_lock = dev_sz = dev_norm = 0.0
    for k in ks:
        for ang in np.linspace(0, 2 * np.pi, 13, endpoint=False):
            kk = np.array([k * np.cos(ang), k * np.sin(ang), 0.0])
            res = model.spin_at(kk)
            for ib in range(2):
                s = res["spins"][ib]
                s_an = 0.5 * np.array([kk[1], -kk[0], 0.0]) / k
                sign = 1.0 if np.dot(s, s_an) >= 0 else -1.0
                dev_lock = max(dev_lock, float(np.linalg.norm(s - sign * s_an)))
                dev_sz = max(dev_sz, abs(float(s[2])))
                dev_norm = max(dev_norm, abs(float(np.linalg.norm(s)) - 0.5))
    rep.add("rashba_spin_momentum_locking", PASS if dev_lock < 1e-10 else FAIL,
            f"max |S - S_analytic| = {dev_lock:.2e} (tol 1e-10)")
    rep.add("rashba_Sz_zero", PASS if dev_sz < 1e-12 else FAIL,
            f"max |S_z| = {dev_sz:.2e}")
    rep.add("rashba_spin_norm_half", PASS if dev_norm < 1e-10 else FAIL,
            f"max ||S| - 1/2| = {dev_norm:.2e} (hbar units)")

    # 3) PST metric of a single band over a circle: (1/2,1/2,0)
    circ = np.array([[0.2 * np.cos(a), 0.2 * np.sin(a), 0.0]
                     for a in np.linspace(0, 2 * np.pi, 37, endpoint=False)])
    spins = np.array([model.spin_at(k)["spins"][0] for k in circ])
    pst = pst_collinearity(spins)
    expected = np.array([0.5, 0.5, 0.0])
    dev = float(np.max(np.abs(pst["eigenvalues"] - expected)))
    rep.add("rashba_pst_metric_inplane", PASS if dev < 1e-8 else FAIL,
            f"M eigenvalues = {np.round(pst['eigenvalues'], 6)}, "
            f"expected (0.5,0.5,0) (max dev {dev:.1e})")

    # 4) Hermiticity
    herm = max(float(np.max(np.abs(model.Hk(k) - model.Hk(k).conj().T)))
               for k in circ)
    rep.add("rashba_hermiticity", PASS if herm < 1e-14 else FAIL,
            f"max|H-H^dag| = {herm:.1e}")

    fig = None
    if outdir:
        fig = _plot(model, outdir)
    return dict(model=model, report=rep, figure=fig, pst=pst)


def _plot(model, outdir):
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    ks = np.linspace(0, 0.30, 121)
    ps.enforce_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4),
                             constrained_layout=True)
    ax = axes[0]
    for ib, col in ((0, "black"), (1, "#c0392b")):
        e = [model.spin_at([k, 0, 0])["energies"][ib] for k in ks]
        ax.plot(ks, e, color=col, lw=1.6,
                label=f"band {ib} (analytic)")
    ax.set_xlabel("$k_x$ (1/Angstrom)")
    ax.set_ylabel("E (eV)")
    ax.legend(fontsize=8, frameon=False)
    ax = axes[1]
    X, Y = np.meshgrid(np.linspace(-0.3, 0.3, 15), np.linspace(-0.3, 0.3, 15))
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            s = model.spin_at([X[i, j], Y[i, j], 0.0])["spins"][0]
            U[i, j], V[i, j] = s[0], s[1]
    ax.quiver(X, Y, U, V, color="black", pivot="middle", width=0.004)
    ax.set_xlabel("$k_x$ (1/Angstrom)")
    ax.set_aspect("equal")
    out = os.path.join(outdir, "benchmark_rashba.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out
