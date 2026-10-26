# -*- coding: utf-8 -*-
"""Benchmark 3: Rashba + Dresselhaus 2-band model.

    H = (hbar^2 k^2/2m) I
        + alpha_R (k_y sigma_x - k_x sigma_y)
        + beta_D  (k_x sigma_x - k_y sigma_y)

Analytic expectations:
  * spin-orbit field d(k) = (alpha ky + beta kx, -alpha kx - beta ky, 0),
    band splitting Delta E = 2 |d|;
  * at alpha = beta the spin direction is *momentum independent*:
    <S> = +- (1/2) (1,-1,0)/sqrt(2)  ->  persistent spin texture (PST);
    PST metric eigenvalues (1,0,0), principal axis (1,-1,0)/sqrt(2);
  * at alpha != beta the texture is non-collinear (collinearity < 1).
"""
from __future__ import annotations

import numpy as np

from .. import plot_style as ps

from ..core.spin_operator import SPIN_HALF, SIGMA
from ..core.validation import ValidationReport, PASS, FAIL
from ..observables.spin_texture import spin_texture
from ..observables.pst import pst_collinearity, axis_angle


class RashbaDresselhausModel:
    kind = "analytic"
    name = "Rashba+Dresselhaus"

    def __init__(self, hbar2_over_2m=0.5, alpha=0.1, beta=0.1):
        self.c = float(hbar2_over_2m)
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.nw = 2

    def d_vector(self, k):
        kx, ky = k[0], k[1]
        return np.array([self.alpha * ky + self.beta * kx,
                         -self.alpha * kx - self.beta * ky, 0.0])

    def Hk(self, k):
        k = np.asarray(k, float).reshape(3)
        kx, ky = k[0], k[1]
        d = self.d_vector(k)
        H = self.c * (kx * kx + ky * ky) * np.eye(2, dtype=complex)
        H += d[0] * SIGMA[0] + d[1] * SIGMA[1]
        return H

    def Hk_batch(self, kpts):
        kpts = np.atleast_2d(np.asarray(kpts, float))
        return np.array([self.Hk(k) for k in kpts])

    def static_bands(self, kpts):
        return np.linalg.eigvalsh(self.Hk_batch(kpts))

    def spin_ops(self):
        return SPIN_HALF

    def spin_at(self, k):
        return spin_texture(self.Hk(k), SPIN_HALF)


def run(outdir=None):
    rep = ValidationReport("Benchmark 3: Rashba+Dresselhaus, PST")
    alpha = beta = 0.1
    model = RashbaDresselhausModel(alpha=alpha, beta=beta)
    ks = np.linspace(0.02, 0.30, 57)

    # 1) splitting = 2|d|
    err = 0.0
    for k in ks:
        for ang in np.linspace(0, 2 * np.pi, 13, endpoint=False):
            kk = np.array([k * np.cos(ang), k * np.sin(ang), 0.0])
            ev = np.linalg.eigvalsh(model.Hk(kk))
            err = max(err, abs(float(ev[1] - ev[0]) - 2 * np.linalg.norm(
                model.d_vector(kk))))
    rep.add("rd_splitting_analytic", PASS if err < 1e-10 else FAIL,
            f"max |DeltaE - 2|d(k)|| = {err:.2e} eV")

    # 2) alpha=beta: collinear spin direction (1,-1,0)/sqrt2
    axis_expect = np.array([1.0, -1.0, 0.0]) / np.sqrt(2)
    dev = 0.0
    for k in ks:
        for ang in np.linspace(0, 2 * np.pi, 13, endpoint=False):
            kk = np.array([k * np.cos(ang), k * np.sin(ang), 0.0])
            for ib in range(2):
                s = model.spin_at(kk)["spins"][ib]
                # H = d.sigma: <sigma> = +-d/|d|; at alpha=beta the spin
                # direction is +-(1,-1,0)/sqrt2 with the sign set by
                # sign(kx+ky).  Collinearity means S is parallel to the axis
                # up to an overall sign, so compare with both orientations.
                s_an = 0.5 * axis_expect
                dev = max(dev, float(min(np.linalg.norm(s - s_an),
                                         np.linalg.norm(s + s_an))))
    rep.add("rd_alpha_eq_beta_collinear", PASS if dev < 1e-10 else FAIL,
            f"max |S - +-(1/2)(1,-1,0)/sqrt2| = {dev:.2e}")

    # 3) PST metric at alpha=beta over full circle
    circ = np.array([[0.2 * np.cos(a), 0.2 * np.sin(a), 0.0]
                     for a in np.linspace(0, 2 * np.pi, 49, endpoint=False)])
    spins = np.array([model.spin_at(k)["spins"][0] for k in circ])
    pst = pst_collinearity(spins)
    dev_ev = float(np.max(np.abs(pst["eigenvalues"] - np.array([1.0, 0, 0]))))
    dev_ax = axis_angle(pst["pst_axis"], axis_expect)
    rep.add("rd_pst_collinearity_one", PASS if dev_ev < 1e-8 else FAIL,
            f"M eigenvalues = {np.round(pst['eigenvalues'], 8)} -> (1,0,0)")
    rep.add("rd_pst_axis", PASS if dev_ax < 1e-8 else FAIL,
            f"PST axis = {np.round(pst['pst_axis'], 8)}, "
            f"angle to (1,-1,0)/sqrt2 = {np.degrees(dev_ax):.2e} deg")

    # 4) alpha != beta -> not collinear
    model2 = RashbaDresselhausModel(alpha=0.1, beta=0.03)
    spins2 = np.array([model2.spin_at(k)["spins"][0] for k in circ])
    pst2 = pst_collinearity(spins2)
    rep.add("rd_alpha_ne_beta_noncollinear", PASS if pst2["collinearity"] < 0.99
            else FAIL,
            f"collinearity(alpha!=beta) = {pst2['collinearity']:.4f} < 1")

    fig = None
    if outdir:
        fig = _plot(model, model2, outdir)
    return dict(model=model, report=rep, figure=fig, pst=pst, pst2=pst2)


def _plot(model, model2, outdir):
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    ps.enforce_style()
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3),
                             constrained_layout=True)
    ks = np.linspace(0, 0.3, 121)
    ax = axes[0]
    for ib, col in ((0, "black"), (1, "#c0392b")):
        ax.plot(ks, [model.spin_at([k, 0, 0])["energies"][ib] for k in ks],
                color=col, lw=1.6)
    ax.set_xlabel("$k_x$ (1/Angstrom)")
    ax.set_ylabel("E (eV)")
    for ax, mod, title in ((axes[1], model, "(b) spin texture alpha=beta (PST)"),
                           (axes[2], model2, "(c) alpha!=beta (non-collinear)")):
        X, Y = np.meshgrid(np.linspace(-0.3, 0.3, 15),
                           np.linspace(-0.3, 0.3, 15))
        U = np.zeros_like(X)
        V = np.zeros_like(Y)
        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                s = mod.spin_at([X[i, j], Y[i, j], 0.0])["spins"][0]
                U[i, j], V[i, j] = s[0], s[1]
        ax.quiver(X, Y, U, V, color="black", pivot="middle", width=0.004)
        ax.set_xlabel("$k_x$ (1/Angstrom)")
        ax.set_aspect("equal")
    out = os.path.join(outdir, "benchmark_rashba_dresselhaus.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out
