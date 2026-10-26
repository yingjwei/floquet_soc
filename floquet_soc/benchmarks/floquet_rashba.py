# -*- coding: utf-8 -*-
"""Stage-B benchmark: driven 2-band Rashba model.

Validates the Floquet machinery (stage B) on an analytically controllable
model before touching PbNV:

  B1  zero drive: quasienergies reproduce static bands mod hbar*omega;
      periodic-averaged spin texture reproduces the static spin texture;
  B2  photon-sector weight normalization  sum_m w_{alpha m} = 1;
  B3  Sambe Hamiltonian and direct time evolution agree;
  B4  high-frequency effective Hamiltonian converges to the exact
      quasienergies as omega grows (error decreases with omega);
  B5  spin-resolved spectral density: the signed weights satisfy
      sum_m s_{alpha m,i} = S_bar_{alpha,i} and the broadened spectrum
      integrates back to the spin expectation;
  B6  micromotion matters: report |S_bar - S_t0| (t=0 approximation error).

The drive for continuum models is minimal coupling H(k + A(t)) (length gauge),
with gamma = A0 * a_ref and a_ref = 1 Angstrom, i.e. A0 = gamma 1/Angstrom.
"""
from __future__ import annotations

import numpy as np

from .. import plot_style as ps

from ..core.spin_operator import SPIN_HALF
from ..core.validation import ValidationReport, PASS, FAIL
from ..drives.peierls import Drive
from ..solvers.evolution import floquet_modes
from ..solvers.sambe import sambe_vs_evolution
from ..solvers.effective import effective_split_error
from ..observables.floquet_spin import (floquet_spin_observables,
                                        spin_spectral_density)
from ..observables.spin_texture import spin_texture
from .rashba import RashbaModel


def _drive(hw, gamma, pol="circular", sign=-1):
    return Drive(hw=hw, gamma=gamma, a_ref=1.0, polarization=pol,
                 e1=(1, 0, 0), e2=(0, 1, 0), sign=sign)


def run(outdir=None):
    rep = ValidationReport("Stage-B benchmark: driven Rashba model")
    model = RashbaModel(hbar2_over_2m=0.5, alpha=0.1)
    ops = SPIN_HALF
    k = np.array([0.10, 0.05, 0.0])
    hw = 1.0

    # ---------------- B1: zero drive limits
    d0 = _drive(hw, 0.0)
    fl0 = floquet_modes(model, d0, k, n_step=200, want_wspec=True, m_sp=12)
    stat = np.linalg.eigvalsh(model.Hk(k))
    stat_fold = (stat + hw / 2) % hw - hw / 2
    dev_q = float(np.max(np.abs(np.sort(fl0["eps"]) - np.sort(stat_fold))))
    rep.add("B1_zero_drive_quasienergy", PASS if dev_q < 1e-10 else FAIL,
            f"max|eps - E_static mod hw| = {dev_q:.2e} eV (tol 1e-10)")

    obs0 = floquet_spin_observables(model, d0, k, ops, n_step=200, m_sp=12)
    stat_spin = spin_texture(model.Hk(k), ops)["spins"]
    dev_s = float(np.max(np.abs(obs0["S_bar"] - stat_spin)))
    rep.add("B1_zero_drive_spin_texture", PASS if dev_s < 1e-10 else FAIL,
            f"max|S_bar - S_static| = {dev_s:.2e} (tol 1e-10)")

    # ---------------- B2: photon weight normalization
    d2 = _drive(hw, 0.05)
    fl2 = floquet_modes(model, d2, k, n_step=300, want_wspec=True, m_sp=12)
    ws = fl2["wspec"]
    dev_w = float(np.max(np.abs(ws["w_sum"] - 1.0)))
    rep.add("B2_photon_weight_normalization", PASS if dev_w < 1e-8 else FAIL,
            f"max|sum_m w - 1| = {dev_w:.2e} (m_sp=12)")

    # ---------------- B3: Sambe vs evolution
    res = sambe_vs_evolution(model, d2, k, n_harm=3, n_step_sambe=400,
                             n_step_evo=300)
    dev_sam = res["max_mismatch"]
    rep.add("B3_sambe_vs_evolution", PASS if dev_sam < 1e-5 else FAIL,
            f"max quasienergy mismatch Sambe vs evolution = {dev_sam:.2e} eV "
            f"(n_harm=3, tol 1e-5)")

    # ---------------- B4: high-frequency expansion convergence
    kk = np.array([0.12, 0.04, 0.0])
    errs = {}
    for om in (1.0, 2.0, 4.0, 8.0):
        dd = _drive(om, 0.01)
        eps = floquet_modes(model, dd, kk, n_step=400)["eps"]
        r = effective_split_error(model, dd, kk, eps, order=2, n_harm=3,
                                  n_step=800, check_validity=False)
        errs[om] = r["max_abs"]
    om_sorted = sorted(errs)
    decreasing = all(errs[om_sorted[i + 1]] < errs[om_sorted[i]]
                     for i in range(len(om_sorted) - 1))
    ratio = errs[om_sorted[-1]] / max(errs[om_sorted[0]], 1e-30)
    rep.add("B4_effective_convergence",
            PASS if (decreasing and errs[8.0] < 1e-4) else FAIL,
            f"max|eps - eig(H_eff)| for w=1,2,4,8 eV: "
            f"{[float(f'{errs[o]:.2e}') for o in om_sorted]}; decreasing={decreasing}, "
            f"ratio(8/1)={ratio:.2e} (tol: final <1e-4 eV)")

    # ---------------- B5: spin-resolved spectral density identity
    d5 = _drive(1.0, 0.05)
    obs = floquet_spin_observables(model, d5, k, ops, n_step=300, m_sp=12)
    ident = float(np.max(np.abs(obs["spin_weights"].sum(axis=0)
                                - obs["S_bar"])))
    rep.add("B5_spin_weight_identity", PASS if ident < 1e-12 else FAIL,
            f"max|sum_m s_am,i - S_bar| = {ident:.2e}")
    E = np.linspace(-8.0, 8.0, 16001)
    spec = spin_spectral_density(obs, E, sigma=2e-3)
    dE = E[1] - E[0]
    integ = spec["A_signed"].sum(axis=1) * dE
    target = obs["S_bar"].sum(axis=0)
    dev_int = float(np.max(np.abs(integ - target)))
    rep.add("B5_spin_spectrum_integrates_to_spin",
            PASS if dev_int < 2e-4 else FAIL,
            f"max|∫A^S dE - sum_alpha S_bar| = {dev_int:.2e} "
            f"(sigma=2 meV, grid +-8 eV)")

    # ---------------- B6: micromotion contribution
    d6 = _drive(1.0, 0.10)
    obs6 = floquet_spin_observables(model, d6, k, ops, n_step=300, m_sp=12)
    mm = float(np.max(np.abs(obs6["S_bar"] - obs6["S_t0"])))
    rep.add("B6_micromotion_reported", PASS if np.isfinite(mm) else FAIL,
            f"max|S_bar - S_t0| = {mm:.3e} (difference between the correct "
            f"periodic average and the t=0 approximation)")

    fig = None
    if outdir:
        fig = _plot(model, ops, outdir)
    return dict(model=model, report=rep, figure=fig, errors=errs,
                micromotion=mm)


def _plot(model, ops, outdir):
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    k = np.array([0.10, 0.05, 0.0])
    d0 = _drive(1.0, 0.0)
    d1 = _drive(1.0, 0.10)
    obs = floquet_spin_observables(model, d1, k, ops, n_step=300, m_sp=12)
    stat = spin_texture(model.Hk(k), ops)["spins"]

    ps.enforce_style()
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2),
                             constrained_layout=True)
    ang = np.linspace(0, 2 * np.pi, 96, endpoint=False)
    R = 0.12
    kk = np.array([[R * np.cos(a), R * np.sin(a), 0.0] for a in ang])
    Sb = np.array([floquet_spin_observables(model, d1, q, ops,
                                            n_step=300, m_sp=12)["S_bar"][0]
                   for q in kk])
    Ss = np.array([spin_texture(model.Hk(q), ops)["spins"][0] for q in kk])
    ax = axes[0]
    ax.quiver(kk[:, 0], kk[:, 1], Ss[:, 0], Ss[:, 1], color="black",
              pivot="middle", width=0.006)
    ax.set_xlabel("$k_x$ (1/Angstrom)")
    ax.set_ylabel("$k_y$ (1/Angstrom)")
    ax.set_aspect("equal")
    ax = axes[1]
    ax.quiver(kk[:, 0], kk[:, 1], Sb[:, 0], Sb[:, 1], color="#c0392b",
              pivot="middle", width=0.006)
    ax.set_xlabel("$k_x$ (1/Angstrom)")
    ax.set_aspect("equal")
    ax = axes[2]
    mm = np.linalg.norm(Sb - Ss, axis=1)
    ax.plot(ang * 180 / np.pi, mm * 1000, color="black", lw=1.4)
    ax.set_xlabel("k angle (deg)")
    ax.set_ylabel("|S_bar - S_static| (10$^{-3}$)")
    out = os.path.join(outdir, "benchmark_floquet_rashba.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out
