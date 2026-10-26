# -*- coding: utf-8 -*-
"""Benchmark 5: Rashba ferromagnet (exchange vs SOC).

H(k) = c k^2 I + alpha (k_y sigma_x - k_x sigma_y) + M sigma_z

Checks:
  * at k=0 the splitting is exactly 2M (exchange), independent of alpha;
  * at finite k the splitting depends on alpha and M (SOC + exchange);
  * the k=0 spin is +/- z with |<S>| = 1/2;
  * gamma=0 Floquet modes recover the static spectrum.
"""
from __future__ import annotations

import numpy as np

from ..core.spin_operator import SPIN_HALF, SIGMA
from ..core.validation import ValidationReport, PASS, FAIL
from ..drives.peierls import Drive
from ..observables.spin_texture import spin_texture
from ..solvers.evolution import floquet_modes


class MagneticRashbaModel:
    kind = "analytic"
    name = "RashbaFM"
    nw = 2
    a_ref = 1.0

    def __init__(self, c=0.5, alpha=0.1, M=0.2):
        self.c, self.alpha, self.M = float(c), float(alpha), float(M)

    def Hk(self, k):
        k = np.asarray(k, float).reshape(3)
        kx, ky = k[0], k[1]
        H = self.c * (kx * kx + ky * ky) * np.eye(2, dtype=complex)
        H += self.alpha * (ky * SIGMA[0] - kx * SIGMA[1])
        H += self.M * SIGMA[2]
        return H

    def spin_ops(self):
        return SPIN_HALF


def run(outdir=None):
    rep = ValidationReport("Benchmark 5: Rashba ferromagnet (exchange vs SOC)")
    model = MagneticRashbaModel(c=0.5, alpha=0.1, M=0.2)
    k0 = np.zeros(3)
    E0 = np.linalg.eigvalsh(model.Hk(k0))
    rep.add("exchange_splitting_at_Gamma",
            PASS if abs((E0[1] - E0[0]) - 2 * model.M) < 1e-12 else FAIL,
            f"Delta E(Gamma) = {E0[1]-E0[0]:.6f} eV, 2M = "
            f"{2*model.M:.6f} eV (tol 1e-12)")
    s0 = spin_texture(model.Hk(k0), SPIN_HALF)["spins"]
    rep.add("Gamma_spin_is_out_of_plane",
            PASS if (abs(s0[0, 2]) > 0.499 and abs(s0[1, 2]) > 0.499
                     and np.max(np.abs(s0[:, :2])) < 1e-12) else FAIL,
            f"|<S_z>| = {abs(s0[:,2]).min():.6f}, "
            f"max|S_xy| = {np.max(np.abs(s0[:,:2])):.2e}")

    k1 = np.array([0.10, 0.0, 0.0])
    E1 = np.linalg.eigvalsh(model.Hk(k1))
    D1 = float(E1[1] - E1[0])
    rep.add("finite_k_splitting_has_soc",
            PASS if abs(D1 - 2 * model.M) > 1e-6 else FAIL,
            f"Delta E(k=0.1) = {D1:.6f} eV differs from 2M by "
            f"{D1-2*model.M:.2e} eV (SOC contribution)")

    drive = Drive(hw=0.4, gamma=0.0, a_ref=1.0, polarization="circular")
    fl = floquet_modes(model, drive, k1, n_step=120, want_wspec=True, m_sp=12)
    eps_sorted = np.sort(fl["eps_u"])
    d = (eps_sorted - E1 + drive.hw / 2) % drive.hw - drive.hw / 2
    err = float(np.max(np.abs(d)))
    rep.add("gamma0_static_limit", PASS if err < 1e-10 else FAIL,
            f"max folded |eps_F - E_static| at k=0.1 = {err:.2e} eV "
            f"(tol 1e-10, modulo hw={drive.hw:g})")
    return dict(report=rep, model=model, floquet=fl)


def _unused():
    return np
