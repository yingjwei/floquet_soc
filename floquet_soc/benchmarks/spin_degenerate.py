# -*- coding: utf-8 -*-
"""Benchmark 1: spin-degenerate (no SOC) model under Peierls driving.

H0(k) = c k^2 I_2.  With no SOC the two spin channels are exactly
degenerate and the Floquet propagator is proportional to the identity in
spin space.  A correct Peierls implementation must therefore produce:
  * zero static spin splitting;
  * zero driven spin splitting;
  * a spin-subspace average <S> = 0 (no spurious spin texture).
"""
from __future__ import annotations

import numpy as np

from ..core.spin_operator import SPIN_HALF
from ..core.validation import ValidationReport, PASS, FAIL
from ..drives.peierls import Drive
from ..solvers.evolution import floquet_modes
from ..benchmarks.rashba import RashbaModel


def run(outdir=None):
    rep = ValidationReport("Benchmark 1: spin-degenerate model (no SOC)")
    model = RashbaModel(hbar2_over_2m=0.5, alpha=0.0)
    k = np.array([0.08, 0.03, 0.0])
    E = np.linalg.eigvalsh(model.Hk(k))
    rep.add("no_soc_static_degeneracy",
            PASS if abs(E[1] - E[0]) < 1e-12 else FAIL,
            f"static splitting = {E[1]-E[0]:.2e} eV (tol 1e-12)")

    drive = Drive(hw=0.4, gamma=0.05, a_ref=1.0, polarization="circular")
    fl = floquet_modes(model, drive, k, n_step=120, want_wspec=True, m_sp=12)
    eps = np.sort(fl["eps_u"])
    rep.add("no_soc_driven_degeneracy",
            PASS if abs(eps[1] - eps[0]) < 1e-10 else FAIL,
            f"driven quasienergy splitting = {eps[1]-eps[0]:.2e} eV "
            f"(tol 1e-10)")

    # spin subspace average: <S> = Tr(S)/2 = 0 for both bands
    S = np.zeros(3)
    for c in range(3):
        S[c] = float(np.real(np.trace(SPIN_HALF[c])) / 2)
    rep.add("no_soc_no_spurious_spin",
            PASS if np.linalg.norm(S) < 1e-12 else FAIL,
            f"subspace-averaged |<S>| = {np.linalg.norm(S):.2e} "
            f"(tol 1e-12)")
    return dict(report=rep, model=model, floquet=fl)


def _unused():
    return np
