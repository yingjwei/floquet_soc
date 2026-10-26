# -*- coding: utf-8 -*-
"""A4: quasienergy fold unification + band continuity along a k path.

Before this batch ``compute_path`` re-sorted the Floquet modes by quasienergy at
every k point independently, so a band crossing (or a state drifting across the
first-Floquet-zone edge) swapped column indices and produced non-physical +/- hw
jumps in ``eps``/``Ebar``.  The fix tracks each state across k by eigenvector
overlap and then unwraps the folded quasienergy.
"""
from __future__ import annotations

import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.solvers import fold_quasienergy                     # noqa: E402
from floquet_soc.observables.floquet_bands import (                  # noqa: E402
    _match_columns, _unwrap_along_path, compute_path)
from floquet_soc.benchmarks.rashba import RashbaModel                # noqa: E402
from floquet_soc.drives.peierls import Drive                         # noqa: E402


class TestFoldQuasienergy(unittest.TestCase):

    def test_fold_is_in_first_zone(self):
        hw = 0.4
        vals = fold_quasienergy(np.array([0.0, 0.2, -0.2, 0.4, -0.4]), hw)
        self.assertTrue(np.all((vals >= -hw / 2) & (vals < hw / 2)))

    def test_fold_is_modulo_hw(self):
        hw = 0.4
        x = 0.13
        self.assertAlmostEqual(fold_quasienergy(x, hw),
                               fold_quasienergy(x + hw, hw))
        self.assertAlmostEqual(fold_quasienergy(x, hw),
                               fold_quasienergy(x - hw, hw))


class TestMatchColumns(unittest.TestCase):

    def test_swapped_columns(self):
        prev = np.eye(2)
        cur = np.array([[0.0, 1.0], [1.0, 0.0]])   # columns swapped
        self.assertEqual(list(_match_columns(prev, cur)), [1, 0])

    def test_identity(self):
        prev = np.eye(3)
        self.assertEqual(list(_match_columns(prev, prev)), [0, 1, 2])


class TestUnwrapAlongPath(unittest.TestCase):

    def test_zone_edge_jump_is_removed(self):
        hw = 1.0
        tracked = np.array([[0.4], [-0.4]])         # -0.4 is 0.6 - hw
        out = _unwrap_along_path(tracked, hw)
        self.assertTrue(np.allclose(out[:, 0], [0.4, 0.6]),
                        f"got {out[:, 0]}")

    def test_columns_are_independent(self):
        hw = 1.0
        tracked = np.array([[0.4, -0.4], [-0.4, -0.2]])
        out = _unwrap_along_path(tracked, hw)
        self.assertTrue(np.allclose(out[:, 0], [0.4, 0.6]))
        self.assertTrue(np.allclose(out[:, 1], [-0.4, -0.2]))


class TestComputePathContinuity(unittest.TestCase):

    def test_path_has_continuous_unwrapped_eps(self):
        model = RashbaModel(hbar2_over_2m=0.5, alpha=0.15)
        drive = Drive(hw=1.0, gamma=0.2, a_ref=1.0)
        kpts = np.array([[x, 0.0, 0.0] for x in np.linspace(0.0, 0.25, 9)])
        res = compute_path(model, drive, kpts, n_step=60, m_sp=4, sigma=0.02)
        hw = res["hw"]
        self.assertTrue(np.all((res["eps"] >= -hw / 2) & (res["eps"] < hw / 2)),
                        "eps must stay folded")
        self.assertTrue(np.all(np.isfinite(res["eps_unwrapped"])))
        # no non-physical +/- hw jump survives along the tracked path
        jumps = np.abs(np.diff(res["eps_unwrapped"], axis=0)).max()
        self.assertLess(jumps, hw / 2 + 1e-6,
                        f"unwrapped eps still has a {jumps:.3f} jump")


if __name__ == "__main__":
    unittest.main()
