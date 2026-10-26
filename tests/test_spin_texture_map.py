# -*- coding: utf-8 -*-
"""Tests for the generic spin-texture map program and DFT cross-check."""
import os
import shutil
import tempfile
import unittest

import numpy as np

import sys
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.observables.spin_texture_map import (        # noqa: E402
    resolve_band, full_kpts, window_kpts, compute_map)
from floquet_soc.validation.dft_spin import (                 # noqa: E402
    read_spin_texture_dat, compare, save_and_plot)
from floquet_soc.observables.floquet_points import resolve_points  # noqa: E402


class FakeModel:
    """Two-band model with a known spin texture and a lattice."""

    nw = 2

    def __init__(self, gap=2.0):
        self.gap = gap
        self.lattice = np.diag([3.0, 3.0, 20.0])

    def Hk(self, k):
        return np.diag([-0.5 * self.gap, 0.5 * self.gap])


def fake_spin_at(k):
    # Pauli matrices (spin expectation units)
    return np.array([[[0, 1], [1, 0]],
                     [[0, -1j], [1j, 0]],
                     [[1, 0], [0, -1]]], dtype=complex)


class TestSpinTextureMap(unittest.TestCase):
    def test_resolve_band(self):
        m = FakeModel()
        cfg = dict(num_electrons=1)
        self.assertEqual(resolve_band(cfg, m, None), 0)
        self.assertEqual(resolve_band(cfg, m, "vbm"), 0)
        self.assertEqual(resolve_band(cfg, m, "cbm"), 1)
        self.assertEqual(resolve_band(cfg, m, 0), 0)
        with self.assertRaises(ValueError):
            resolve_band(cfg, m, 7)

    def test_grids(self):
        kf = full_kpts({}, n=7, kmax=1.0)
        self.assertEqual(kf.shape, (49, 3))
        lat = np.diag([3.0, 3.0, 20.0])
        kw, xy = window_kpts(lat, np.zeros(3), 0.1, n=5)
        self.assertEqual(kw.shape, (25, 3))
        self.assertEqual(xy.shape, (25, 2))
        self.assertAlmostEqual(float(xy[:, 0].min()), -0.1)
        self.assertAlmostEqual(float(xy[:, 1].max()), 0.1)

    def test_compute_map_and_mask(self):
        m = FakeModel(gap=2.0)
        kf = full_kpts({}, n=3, kmax=0.5)
        tex = compute_map(m, fake_spin_at, 1, kf)
        self.assertEqual(tex["energies"].shape, (9, 2))
        self.assertEqual(tex["spins_band"].shape, (9, 3))
        self.assertTrue(np.allclose(np.abs(tex["spins_band"][:, 2]), 1.0))
        self.assertTrue(np.allclose(tex["spins_band"][:, :2], 0.0))
        # degenerate model: gap=0, in-plane arrows suppressed by mask
        m0 = FakeModel(gap=0.0)
        tex0 = compute_map(m0, fake_spin_at, 1, kf, mask_split=1e-3)
        self.assertTrue(np.all(tex0["split"] < 1e-3))


class TestFloquetPoints(unittest.TestCase):
    def test_resolve_labels_and_triples(self):
        cfg = dict(kpath=[["K", [1 / 3, 1 / 3, 0.0]],
                          ["M", [0.5, 0.0, 0.0]]])
        pts = resolve_points(cfg, "K,Kp")
        labs = [p[0] for p in pts]
        self.assertEqual(labs[0], "K")
        self.assertAlmostEqual(pts[1][1][0], 2 / 3, places=6)
        self.assertAlmostEqual(pts[1][1][1], 2 / 3, places=6)
        tri = resolve_points(cfg, "0.5,0,0;0,0,0")
        self.assertEqual(len(tri), 2)
        self.assertTrue(np.allclose(tri[0][1], [0.5, 0, 0]))


class TestDftSpin(unittest.TestCase):
    def test_read_compare_save(self):
        txt = """#kx ky Energy Sx Sy Sz
#Grid-Size: 3 3
#Band-Index: 2
0.0 0.0 1.0 0.0 0.0 -1.0
0.1 0.0 1.0 0.0 0.0 -1.0
0.0 0.1 1.0 0.0 0.0 -1.0
"""
        td = tempfile.mkdtemp(prefix="floquet_soc_dftspin_")
        self.addCleanup(shutil.rmtree, td, True)
        p = os.path.join(td, "dft.dat")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(txt)
        d = read_spin_texture_dat(p)
        self.assertEqual(d["band_index"], 2)
        self.assertEqual(d["grid"], [3, 3])
        m = FakeModel()
        res = compare(m, fake_spin_at, d, band=1, mask_split=0.01)
        self.assertEqual(res["band_wannier"], 1)
        self.assertLess(res["energy"]["rms_resid_ev"], 1e-8)
        self.assertGreater(res["scale"]["global_fit"], 0.9)
        # numbers only: the figures are a user-facing artifact, not a test
        # assertion, and must never be left inside the source tree
        out = save_and_plot(res, m, fake_spin_at, td, prefix="cmp",
                            plots=False)
        self.assertTrue(os.path.exists(out["json"]))
        self.assertTrue(os.path.exists(out["npz"]))
        self.assertNotIn("texture", out)


if __name__ == "__main__":
    unittest.main(verbosity=2)
