# -*- coding: utf-8 -*-
"""Stage C tests: SOC splitting, effective-term fitting, PST metrics."""
import json
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.cli import main as cli_main
from floquet_soc.core.spin_operator import SPIN_HALF, SIGMA
from floquet_soc.observables.soc_splitting import (
    adjacent_splitting, spin_paired_splitting, fit_linear_model,
    soc_off_difference)
from floquet_soc.observables.pst import (pst_collinearity, analyse_texture,
                                        axis_angle)
from floquet_soc.observables.spin_texture import spin_texture
from floquet_soc.benchmarks.rashba import RashbaModel
from floquet_soc.benchmarks.rashba_dresselhaus import RashbaDresselhausModel

TMP = os.path.join(HERE, "outputs", "_tests_c")


class TestSocSplitting(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = RashbaModel(hbar2_over_2m=0.5, alpha=0.1)
        ax = np.linspace(-0.20, 0.20, 13)
        cls.kpts = np.array([[x, y, 0.0] for x in ax for y in ax])
        cls.kmag = np.linalg.norm(cls.kpts[:, :2], axis=1)
        cls.E = cls.model.static_bands(cls.kpts)
        cls.S = np.array([spin_texture(cls.model.Hk(k), SPIN_HALF)["spins"]
                          for k in cls.kpts])

    def test_adjacent_splitting_analytic(self):
        d = adjacent_splitting(self.E)[:, 0]
        dev = float(np.max(np.abs(d - 2 * self.model.alpha * self.kmag)))
        self.assertLess(dev, 1e-12)

    def test_spin_paired_splitting(self):
        pairs = spin_paired_splitting(self.E, self.S)
        self.assertGreater(len(pairs), 100)
        gaps = np.array([p["gap"] for p in pairs])
        kmag = np.array([self.kmag[p["k"]] for p in pairs])
        dev = float(np.max(np.abs(gaps - 2 * self.model.alpha * kmag)))
        self.assertLess(dev, 1e-10)
        anti = np.array([p["spin_antiparallel"] for p in pairs])
        self.assertGreater(float(np.min(anti)), 0.99)

    def test_fit_linear_model_recovers_rashba(self):
        def basis(k):
            return [k[0] ** 2 + k[1] ** 2, k[1], -k[0]]
        ops = [np.eye(2, dtype=complex), SIGMA[0], SIGMA[1]]
        res = fit_linear_model(self.model.Hk, self.kpts,
                               [lambda k, i=i: basis(k)[i] for i in range(3)],
                               ops)
        c = res["coefficients"]
        self.assertAlmostEqual(c[0], 0.5, places=10)
        self.assertAlmostEqual(c[1], 0.1, places=10)
        self.assertAlmostEqual(c[2], 0.1, places=10)
        self.assertLess(res["rms"], 1e-12)

    def test_soc_off_requires_matched_basis(self):
        self.assertEqual(soc_off_difference(None, {})["status"],
                         "NEEDS_INPUT")
        a = {"files": {"win": "a", "hr": "b", "centres": "c"}}
        self.assertEqual(soc_off_difference(a, {})["status"], "NEEDS_INPUT")


class TestPST(unittest.TestCase):
    def test_rashba_pst_eigenvalues(self):
        m = RashbaModel(alpha=0.1)
        circ = np.array([[0.2 * np.cos(a), 0.2 * np.sin(a), 0.0]
                         for a in np.linspace(0, 2 * np.pi, 49, endpoint=False)])
        S = np.array([spin_texture(m.Hk(k), SPIN_HALF)["spins"][0]
                      for k in circ])
        pst = pst_collinearity(S)
        self.assertLess(float(np.max(np.abs(
            pst["eigenvalues"] - np.array([0.5, 0.5, 0.0])))), 1e-8)

    def test_rd_pst_collinear_axis(self):
        m = RashbaDresselhausModel(alpha=0.1, beta=0.1)
        circ = np.array([[0.2 * np.cos(a), 0.2 * np.sin(a), 0.0]
                         for a in np.linspace(0, 2 * np.pi, 49, endpoint=False)])
        E = np.array([np.linalg.eigvalsh(m.Hk(k)) for k in circ])
        S = np.array([spin_texture(m.Hk(k), SPIN_HALF)["spins"] for k in circ])
        res = analyse_texture(E, S, bands=[0])
        self.assertGreater(res["n_states"], 40)
        self.assertAlmostEqual(res["collinearity"], 1.0, places=8)
        axis = np.array([1.0, -1.0, 0.0]) / np.sqrt(2)
        self.assertLess(axis_angle(res["pst_axis"], axis), 1e-8)

    def test_rd_noncollinear(self):
        m = RashbaDresselhausModel(alpha=0.1, beta=0.03)
        circ = np.array([[0.2 * np.cos(a), 0.2 * np.sin(a), 0.0]
                         for a in np.linspace(0, 2 * np.pi, 49, endpoint=False)])
        E = np.array([np.linalg.eigvalsh(m.Hk(k)) for k in circ])
        S = np.array([spin_texture(m.Hk(k), SPIN_HALF)["spins"] for k in circ])
        res = analyse_texture(E, S, bands=[0])
        self.assertLess(res["collinearity"], 0.99)


class TestStageCCLI(unittest.TestCase):
    def test_soc_analytic(self):
        rc = cli_main(["soc", "--config",
                       os.path.join(HERE, "configs", "rashba.json"),
                       "--nk", "2", "--outdir", os.path.join(TMP, "soc")])
        self.assertEqual(rc, 0)
        self.assertTrue(os.path.exists(os.path.join(TMP, "soc",
                                                    "soc_splitting.npz")))

    def test_pst_analytic(self):
        rc = cli_main(["pst", "--config",
                       os.path.join(HERE, "configs", "rashba_dresselhaus.json"),
                       "--outdir", os.path.join(TMP, "pst")])
        self.assertEqual(rc, 0)
        self.assertTrue(os.path.exists(os.path.join(TMP, "pst", "pst.npz")))
        d = np.load(os.path.join(TMP, "pst", "pst.npz"))
        self.assertGreater(d["collinearity"], 0.99)

    def test_soc_pbnv(self):
        from floquet_soc.cli import load_config
        cfg = load_config(os.path.join(HERE, "configs", "pbnv.json"))
        for key in ("chk", "eig"):
            if not os.path.exists(cfg["files"][key]):
                self.skipTest("PbNV checkpoint data not available")
        rc = cli_main(["soc", "--config",
                       os.path.join(HERE, "configs", "pbnv.json"),
                       "--nk", "5", "--outdir", os.path.join(TMP, "pbnv")])
        self.assertEqual(rc, 0)


class TestLinearRashbaFit(unittest.TestCase):
    """Unit tests for the linear Rashba fit ported from the user's script."""

    @staticmethod
    def _synthetic_pair(alpha=0.1, a=0.5, e0=1.0, kmax=0.12, nk=121):
        k = np.linspace(-kmax, kmax, nk)
        lower = a * k * k + e0 - alpha * np.abs(k)
        upper = a * k * k + e0 + alpha * np.abs(k)
        return k, np.stack([lower, upper])

    def test_recovers_linear_rashba_alpha(self):
        from floquet_soc.observables.rashba_fit import fit_gap_edge_pair
        k, e = self._synthetic_pair(alpha=0.1)
        r = fit_gap_edge_pair(k, e)
        self.assertEqual(r["status"], "LINEAR")
        self.assertAlmostEqual(r["alpha_r"], 0.1, delta=1e-8)
        self.assertAlmostEqual(r["alpha_from_split"], 0.1, delta=1e-8)
        self.assertGreater(r["r2_split"], 0.999)

    def test_rejects_degenerate_parabola(self):
        from floquet_soc.observables.rashba_fit import fit_gap_edge_pair
        k = np.linspace(-0.12, 0.12, 121)
        e = 0.5 * k * k + 1.0
        r = fit_gap_edge_pair(k, np.stack([e, e]))
        # no genuine +-alpha splitting: never LINEAR; the original script's
        # fallback window may still be reported as FALLBACK
        self.assertNotEqual(r["status"], "LINEAR")

    def test_quadratic_splitting_is_not_linear(self):
        from floquet_soc.observables.rashba_fit import fit_gap_edge_pair
        k = np.linspace(-0.12, 0.12, 121)
        lower = 1.0 + 0.5 * k * k - 5.0 * k * k
        upper = 1.0 + 0.5 * k * k + 5.0 * k * k
        r = fit_gap_edge_pair(k, np.stack([lower, upper]))
        self.assertNotEqual(r["status"], "LINEAR")


if __name__ == "__main__":
    unittest.main(verbosity=2)
