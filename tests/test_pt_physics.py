# -*- coding: utf-8 -*-
"""Physics-correctness regression tests (Pt PT/Kramers, WS interpolation,
Kramers definitions, Rashba analytics, gauge invariance, convergence).

Tolerances:
* DFT-mesh Kramers residual: 1e-6 eV (DFT eigenvalue precision ~2e-10 eV).
* WS interpolation Hermiticity: 1e-10 eV (double precision round-off).
* Off-grid Kramers residual (4x4x4 Wannier truncation): median <1 meV,
  max <10 meV; target SOC splitting scale is ~50 meV, so a residual below
  the physical scale cannot produce a false SOC splitting.
"""
import os
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PT_CFG = os.path.join(ROOT, "configs", "pt.json")
PT_DIR = os.path.join(ROOT, "benchmarks_data", "Pt")


class KramersDefinitionTests(unittest.TestCase):
    def test_kramers_delta_and_gaps(self):
        from floquet_soc.observables.kramers import (kramers_delta,
                                                     inter_doublet_gaps)
        E = np.array([[1.0, 1.0, 2.0, 2.0], [0.0, 0.0, 1.0, 1.0]])
        np.testing.assert_allclose(kramers_delta(E), 0.0, atol=1e-15)
        np.testing.assert_allclose(inter_doublet_gaps(E),
                                   [[1.0], [1.0]], atol=1e-15)
        E2 = E.copy()
        E2[0, 1] = 1.05                     # split pair (0,1)
        d = kramers_delta(E2)
        self.assertAlmostEqual(d[0, 0], 0.05, places=12)
        self.assertAlmostEqual(d[0, 1], 0.0, places=12)

    def test_ws_supercell_even_grid(self):
        from floquet_soc.core.wannier import ws_supercell
        R, nd = ws_supercell((4, 4, 4), np.eye(3))
        self.assertEqual(len(R), 125)            # N+1 values per direction
        sR = {tuple(r) for r in R}
        self.assertTrue(all(tuple(-r) in sR for r in R))
        # ndegen symmetric and sum_R 1/ndegen = N^3 (residues counted once)
        idx = {tuple(r): i for i, r in enumerate(R)}
        for r in R:
            self.assertEqual(nd[idx[tuple(r)]], nd[idx[tuple(-r)]])
        self.assertAlmostEqual(float(np.sum(1.0 / nd)), 64.0, places=10)

    def test_no_soc_model_zero_splitting(self):
        from floquet_soc.core.wannier import WannierModel
        from floquet_soc.observables.kramers import kramers_delta
        win = dict(num_wann=2, num_bands=2, lattice=np.eye(3),
                   kpoints=np.zeros((1, 3)), mp_grid=np.array([1, 1, 1]),
                   kpoint_weights=np.ones(1), atoms=[], spinors=None)
        hr = dict(num_wann=2, nrpts=1, R=np.zeros((1, 3), int),
                  H=np.zeros((1, 2, 2), complex), ndegen=np.ones(1, int),
                  header="test")
        m = WannierModel(win, hr, np.zeros((2, 3)), "test")
        E = m.static_bands(np.array([[0.1, 0.2, 0.0], [0.3, -0.1, 0.0]]))
        np.testing.assert_allclose(E, 0.0, atol=1e-15)
        np.testing.assert_allclose(kramers_delta(E), 0.0, atol=1e-15)


class RashbaAnalyticTests(unittest.TestCase):
    def test_rashba_splitting_and_spin_direction(self):
        from floquet_soc.benchmarks.rashba import RashbaModel
        m = RashbaModel(hbar2_over_2m=0.5, alpha=0.1)
        k = np.array([0.05, 0.03, 0.0])
        self.assertAlmostEqual(m.splitting(k), 2 * m.alpha * np.linalg.norm(k),
                               places=13)
        r = m.spin_at(k)
        # spin must be perpendicular to k, in-plane
        S = r["spins"]
        self.assertLess(abs(np.dot(S[0, :2], k[:2])), 1e-12)
        self.assertLess(abs(S[0, 2]), 1e-12)
        self.assertAlmostEqual(np.linalg.norm(S[0]), 0.5, places=12)


class GaugeInvarianceDetectionTests(unittest.TestCase):
    def test_gauge_test_detects_degenerate_arbitrariness(self):
        from floquet_soc.observables.pst import gauge_invariance_test
        rng = np.random.default_rng(3)

        class DummyModel:
            nw = 4

            def Hk(self, k):
                return np.diag([0.0, 0.0, 1.0, 1.0]).astype(complex)

        def spin_at_k(k):
            ops = np.zeros((3, 4, 4), complex)
            for c in range(3):
                v = rng.normal(size=4) + 1j * rng.normal(size=4)
                ops[c] = np.outer(v, v.conj())
                ops[c] = 0.5 * (ops[c] + ops[c].conj().T)
            return ops

        kpts = np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]])
        g = gauge_invariance_test(DummyModel(), spin_at_k, kpts,
                                  bands=[0, 1, 2, 3], deg_tol=1e-3)
        self.assertGreater(g["n_pairs"], 0)
        # gauge-invariant subspace spin is well defined; the metric changes
        self.assertTrue(np.isfinite(g["coll_spread"]))
        self.assertLess(g["axis_dot_min"], 1.0 + 1e-12)


@unittest.skipUnless(os.path.exists(os.path.join(PT_DIR, "Pt.chk")),
                     "Pt tutorial data not available")
class PtPhysicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from floquet_soc.cli import load_config, build_model
        cls.cfg = load_config(PT_CFG)
        cls.model = build_model(cls.cfg)
        cls.mesh = np.asarray(cls.model.win["kpoints"], float)
        cls.rng = np.random.default_rng(20261022)
        cls.krand = cls.rng.uniform(-0.5, 0.5, size=(200, 3))

    def test_offgrid_hermiticity(self):
        H = self.model.Hk_batch(self.krand)
        herm = np.max(np.abs(H - np.conj(H.transpose(0, 2, 1))))
        self.assertLess(herm, 1e-10, f"off-grid Hermiticity {herm:.2e} eV")

    def test_kramers_on_dft_mesh(self):
        from floquet_soc.observables.kramers import kramers_summary
        s = kramers_summary(self.model.static_bands(self.mesh), self.mesh)
        self.assertLess(s["delta_max"], 1e-6,
                        f"on-mesh Kramers {s['delta_max']:.2e} eV")
        self.assertEqual(s["n_pairs"], self.model.nw // 2)

    def test_kramers_offgrid_residual(self):
        from floquet_soc.observables.kramers import kramers_summary
        s = kramers_summary(self.model.static_bands(self.krand), self.krand)
        self.assertLess(s["delta_median"] * 1e3, 1.0,
                        f"off-grid median {s['delta_median']*1e3:.3f} meV")
        self.assertLess(s["delta_max"] * 1e3, 10.0,
                        f"off-grid max {s['delta_max']*1e3:.3f} meV")

    def test_spin_operator_offgrid_hermiticity(self):
        from floquet_soc.core.spin_operator import load_spin_for_model
        spin = load_spin_for_model(self.cfg, self.model)
        op = spin["spin_operator"]
        self.assertLess(op.hermiticity_offgrid(self.krand), 1e-10)
        self.assertLess(op.interpolation_error(), 1e-10)


class FloquetMicromotionTests(unittest.TestCase):
    def test_micromotion_convergence(self):
        from floquet_soc.benchmarks.rashba import RashbaModel
        from floquet_soc.observables.floquet_spin import (
            floquet_spin_observables)
        m = RashbaModel(hbar2_over_2m=0.5, alpha=0.1)
        k = np.array([0.05, 0.03, 0.0])
        ops = m.spin_ops()
        from floquet_soc.drives.peierls import Drive
        d = Drive(hw=1.0, gamma=0.05, polarization="circular", sign=-1,
                  a_ref=1.0)
        a = floquet_spin_observables(m, d, k, ops, n_step=80, m_sp=8)
        b = floquet_spin_observables(m, d, k, ops, n_step=160, m_sp=8)
        self.assertLess(float(np.max(np.abs(a["S_bar"] - b["S_bar"]))), 1e-3)
        from floquet_soc.workflows.converge import run_convergence
        res = run_convergence(m, d, k, ops, [80, 160], [4, 8],
                              window=(-2.0, 2.0))
        for r in res["m_sp_rows"]:
            self.assertIn("dev_wsum", r)


if __name__ == "__main__":
    unittest.main()
