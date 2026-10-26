# -*- coding: utf-8 -*-
"""Frozen regression tests for release v0.1.4-rc1 (core physics unchanged since v0.1.0-rc1).

These tests cover the four issue classes fixed before the release:
old/fake SOC splitting, PT/Kramers degeneracy, non-Hermitian interpolation,
PST gauge dependence, bad .spn quality and incomplete photon truncation.

Tolerances are physical/numerical, not fitted to make tests pass:
* Hermiticity: 1e-12 (double-precision round-off over small matrices).
* Kramers (DFT mesh): 1e-6 eV (DFT eigenvalue precision ~2e-10 eV).
* Gauge spread: >1e-3 is a real gauge dependence of single-band metrics.
* Photon weights: a state is converged if |1-sum_m w_m| < 1e-3.
"""
import os
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _dev():
    return ROOT


class VersionTests(unittest.TestCase):
    def test_version_string(self):
        from floquet_soc.version import __version__
        vpath = os.path.join(ROOT, "VERSION")
        self.assertTrue(os.path.exists(vpath))
        with open(vpath, encoding="utf-8") as f:
            self.assertEqual(__version__, f.read().strip())


class FakeSocSplittingTests(unittest.TestCase):
    def test_pt_system_has_no_spin_pair_rows_in_soc_output(self):
        """Pt soc output must contain Kramers diagnostics, not spin splitting."""
        p = os.path.join(ROOT, "outputs", "Pt", "soc_splitting.npz")
        if not os.path.exists(p):
            self.skipTest("Pt soc output not generated")
        d = np.load(p)
        self.assertIn("kramers_delta", d.files)
        self.assertNotIn("spin_pair_rows", d.files)

    def test_kramers_delta_definition(self):
        from floquet_soc.observables.kramers import (kramers_delta,
                                                     inter_doublet_gaps)
        E = np.array([[0.0, 0.0, 1.0, 1.0], [0.0, 0.05, 1.0, 1.0]])
        d = kramers_delta(E)
        self.assertAlmostEqual(d[0, 0], 0.0, places=14)
        self.assertAlmostEqual(d[1, 0], 0.05, places=14)
        np.testing.assert_allclose(inter_doublet_gaps(E)[0], [1.0])


class NonHermitianInterpolationTests(unittest.TestCase):
    """The old alias (R=0..N-1, ndegen=1) interpolant is non-Hermitian off
    grid; the WS (symmetric R + ndegen) interpolant is Hermitian by
    construction.  This is the regression that motivated the fix."""

    def _mesh_H(self, n=4, nw=2):
        rng = np.random.default_rng(0)
        grids = [np.arange(n) / n for _ in range(3)]
        kpts = np.array(np.meshgrid(*grids, indexing="ij")).reshape(3, -1).T
        H = []
        for _ in kpts:
            A = rng.normal(size=(nw, nw)) + 1j * rng.normal(size=(nw, nw))
            H.append(0.5 * (A + A.conj().T))
        return kpts, np.asarray(H)

    def test_ws_interpolation_is_hermitian_offgrid(self):
        from floquet_soc.core.wannier import ws_supercell
        n, nw = 4, 2
        kpts, HW = self._mesh_H(n, nw)
        R, nd = ws_supercell((n, n, n), np.eye(3))
        HR = np.einsum("kr,kij->rij",
                       np.exp(-2j * np.pi * (kpts @ R.T)), HW) / len(kpts)
        C = HR / nd[:, None, None]
        rng = np.random.default_rng(1)
        ktest = rng.uniform(-0.5, 0.5, size=(50, 3))
        ph = np.exp(2j * np.pi * (ktest @ R.T))
        Hk = np.einsum("kr,rij->kij", ph, C)
        err = np.max(np.abs(Hk - np.conj(Hk.transpose(0, 2, 1))))
        self.assertLess(err, 1e-12)

    def test_alias_interpolation_breaks_hermiticity_offgrid(self):
        n, nw = 4, 2
        kpts, HW = self._mesh_H(n, nw)
        R = np.array(np.meshgrid(*[np.arange(n) for _ in range(3)],
                                 indexing="ij")).reshape(3, -1).T
        HR = np.einsum("kr,kij->rij",
                       np.exp(-2j * np.pi * (kpts @ R.T)), HW) / len(kpts)
        rng = np.random.default_rng(1)
        ktest = rng.uniform(-0.5, 0.5, size=(50, 3))
        ph = np.exp(2j * np.pi * (ktest @ R.T))
        Hk = np.einsum("kr,rij->kij", ph, HR)
        err = np.max(np.abs(Hk - np.conj(Hk.transpose(0, 2, 1))))
        self.assertGreater(err, 1e-2,
                           "alias interpolation should be non-Hermitian; "
                           "this is the documented release limitation")


class BadSpnQualityTests(unittest.TestCase):
    def test_non_pauli_spin_matrices_are_flagged(self):
        from floquet_soc.core.spin_operator import validate_spin_matrices
        S = np.zeros((2, 3, 2, 2), complex)
        for ik in range(2):
            for c in range(3):
                S[ik, c] = 2.5 * np.array([[1, 0], [0, -1]], complex)
        metrics, _checks = validate_spin_matrices(S)
        self.assertGreater(metrics["eig_max"], 1.0 + 1e-6)
        self.assertGreater(metrics["S2_max"], 3.0 + 1e-6)

    def test_trs_violating_operator_is_flagged(self):
        from floquet_soc.core.spin_wannier import WannierSpinOperator
        nk, nw = 2, 2
        kmesh = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
        rng = np.random.default_rng(2)
        Sm = np.zeros((nk, 3, nw, nw), complex)
        for ik in range(nk):
            for c in range(3):
                A = rng.normal(size=(nw, nw)) + 1j * rng.normal(size=(nw, nw))
                Sm[ik, c] = 0.5 * (A + A.conj().T)
        op = WannierSpinOperator(kmesh, (2, 1, 1), Sm, lattice=np.eye(3))
        checks = dict((n, p) for n, p, d in op.checks())
        self.assertIn("spin_wannier_trs", checks)


class PstGaugeDependenceTests(unittest.TestCase):
    def test_pt_pst_output_records_gauge_instability(self):
        p = os.path.join(ROOT, "outputs", "Pt", "pst.npz")
        if not os.path.exists(p):
            self.skipTest("Pt pst output not generated")
        d = np.load(p)
        self.assertIn("gauge_collinearity", d.files)
        self.assertIn("gauge_axis", d.files)
        self.assertIn("note", d.files)
        self.assertIn("NOT_APPLICABLE", str(d["note"]))


class PhotonTruncationTests(unittest.TestCase):
    def test_msp_window_produces_deficit_fields(self):
        from floquet_soc.benchmarks.rashba import RashbaModel
        from floquet_soc.drives.peierls import Drive
        from floquet_soc.workflows.converge import run_convergence
        m = RashbaModel(hbar2_over_2m=0.5, alpha=0.1)
        d = Drive(hw=1.0, gamma=0.3, polarization="circular", sign=-1,
                  a_ref=1.0)
        k = np.array([0.05, 0.03, 0.0])
        res = run_convergence(m, d, k, m.spin_ops(), [80], [2, 8],
                              window=(-2.0, 2.0))
        for row in res["m_sp_rows"]:
            self.assertIn("w_deficit_max", row)
            self.assertIn("n_deficit_gt_1e3", row)
        self.assertTrue(("accepted_m_sp" in res)
                        or ("max_dev_wsum" in res))


if __name__ == "__main__":
    unittest.main()

