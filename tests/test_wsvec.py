# -*- coding: utf-8 -*-
"""Regression tests for Wannier90 use_ws_distance (wsvec) interpolation.

v0.1.6-rc1: with use_ws_distance=.true. and write_ndegen_applied absent, the
R list in wannier90_hr.dat is folded; wannier90_wsvec.dat gives pair-dependent
supercell shifts that must be applied before the plain Fourier transform
(Wannier90 ws_apply_ndegen).  The framework now expands H(k), the Wannier-basis
spin operator and the Peierls-driven Hamiltonian consistently.
"""
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.core.wannier import (read_wsvec, WannierModel,
                                      WannierInputError)        # noqa: E402
from floquet_soc.core.spin_operator import load_spin_for_model    # noqa: E402
from floquet_soc.drives.peierls import peierls_Hk               # noqa: E402

MOSSE = os.path.join(os.path.dirname(HERE), "MoSSe")
FILES = dict(
    win=os.path.join(MOSSE, "wannier90.win"),
    hr=os.path.join(MOSSE, "wannier90_hr.dat"),
    chk=os.path.join(MOSSE, "wannier90.chk"),
    eig=os.path.join(MOSSE, "wannier90.eig"),
    spn=os.path.join(MOSSE, "wannier90.spn"),
    centres=os.path.join(MOSSE, "wannier90_centres.xyz"),
    wsvec=os.path.join(MOSSE, "wannier90_wsvec.dat"),
)
HAVE = all(os.path.exists(v) for v in FILES.values())


@unittest.skipUnless(HAVE, "MoSSe wsvec data not available")
class TestWannierWSVec(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = WannierModel.from_files(
            FILES["win"], FILES["hr"], FILES["centres"],
            wsvec_path=FILES["wsvec"])
        cls.model_plain = WannierModel.from_files(
            FILES["win"], FILES["hr"], FILES["centres"])

    def test_wsvec_flags_and_statistics(self):
        wv = read_wsvec(FILES["wsvec"])
        self.assertTrue(wv["use_ws_distance"])
        self.assertFalse(wv["write_ndegen_applied"])
        self.assertEqual(int(wv["blocks"]), 225 * 34 * 34)
        self.assertEqual(int(wv["n_wsvec"].max()), 3)
        self.assertGreater(int(np.sum(wv["n_wsvec"] > 1)), 0)
        self.assertGreater(int(np.sum(np.any(wv["shifts_flat"] != 0,
                                             axis=1))), 0)

    def test_model_expansion(self):
        self.assertIsNotNone(self.model.HR_full)
        self.assertIsNone(self.model_plain.HR_full)
        self.assertGreater(len(self.model.Rlist_full), len(self.model.Rlist))
        # on-grid: the two conventions are identical (shift is a supercell
        # vector, exp(2*pi*i k.shift)=1 on the mesh)
        kmesh = np.asarray(self.model.win["kpoints"], float)
        dH = max(float(np.max(np.abs(self.model.Hk(k)
                                   - self.model_plain.Hk(k))))
                 for k in kmesh[::3])
        self.assertLess(dH, 1e-10)
        # off-grid: the conventions differ (this is the bug being guarded)
        k = np.array([0.02, 0.012, 0.0])
        dOff = float(np.max(np.abs(self.model.Hk(k)
                                   - self.model_plain.Hk(k))))
        self.assertGreater(dOff, 1e-4)
        H = self.model.Hk(k)
        self.assertLess(float(np.max(np.abs(H - H.conj().T))), 1e-10)

    def test_peierls_uses_same_vectors(self):
        k = np.array([0.05, 0.03, 0.0])
        self.assertLess(float(np.max(np.abs(
            peierls_Hk(self.model, k, np.zeros(3)) - self.model.Hk(k)))),
            1e-12)
        A = np.array([0.0032, -0.0017, 0.0])
        H = peierls_Hk(self.model, k, A)
        self.assertLess(float(np.max(np.abs(H - H.conj().T))), 1e-10)
        # the plain convention gives a different driven Hamiltonian
        Hp = peierls_Hk(self.model_plain, k, A)
        self.assertGreater(float(np.max(np.abs(H - Hp))), 1e-5)

    def test_spin_operator_wsvec(self):
        cfg = dict(files=FILES, spin=dict(source="spn", file=FILES["spn"],
                                         chk=FILES["chk"], eig=FILES["eig"]))
        spin = load_spin_for_model(cfg, self.model)
        op = spin["spin_operator"]
        self.assertEqual(op.R_source, "wannier90_wsvec.dat")
        self.assertTrue(op.full_ws)
        self.assertLess(op.hermiticity_offgrid(), 1e-12)
        self.assertLess(op.interpolation_error(), 1e-10)
        names = [n for n, p, d in spin["checks"]]
        self.assertIn("spin_wannier_interpolation_set", names)


@unittest.skipUnless(HAVE, "MoSSe wsvec data not available")
class TestWSVecPhysicsConvention(unittest.TestCase):
    def _alpha(self, model):
        lat = np.asarray(model.lattice, float)
        B = 2 * np.pi * np.linalg.inv(lat).T
        K = np.array([1.0 / 3.0, 1.0 / 3.0, 0.0]) @ B
        u = K / np.linalg.norm(K)
        ts = np.linspace(2e-3, 3.0e-2, 12)
        kf = np.array([(t * u) @ np.linalg.inv(B) for t in ts])
        E = np.array([np.linalg.eigvalsh(model.Hk(k)) for k in kf])
        dE = E[:, 25] - E[:, 24]
        A = 2 * ts
        return float(np.sum(A * dE) / np.sum(A * A))

    def test_rashba_alpha_convention(self):
        ws = WannierModel.from_files(FILES["win"], FILES["hr"],
                                     FILES["centres"],
                                     wsvec_path=FILES["wsvec"])
        plain = WannierModel.from_files(FILES["win"], FILES["hr"],
                                        FILES["centres"])
        a_ws, a_plain = self._alpha(ws) * 1000, self._alpha(plain) * 1000
        # regression values from the independent audit (v0.1.6-rc1)
        self.assertAlmostEqual(a_ws, 44.3, delta=3.0)
        self.assertAlmostEqual(a_plain, 64.8, delta=3.0)
        self.assertGreater(a_plain - a_ws, 10.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
