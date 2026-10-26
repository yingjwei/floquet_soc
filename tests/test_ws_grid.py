# -*- coding: utf-8 -*-
"""Regression tests for the Wannier90 WS supercell (odd/even mp_grid) fix.

v0.1.5-rc1: ``ws_supercell`` used ``np.arange(-N // 2, ...)`` which is wrong
for odd N (Python floor division).  For the MoSSe 15x15x1 mesh this produced
512 R points with ndegen 2/4/8 instead of 225 points with ndegen 1, making the
off-grid spin-operator interpolation non-Hermitian.
"""
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.core.wannier import ws_supercell, read_wsvec, WannierModel  # noqa: E402
from floquet_soc.core.spin_operator import load_spin_for_model    # noqa: E402
from floquet_soc.core.chk import read_eig                         # noqa: E402

MOSSE = os.path.join(os.path.dirname(HERE), "MoSSe")
MOSSE_FILES = dict(
    win=os.path.join(MOSSE, "wannier90.win"),
    hr=os.path.join(MOSSE, "wannier90_hr.dat"),
    chk=os.path.join(MOSSE, "wannier90.chk"),
    eig=os.path.join(MOSSE, "wannier90.eig"),
    spn=os.path.join(MOSSE, "wannier90.spn"),
    centres=os.path.join(MOSSE, "wannier90_centres.xyz"),
    wsvec=os.path.join(MOSSE, "wannier90_wsvec.dat"),
)
HAVE_MOSSE = all(os.path.exists(v) for v in MOSSE_FILES.values())


def _lattice():
    return np.array([[3.2514695, 0.0, 0.0],
                     [-1.6257348, 2.8158552, 0.0],
                     [0.0, 0.0, 24.1702538]])


def _set_of_rows(R):
    return {tuple(int(x) for x in row) for row in np.asarray(R)}


class TestWsSupercell(unittest.TestCase):
    def test_odd_grid_15x15(self):
        R, nd = ws_supercell([15, 15, 1], _lattice())
        self.assertEqual(R.shape, (225, 3))
        self.assertEqual(R[:, 0].min(), -7)
        self.assertEqual(R[:, 0].max(), 7)
        self.assertEqual(R[:, 1].min(), -7)
        self.assertEqual(R[:, 1].max(), 7)
        self.assertEqual(R[:, 2].min(), 0)
        self.assertEqual(R[:, 2].max(), 0)
        self.assertTrue(np.all(nd == 1))
        rows = _set_of_rows(R)
        for r in rows:                                   # closed under -R
            self.assertIn((-r[0], -r[1], -r[2]), rows)

    def test_even_grid_16x16(self):
        R, nd = ws_supercell([16, 16, 1], _lattice())
        self.assertEqual(R.shape, (17 * 17, 3))
        self.assertEqual(R[:, 0].min(), -8)
        self.assertEqual(R[:, 0].max(), 8)
        self.assertTrue(np.all(nd >= 1))
        self.assertEqual(int(nd.max()), 4)
        self.assertAlmostEqual(float(np.sum(1.0 / nd)), 16.0 * 16.0, places=10)
        rows = _set_of_rows(R)
        for r in rows:
            self.assertIn((-r[0], -r[1], -r[2]), rows)

    def test_single_kpoint_dim(self):
        R, nd = ws_supercell([1, 1, 1], _lattice())
        self.assertEqual(R.shape, (1, 3))
        self.assertEqual(tuple(int(x) for x in R[0]), (0, 0, 0))
        self.assertEqual(int(nd[0]), 1)


@unittest.skipUnless(HAVE_MOSSE, "MoSSe data not available")
class TestMoSSeWsSet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = WannierModel.from_files(MOSSE_FILES["win"],
                                            MOSSE_FILES["hr"],
                                            MOSSE_FILES["centres"])

    def test_wsvec_matches_hr_set(self):
        wsvec = read_wsvec(MOSSE_FILES["wsvec"], keep_vectors=True)
        hr_R = _set_of_rows(self.model.Rlist)
        self.assertEqual(_set_of_rows(wsvec["R"]), hr_R)
        nw = self.model.nw
        self.assertEqual(int(wsvec["blocks"]), len(hr_R) * nw * nw)
        grid = np.asarray(self.model.win["mp_grid"], int)
        for vs in wsvec["vectors"]:
            self.assertTrue(np.all(np.mod(vs, grid) == 0))
        self.assertLessEqual(int(wsvec["n_wsvec"].max()), 8)

    def test_Hk_offgrid_hermiticity(self):
        ktest = np.array([[0.01, 0.02, 0.0], [0.3, 0.1, 0.0],
                          [-0.2, 0.25, 0.0], [1.0 / 3, 1.0 / 3, 0.0]])
        for k in ktest:
            H = self.model.Hk(k)
            self.assertLess(float(np.max(np.abs(H - H.conj().T))), 1e-10)

    def test_from_chk_eig_reproduces_hr_on_mesh(self):
        # from_chk_eig uses ws_supercell; on the DFT mesh it must reproduce
        # the chk+hr Hamiltonian (this is the regression for the odd-grid fix)
        m2 = WannierModel.from_chk_eig(MOSSE_FILES["chk"],
                                       MOSSE_FILES["eig"],
                                       MOSSE_FILES["win"])
        from floquet_soc.core.chk import read_chk
        chk = read_chk(MOSSE_FILES["chk"])
        dmax = max(float(np.max(np.abs(self.model.Hk(k) - m2.Hk(k))))
                   for k in chk["kpt_red"])
        self.assertLess(dmax, 1e-3)


@unittest.skipUnless(HAVE_MOSSE, "MoSSe data not available")
class TestMoSSeSpinOffgrid(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = WannierModel.from_files(MOSSE_FILES["win"],
                                            MOSSE_FILES["hr"],
                                            MOSSE_FILES["centres"])
        cls.cfg = dict(files=MOSSE_FILES,
                       spin=dict(source="spn", file=MOSSE_FILES["spn"],
                                 chk=MOSSE_FILES["chk"],
                                 eig=MOSSE_FILES["eig"]))

    def test_spin_operator_offgrid_hermitian(self):
        spin = load_spin_for_model(self.cfg, self.model)
        self.assertIn("spin_operator", spin)
        op = spin["spin_operator"]
        self.assertTrue(op.ws_interpolation)
        # the operator must use the same Wannier90 WS R set as H(k)
        self.assertEqual(op.R_source, "wannier90_hr.dat")
        self.assertEqual(_set_of_rows(op.Rvecs),
                         _set_of_rows(self.model.Rlist))
        # 15x15x1 is an odd grid: this was ~1.0 before the v0.1.5 fix
        self.assertLess(op.hermiticity_offgrid(), 1e-10)
        self.assertLess(op.interpolation_error(), 1e-10)


PBNV = os.path.join(os.path.dirname(HERE), "PbNV")
PBNV_FILES = dict(
    win=os.path.join(PBNV, "wannier90.win"),
    hr=os.path.join(PBNV, "wannier90_hr.dat"),
    chk=os.path.join(PBNV, "wannier90.chk"),
    eig=os.path.join(PBNV, "wannier90.eig"),
    spn=os.path.join(PBNV, "wannier90.spn"),
    centres=os.path.join(PBNV, "wannier90_centres.xyz"),
)
HAVE_PBNV = all(os.path.exists(v) for v in PBNV_FILES.values())


@unittest.skipUnless(HAVE_PBNV, "PbNV data not available")
class TestPbNVEvenGridSpin(unittest.TestCase):
    def test_even_grid_spin_interpolation(self):
        # PbNV: 6x6 mesh -> Wannier90 WS set has 49 R vectors with
        # ndegen 1/2/4 (even-N folding), so len(R) != nk must be allowed.
        model = WannierModel.from_files(PBNV_FILES["win"], PBNV_FILES["hr"],
                                        PBNV_FILES["centres"])
        cfg = dict(files=PBNV_FILES,
                   spin=dict(source="spn", file=PBNV_FILES["spn"],
                             chk=PBNV_FILES["chk"], eig=PBNV_FILES["eig"]))
        spin = load_spin_for_model(cfg, model)
        op = spin["spin_operator"]
        self.assertEqual(op.R_source, "wannier90_hr.dat")
        self.assertGreater(len(op.Rvecs), len(op.kmesh))
        self.assertGreaterEqual(int(op.ndegen.max()), 2)
        self.assertLess(op.hermiticity_offgrid(), 1e-10)
        self.assertLess(op.interpolation_error(), 1e-10)


if __name__ == "__main__":
    unittest.main(verbosity=2)
