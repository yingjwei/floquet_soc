# -*- coding: utf-8 -*-
"""B/B8: the H(k) Hermiticity guard and the folded/expanded split.

Two distinct failure modes are covered:

* folded ``_hr.dat`` with an incomplete R/-R set  -> the H(R) residual is
  *undefined* (``None``), never a false 0.0;
* one-sided ``wsvec.dat`` block set -> ``H(k)`` is genuinely non-Hermitian and
  must be refused, while ``hermiticity_errors()`` still reports the residual
  (it must not raise before it can be measured).
"""
from __future__ import annotations

import os
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.core.wannier import (WannierModel, HermiticityError,   # noqa: E402
                                     assert_hermitian)


def _win(nw=2):
    return dict(num_wann=nw, num_bands=nw, spinors=True,
                lattice=np.eye(3) * 4.0, mp_grid=np.array([2, 2, 1]),
                kpoints=np.zeros((1, 3)))


def _hr(Rlist, blocks, ndegen=None):
    Rlist = np.asarray(Rlist, int)
    return dict(num_wann=blocks[0].shape[0], nrpts=len(Rlist), R=Rlist,
                H=np.asarray(blocks, complex),
                ndegen=np.ones(len(Rlist), int) if ndegen is None
                else np.asarray(ndegen, int))


def _toy(nw=2):
    """On-site + symmetric hopping on R=(+-1,0,0): Hermitian everywhere."""
    R = [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]
    h0 = np.array([[0.0, 0.25], [0.25, 0.3]], complex)
    h1 = np.array([[0.1, 0.2j], [-0.2j, -0.1]], complex)
    return _hr(R, [h0, h1, h1.conj().T])


class TestHermiticityGuard(unittest.TestCase):

    def test_good_model_passes(self):
        m = WannierModel(_win(), _toy(), tau=np.zeros((2, 3)))
        self.assertLess(assert_hermitian(m.Hk([0.13, 0.07, 0.0])), 1e-12)
        err = m.hermiticity_errors()
        self.assertEqual(err["n_unpaired_R"], 0)
        self.assertIsNotNone(err["HR_hermiticity_folded"])
        self.assertLess(err["Hk_hermiticity_expanded"], 1e-12)

    def test_unpaired_R_reports_none_not_zero(self):
        R = [(0, 0, 0), (1, 0, 0)]                    # -R is missing
        h0 = np.array([[0.0, 0.0], [0.0, 0.3]], complex)
        h1 = np.array([[0.1, 0.2j], [-0.2j, -0.1]], complex)
        m = WannierModel(_win(), _hr(R, [h0, h1]), tau=np.zeros((2, 3)))
        err = m.hermiticity_errors()
        self.assertEqual(err["n_unpaired_R"], 1)      # (1,0,0) has no -R
        self.assertIsNone(err["HR_hermiticity_folded"],
                          "unpaired R must not be reported as residual 0.0")
        # H(k) is the fallback judge: it must see the asymmetry
        self.assertGreater(err["Hk_hermiticity_expanded"], 1e-3)
        with self.assertRaises(HermiticityError):
            m.Hk([0.137, 0.083, 0.021])

    def test_one_sided_wsvec_is_refused(self):
        """A wsvec block set that fills only (i,j) and never (j,i)."""
        win = _win()
        win["num_wann"] = 2
        hr = _toy()
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "wannier90_wsvec.dat")
            with open(path, "w", encoding="utf-8") as f:
                f.write(" use_ws_distance=.true.\n")
                # R=(0,0,0), pair (1,2) -> shift +(1,0,0) only; the conjugate
                # pair (2,1) is deliberately absent
                f.write("   0   0   0   1   2\n")
                f.write("1\n")
                f.write("   1   0   0\n")
            m = WannierModel(win, hr, tau=np.zeros((2, 3)), wsvec_path=path)
            self.assertIsNotNone(m.HR_full, "wsvec must have been applied")
            err = m.hermiticity_errors()
            self.assertGreater(err["Hk_hermiticity_expanded"], 1e-3)
            with self.assertRaises(HermiticityError):
                m.Hk([0.137, 0.083, 0.021])
            with self.assertRaises(HermiticityError):
                m.static_bands(np.array([[0.0, 0.0, 0.0]]))

    def test_offmesh_probe_is_needed(self):
        """On-mesh probes can be blind; the default probe set must not be."""
        err = WannierModel(_win(), _toy(), tau=np.zeros((2, 3))).hermiticity_errors()
        kt = err["ktest"]
        onmesh = kt[(np.abs(np.mod(kt * np.array([2, 2, 1]), 1.0)) < 1e-12)
                    .all(axis=1)]
        self.assertGreater(len(kt) - len(onmesh), 0,
                           "at least one probe must be off the mp_grid mesh")


if __name__ == "__main__":
    unittest.main()
