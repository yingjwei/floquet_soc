# -*- coding: utf-8 -*-
"""Result-correctness regressions (commercial-hardening pass).

* A1: ``floquet_points`` used to report ``overlap`` against ``evec[j]`` (the
  j-th *row*) while selecting ``j`` with the j-th *column*.  Selection was
  right, the diagnostic value was wrong.
* A3: the Peierls Bessel series was cut at ``pmax`` with no check on the
  discarded weight; a Hermiticity check cannot see it, so we assert the
  explicit :class:`TruncationWarning` fires.
"""
from __future__ import annotations

import os
import sys
import unittest
import warnings

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.benchmarks.rashba import RashbaModel              # noqa: E402
from floquet_soc.core.spin_operator import SPIN_HALF              # noqa: E402
from floquet_soc.observables.floquet_points import run_points     # noqa: E402
from floquet_soc.observables.floquet_points import make_drive     # noqa: E402
from floquet_soc.observables.floquet_spin import (                # noqa: E402
    floquet_spin_observables)
from floquet_soc.drives.peierls import (TruncationWarning,        # noqa: E402
                                        peierls_fourier_components)
from floquet_soc.core.wannier import WannierModel                  # noqa: E402


class TestFloquetPointsOverlap(unittest.TestCase):
    """A1: the reported overlap must be the column overlap at the argmax."""

    def test_overlap_matches_the_column_selected_by_argmax(self):
        model = RashbaModel(hbar2_over_2m=0.5, alpha=0.1)
        spin = lambda k=None: SPIN_HALF                       # noqa: E731
        cfg = {"num_electrons": 1}
        k = np.array([0.05, 0.0, 0.0])
        hw, gamma, pol = 1.0, 0.05, "circ+"
        n_step, m_sp = 12, 2
        res = run_points(model, spin, cfg, [("P", k)], hw, [gamma], [pol],
                         n_step=n_step, m_sp=m_sp)
        targets = res["points"]["P"]["cases"][f"g{gamma}_{pol}"]["targets"]

        # independent recomputation of the *correct* column overlap
        drive = make_drive(hw, gamma, pol, getattr(model, "a_ref", 1.0))
        obs = floquet_spin_observables(model, drive, k, spin(k),
                                       n_step=n_step, m_sp=m_sp)
        E, U = np.linalg.eigh(model.Hk(k))
        for b, name in res["band_names"].items():
            match = U[:, b].conj() @ obs["evec"]
            j = int(np.argmax(np.abs(match)))
            column_ov = float(np.abs(match[j]))
            self.assertAlmostEqual(targets[name]["overlap"], column_ov,
                                   places=12,
                                   msg=f"band {name}: overlap must be the "
                                       f"column overlap at j={j}")
            self.assertEqual(targets[name]["band"], b)


def _toy_wannier():
    win = dict(num_wann=2, num_bands=2, spinors=True,
               lattice=np.eye(3) * 4.0, mp_grid=np.array([2, 2, 1]),
               kpoints=np.zeros((1, 3)))
    h0 = np.array([[0.0, 0.25], [0.25, 0.3]], complex)
    h1 = np.array([[0.1, 0.2j], [-0.2j, -0.1]], complex)
    hr = dict(num_wann=2, nrpts=3,
              R=np.array([(0, 0, 0), (1, 0, 0), (-1, 0, 0)]),
              H=np.asarray([h0, h1, h1.conj().T], complex),
              ndegen=np.ones(3, int))
    return WannierModel(win, hr, tau=np.zeros((2, 3)))


class TestPeierlsTruncationGuard(unittest.TestCase):
    """A3: a large dimensionless argument g must trip the explicit warning."""

    def test_large_g_warns(self):
        model = _toy_wannier()
        drive = make_drive(hw=1.0, gamma=10.0, pol="circ+", a_ref=4.0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            peierls_fourier_components(model, drive, [0.0, 0.0, 0.0], pmax=3)
        self.assertTrue(any(issubclass(w.category, TruncationWarning)
                            for w in caught),
                        "a huge g_max must emit TruncationWarning")

    def test_small_g_is_silent(self):
        model = _toy_wannier()
        drive = make_drive(hw=1.0, gamma=0.001, pol="circ+", a_ref=4.0)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            peierls_fourier_components(model, drive, [0.0, 0.0, 0.0], pmax=3)
        self.assertFalse(any(issubclass(w.category, TruncationWarning)
                             for w in caught))

    def test_warning_is_deduplicated_per_process(self):
        import floquet_soc.drives.peierls as peierls_mod
        peierls_mod._warned_gmax = 0.0          # reset module state
        model = _toy_wannier()
        drive = make_drive(hw=1.0, gamma=10.0, pol="circ+", a_ref=4.0)
        with warnings.catch_warnings(record=True) as first:
            warnings.simplefilter("always")
            peierls_fourier_components(model, drive, [0.0, 0.0, 0.0], pmax=3)
        with warnings.catch_warnings(record=True) as second:
            warnings.simplefilter("always")
            peierls_fourier_components(model, drive, [0.0, 0.0, 0.0], pmax=3)
        self.assertEqual(
            len([w for w in first
                 if issubclass(w.category, TruncationWarning)]), 1)
        self.assertEqual(
            len([w for w in second
                 if issubclass(w.category, TruncationWarning)]), 0,
            "the same g_max must not re-warn on every k point")


if __name__ == "__main__":
    unittest.main()
