# -*- coding: utf-8 -*-
"""Batch E: gap-edge band names from one source, small models degrade to NaN.

Regression: ``target_bands`` returned only indices, so callers re-labelled with
``["VBM-1","VBM","CBM","CBM+1"][:len(bands)]``.  On a 4-band model with
``num_electrons=1`` that produced three bands named VBM-1/VBM/CBM, and
``soc_scan`` then looked up the hard-coded ``"CBM+1"`` -> KeyError.  A 2-band
model crashed earlier, at ``E[ne+1]`` (IndexError).
"""
from __future__ import annotations

import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.observables.floquet_points import target_bands   # noqa: E402


class _ToyModel:
    """nw-band free-electron model with a k-independent spin structure."""

    kind = "analytic"

    def __init__(self, nw=4, gap=0.2):
        self.nw = nw
        self.a_ref = 1.0
        self.lattice = np.eye(3) * 4.0
        self._gap = gap

    def Hk(self, k):
        k = np.asarray(k, float).reshape(3)
        e = np.array([-self._gap / 2 - 0.1 * i for i in range(self.nw // 2)]
                     + [self._gap / 2 + 0.1 * i
                        for i in range(self.nw - self.nw // 2)])
        return np.diag(e).astype(complex) + 0.01 * np.eye(self.nw)


def _spin(nw):
    """nw-generic spin operators: S_z diagonal, S_x between bands 0 and 1."""
    s = np.zeros((3, nw, nw), complex)
    for i in range(nw):
        s[2, i, i] = 0.5 if i % 2 == 0 else -0.5
    if nw >= 2:
        s[0, 0, 1] = s[0, 1, 0] = 0.5
    if nw >= 4:
        s[1, 2, 3] = s[1, 3, 2] = 0.5j
    return lambda k=None: s


class TestTargetBands(unittest.TestCase):

    def test_canonical_quadruple(self):
        bands, names = target_bands({"num_electrons": 2}, _ToyModel(8))
        self.assertEqual(bands, [0, 1, 2, 3])
        self.assertEqual([names[b] for b in bands],
                         ["VBM-1", "VBM", "CBM", "CBM+1"])

    def test_small_model_keeps_honest_names(self):
        # ne=1 -> band 0 IS the VBM; the old code shifted the labels positionally
        # and renamed band 0 to "VBM-1", which is what made "CBM+1" vanish.
        bands, names = target_bands({"num_electrons": 1}, _ToyModel(4))
        self.assertEqual(bands, [0, 1, 2])
        self.assertEqual([names[b] for b in bands], ["VBM", "CBM", "CBM+1"])

    def test_two_band_model_reports_only_what_exists(self):
        bands, names = target_bands({"num_electrons": 1}, _ToyModel(2))
        self.assertEqual(bands, [0, 1])
        self.assertEqual([names[b] for b in bands], ["VBM", "CBM"])

    def test_no_num_electrons_uses_band_labels(self):
        bands, names = target_bands({}, _ToyModel(3))
        self.assertEqual(bands, [0, 1, 2])
        self.assertEqual([names[b] for b in bands], ["band0", "band1", "band2"])


class TestSocScanSmallModel(unittest.TestCase):

    def _run(self, nw, ne):
        from floquet_soc.observables.soc_scan import run_scan
        cfg = {"num_electrons": ne}
        return run_scan(_ToyModel(nw), _spin(nw), cfg, hw=0.1,
                        gammas=[0.0], pols=["circ+"],
                        points="0.1,0.0,0.0;0.0,0.0,0.0",
                        n_step=12, m_sp=4, pst_n=3)

    def test_three_band_ne1_finds_cbm_plus_one(self):
        """The old positional relabelling called band 2 "CBM" and lost CBM+1,
        so this exact case raised KeyError.  With absolute naming band 2 is
        CBM+1 and the conduction pair resolves."""
        res = self._run(nw=3, ne=1)
        row = res["delta_spin"][0]
        self.assertFalse(np.isnan(row["cond_F"]),
                         "CBM/CBM+1 must resolve on a 3-band, ne=1 model")
        self.assertTrue(np.isnan(row["val_F"]),
                        "ne=1 has no VBM-1, so the valence pair is undefined")

    def test_two_band_ne1_reports_nan_not_keyerror(self):
        """CBM+1 = band 2 does not exist on a 2-band model."""
        res = self._run(nw=2, ne=1)
        row = res["delta_spin"][0]
        self.assertTrue(np.isnan(row["cond_F"]),
                        "missing CBM+1 must yield NaN, not KeyError")
        self.assertTrue(np.isnan(row["cond_static"]) or np.isnan(row["delta_cond"]))

    def test_two_band_ne1_does_not_indexerror(self):
        res = self._run(nw=2, ne=1)
        row = res["delta_spin"][0]
        self.assertTrue(np.isnan(row["cond_F"]) or np.isnan(row["cond_static"]))
        self.assertTrue(np.isnan(row["delta_cond"]))


if __name__ == "__main__":
    unittest.main()
