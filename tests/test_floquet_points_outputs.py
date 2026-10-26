# -*- coding: utf-8 -*-
"""H0: ``floquet_points.save_and_plot`` must return an *honest* figure map.

Regression: the return value used to be re-bound unconditionally, which
(a) dropped the ``valley`` key written by the ``if asym:`` block and
(b) advertised ``floquet_points_valley.png`` even when it was never written.
"""
from __future__ import annotations

import os
import tempfile
import unittest

import numpy as np


def _minimal_res():
    case = dict(
        gamma=0.0, pol="circ+",
        targets={"VBM-1": dict(Ebar=0.01, S_bar=np.zeros(3), dS_norm=0.1),
                 "VBM": dict(Ebar=0.00, S_bar=np.zeros(3), dS_norm=0.2)})
    return dict(hw=0.4, bands=[0, 1], band_names={0: "VBM-1", 1: "VBM"},
                points={"K": {"cases": {"g0_circ+": case}}})


class TestSaveAndPlotContract(unittest.TestCase):

    def setUp(self):
        self.outdir = tempfile.mkdtemp(prefix="fp_save_")

    def test_no_valley_key_without_asymmetry(self):
        from floquet_soc.observables.floquet_points import save_and_plot
        out = save_and_plot(_minimal_res(), None, self.outdir)
        self.assertNotIn("valley", out,
                         "no valley figure was written, so it must not be "
                         "advertised")
        self.assertTrue(os.path.exists(out["json"]))
        self.assertTrue(os.path.exists(out["spin"]))

    def test_valley_key_only_when_written(self):
        from floquet_soc.observables.floquet_points import save_and_plot
        asym = {"K|Kp": {"K|Kp": {"VBM-1": dict(dEbar=0.11),
                                  "VBM": dict(dEbar=0.07)}}}
        out = save_and_plot(_minimal_res(), asym, self.outdir)
        self.assertIn("valley", out)
        self.assertTrue(os.path.exists(out["valley"]),
                        "advertised valley figure must exist")

    def test_empty_asym_dict_behaves_like_none(self):
        from floquet_soc.observables.floquet_points import save_and_plot
        out = save_and_plot(_minimal_res(), {}, self.outdir)
        self.assertNotIn("valley", out)


if __name__ == "__main__":
    unittest.main()
