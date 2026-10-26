# -*- coding: utf-8 -*-
"""Smoke tests for the new path/scan workflows."""
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.benchmarks.rashba import RashbaModel       # noqa: E402
from floquet_soc.drives.peierls import Drive                # noqa: E402
from floquet_soc.observables.floquet_bands import (         # noqa: E402
    compute_path, save_and_plot)


class TestFloquetBands(unittest.TestCase):
    def test_compute_and_save(self):
        model = RashbaModel(hbar2_over_2m=0.5, alpha=0.1)
        drive = Drive(hw=0.4, gamma=0.0, a_ref=1.0)
        kpts = np.array([[0.0, 0.0, 0.0], [0.05, 0.0, 0.0],
                         [0.1, 0.0, 0.0]])
        res = compute_path(model, drive, kpts, n_step=20, m_sp=2,
                           ne=None)
        self.assertEqual(res["E"].shape, (3, 2))
        self.assertEqual(res["A"].shape[0], 3)
        self.assertEqual(res["A0"].shape, res["A"].shape)
        self.assertEqual(res["E_grid"].ndim, 1)
        outdir = tempfile.mkdtemp(prefix="floquet_soc_fb_")
        self.addCleanup(shutil.rmtree, outdir, True)
        out = save_and_plot(res, ["Γ", "K"], outdir, prefix="fb")
        self.assertTrue(os.path.exists(out["npz"]))
        self.assertTrue(os.path.exists(out["figure"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
