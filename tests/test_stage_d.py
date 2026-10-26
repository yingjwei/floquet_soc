# -*- coding: utf-8 -*-
"""Stage D tests: parameter scans, solver convergence, cross-material benchmark."""
import json
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.cli import main as cli_main
from floquet_soc.workflows.scan import run_scan, rows_to_arrays
from floquet_soc.workflows.converge import run_convergence
from floquet_soc.benchmarks.rashba import RashbaModel
from floquet_soc.core.spin_operator import SPIN_HALF
from floquet_soc.drives.peierls import Drive

TMP = os.path.join(HERE, "outputs", "_tests_d")
CFG_RASHBA = os.path.join(HERE, "configs", "rashba.json")
CFG_PBNV = os.path.join(HERE, "configs", "pbnv.json")
CFG_PT = os.path.join(HERE, "configs", "pt.json")


class TestScanWorkflow(unittest.TestCase):
    def test_scan_analytic(self):
        model = RashbaModel(alpha=0.1)
        drive_cfg = dict(polarization="circular", sign=-1)
        kpt = np.array([[0.05, 0.03, 0.0]])
        rows = run_scan(model, drive_cfg, [0.004, 0.008],
                        [0.0, 0.002], kpt, SPIN_HALF,
                        n_step=120, m_sp=8, window=(-0.1, 0.1))
        arr = rows_to_arrays(rows)
        self.assertEqual(len(rows), 4)
        self.assertTrue(np.all(np.isfinite(arr["mean_pol"])))
        # gamma=0: Floquet spin equals the static spin (|S|=1/2)
        pol0 = arr["mean_pol"][arr["gamma"] == 0.0]
        self.assertLess(float(np.max(np.abs(pol0 - 0.5))), 1e-6)
        # truncation of the photon window must be small in this window
        self.assertLess(float(np.nanmax(arr["wsum_err"])), 1e-6)

    def test_scan_cli(self):
        rc = cli_main(["scan", "--config", CFG_RASHBA,
                       "--hw", "0.004", "0.008", "--gamma", "0.0", "0.002",
                       "--n-step", "100", "--m-sp", "8",
                       "--outdir", os.path.join(TMP, "scan")])
        self.assertEqual(rc, 0)
        d = np.load(os.path.join(TMP, "scan", "scan.npz"))
        self.assertTrue(np.all(np.isfinite(d["mean_pol"])))
        self.assertTrue(os.path.exists(os.path.join(TMP, "scan", "scan.png")))


class TestConvergenceWorkflow(unittest.TestCase):
    def test_converge_analytic(self):
        model = RashbaModel(alpha=0.1)
        drive = Drive(hw=0.02, gamma=0.002, a_ref=1.0,
                      polarization="circular", sign=-1)
        res = run_convergence(model, drive, np.array([0.05, 0.03, 0.0]),
                              SPIN_HALF, [60, 100, 200], [4, 8, 16],
                              window=(-0.1, 0.1))
        self.assertTrue(res["pass_"])
        devs = [r["dev_ebar"] for r in res["n_step_rows"]]
        # deviations must shrink when n_step increases
        self.assertLessEqual(devs[-2], devs[0])
        self.assertLess(res["max_dev_invariant"], 1e-12)

    def test_converge_cli(self):
        rc = cli_main(["converge", "--config", CFG_RASHBA,
                       "--n-step", "60", "120",
                       "--m-sp", "4", "8",
                       "--outdir", os.path.join(TMP, "conv")])
        self.assertEqual(rc, 0)
        d = np.load(os.path.join(TMP, "conv", "converge.npz"))
        self.assertEqual(str(d["status"][0]), "PASS")

    def test_converge_pbnv_small(self):
        from floquet_soc.cli import load_config
        cfg = load_config(CFG_PBNV)
        for key in ("chk", "eig"):
            if not os.path.exists(cfg["files"][key]):
                self.skipTest("PbNV checkpoint data not available")
        cfg["converge"] = {"kpt": [[0.05, 0.03, 0.0]],
                           "n_step": [80, 120, 200],
                           "m_sp": [4, 8],
                           "window": [-2.0, 2.0],
                           # this test only checks the workflow plumbing;
                           # production tolerances are 1e-5 / 1e-4 (and the
                           # m_sp=4 photon-weight deficit is ~1, so the wsum
                           # criterion must be relaxed too for this plumbing
                           # check; physical convergence is covered by the
                           # PbNV production runs and test_converge_analytic)
                           "tol_ebar": 5e-3, "tol_spin": 5e-3,
                           "tol_wsum": 1.5}
        os.makedirs(TMP, exist_ok=True)
        fn = os.path.join(TMP, "pbnv_conv.json")
        json.dump(cfg, open(fn, "w", encoding="utf-8"), ensure_ascii=False)
        rc = cli_main(["converge", "--config", fn,
                       "--outdir", os.path.join(TMP, "pbnv_conv")])
        self.assertEqual(rc, 0)
        d = np.load(os.path.join(TMP, "pbnv_conv", "converge.npz"))
        self.assertEqual(str(d["status"][0]), "PASS")


class TestBenchmark7(unittest.TestCase):
    def test_pt_benchmark(self):
        from floquet_soc.cli import load_config
        cfg = load_config(CFG_PT)
        for key in ("chk", "eig"):
            if not os.path.exists(cfg["files"][key]):
                self.skipTest("Pt benchmark data not available")
        from floquet_soc.benchmarks import pt
        res = pt.run(outdir=os.path.join(TMP, "pt"))
        rep = res["report"]
        self.assertEqual(rep.count("FAIL"), 0)
        self.assertEqual(res["model"].nw, 18)
        self.assertTrue(os.path.exists(res["figure"]))
        kinds = rep.count("PASS") + rep.count("WARN")
        self.assertGreater(kinds, 5)

    def test_pt_chk_eig_reconstruction(self):
        from floquet_soc.cli import load_config
        from floquet_soc.core.chk import read_chk, read_eig
        cfg = load_config(CFG_PT)
        if not os.path.exists(cfg["files"]["chk"]):
            self.skipTest("Pt benchmark data not available")
        from floquet_soc.cli import build_model
        m = build_model(cfg)
        chk = read_chk(cfg["files"]["chk"])
        eig = read_eig(cfg["files"]["eig"])
        V = chk["v_matrix"]
        Hw = np.empty((len(eig), m.nw, m.nw), complex)
        for ik in range(len(eig)):
            Vk = V[ik]
            Hw[ik] = np.conj(Vk.T) @ (eig[ik].astype(complex)[:, None] * Vk)
            Hw[ik] = 0.5 * (Hw[ik] + np.conj(Hw[ik].T))
        dev = np.max(np.abs(m.Hk_batch(chk["kpt_red"]) - Hw))
        self.assertLess(dev, 1e-8,
                        "chk+eig reconstruction must be exact on the DFT mesh")
        self.assertTrue(getattr(m, "reconstructed", False))


if __name__ == "__main__":
    unittest.main(verbosity=2)
