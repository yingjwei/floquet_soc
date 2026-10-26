# -*- coding: utf-8 -*-
"""Tests for the chk-rotated Wannier spin operator and PbNV spin texture."""
import json
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.cli import load_config, build_model, main as cli_main
from floquet_soc.core.chk import (read_chk, read_eig, checkpoint_checks,
                                  hamiltonian_consistency)
from floquet_soc.core.spin_operator import read_wannier_spn
from floquet_soc.core.spin_wannier import build_wannier_spin_operator

CFG = os.path.join(HERE, "configs", "pbnv.json")
TMP = os.path.join(HERE, "outputs", "_tests_spin")
# Optional data-dependent check against a *formatted* (text) wannier90 .spn
# file.  This fixture is 192 bands x 36 k-points of the real PbNV run and is
# not shipped; point FLOQUET_SOC_FORMATTED_SPN at a copy to enable the test.
FORMATTED_SPN = os.environ.get("FLOQUET_SOC_FORMATTED_SPN", "").strip()


class TestCheckpoint(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config(CFG)
        files = cls.cfg["files"]
        for key in ("win", "hr", "chk", "eig"):
            if not os.path.exists(files[key]):
                raise unittest.SkipTest(f"missing {key}: {files[key]}")
        cls.model = build_model(cls.cfg)
        cls.chk = read_chk(files["chk"])
        cls.eig = read_eig(files["eig"])

    def test_checkpoint_structure(self):
        chk = self.chk
        self.assertEqual(chk["num_bands"], 192)
        self.assertEqual(chk["num_wann"], 156)
        self.assertEqual(chk["num_kpts"], 36)
        self.assertEqual(list(chk["mp_grid"]), [6, 6, 1])
        self.assertEqual(chk["v_matrix"].shape, (36, 192, 156))
        c = checkpoint_checks(chk)
        self.assertLess(c["u_unitarity"], 1e-8)
        self.assertLess(c["v_isometry"], 1e-8)
        self.assertEqual(c["window_inconsistent"], 0)

    def test_hamiltonian_consistency(self):
        hc = hamiltonian_consistency(self.chk, self.eig, self.model)
        self.assertLess(hc["max_H_diff"], 1e-3,
                        "H_W(k) must reproduce H_hr(k)")
        self.assertLess(hc["max_eig_diff"], 1e-3)

    def test_rotation_and_interpolation(self):
        cfg = self.cfg
        op, info = build_wannier_spin_operator(
            cfg["spin"]["file"], cfg["spin"]["chk"], cfg["spin"]["eig"],
            model=self.model)
        checks = dict((n, (p, d)) for n, p, d in info["checks"])
        self.assertTrue(checks["spin_wannier_hermiticity"][0])
        self.assertTrue(checks["spin_wannier_interpolation"][0])
        self.assertTrue(checks["checkpoint_hamiltonian_consistency"][0])
        # TRS quality is data dependent: only require a finite value here
        self.assertIn("spin_wannier_trs", checks)
        self.assertTrue(np.isfinite(info["metrics"]["eig_max"]))

    def test_read_wannier_spn_dispatch(self):
        if not FORMATTED_SPN or not os.path.exists(FORMATTED_SPN):
            self.skipTest("formatted spn fixture not present; set "
                          "FLOQUET_SOC_FORMATTED_SPN to a wannier90 "
                          "*formatted* .spn file for the PbNV run")
        spn = read_wannier_spn(FORMATTED_SPN)
        self.assertEqual(spn["format"], "formatted-spn")
        self.assertEqual(spn["num_bands"], 192)
        self.assertEqual(spn["num_kpts"], 36)
        self.assertEqual(spn["matrices"].shape, (36, 3, 192, 192))
        old = read_wannier_spn(self.cfg["spin"]["file"])
        self.assertEqual(old["format"], "vasp-binary-spn")


class TestPbNVSpinTexture(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config(CFG)
        for key in ("chk", "eig"):
            if not os.path.exists(cls.cfg["files"][key]):
                raise unittest.SkipTest("checkpoint data not available")

    def test_validate_passes(self):
        rc = cli_main(["validate", "--config", CFG,
                       "--outdir", os.path.join(TMP, "validate")])
        self.assertEqual(rc, 0)
        with open(os.path.join(TMP, "validate", "validate.md"),
                  encoding="utf-8") as f:
            text = f.read()
        self.assertIn("spin_operator_wannier_basis", text)
        self.assertNotIn("**BLOCKED** |", text)

    def test_static_spin_texture(self):
        out = os.path.join(TMP, "spin")
        cfg = load_config(CFG)
        cfg["kmesh"] = {"n": 9, "kmax": 0.3}
        cfg["jobs"] = 1
        tmp = os.path.join(TMP, "cfg_static.json")
        json.dump(cfg, open(tmp, "w", encoding="utf-8"))
        rc = cli_main(["spin_static", "--config", tmp, "--outdir", out])
        self.assertEqual(rc, 0)
        npz = np.load(os.path.join(out, "spin_static.npz"))
        self.assertEqual(npz["spins"].shape[2], 3)
        self.assertEqual(npz["spins"].shape[1], 156)
        self.assertTrue(np.all(np.isfinite(npz["spins"])))

    def test_spin_floquet_pbnv(self):
        """Stage-B Floquet spin texture on the real material (small kmesh)."""
        # load_config() resolves `root`/file paths to absolute paths; dump the
        # *resolved* config so copying it into a scratch dir keeps working
        cfg = load_config(CFG)
        cfg["root"] = cfg["_root"]
        cfg["kmesh"] = {"n": 3, "kmax": 0.10}
        cfg["kmesh_spin_floquet"] = {"n": 3, "kmax": 0.10}
        cfg["jobs"] = 1
        cfg["numerics"] = {"n_step": 80, "m_sp": 8}
        tmp_cfg = os.path.join(TMP, "pbnv_small.json")
        os.makedirs(TMP, exist_ok=True)
        json.dump(cfg, open(tmp_cfg, "w", encoding="utf-8"),
                  ensure_ascii=False)
        out = os.path.join(TMP, "spin_floquet")
        rc = cli_main(["spin_floquet", "--config", tmp_cfg,
                       "--outdir", out])
        self.assertEqual(rc, 0)
        npz = np.load(os.path.join(out, "spin_floquet.npz"))
        self.assertEqual(npz["S_bar"].shape, (9, 156, 3))
        self.assertTrue(np.all(np.isfinite(npz["S_bar"])))


if __name__ == "__main__":
    unittest.main(verbosity=2)
