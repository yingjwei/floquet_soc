# -*- coding: utf-8 -*-
"""Stage B tests: Floquet solvers, periodic spin texture, spectral functions,
PbNV engine regression against the previously validated implementation.
"""
import json
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.cli import load_config, build_model, kpath_from_cfg, main as cli_main
from floquet_soc.core.spin_operator import SPIN_HALF
from floquet_soc.core.validation import PASS, FAIL, BLOCKED
from floquet_soc.drives.peierls import Drive
from floquet_soc.solvers.evolution import floquet_modes
from floquet_soc.solvers.sambe import sambe_vs_evolution
from floquet_soc.solvers.effective import effective_split_error
from floquet_soc.observables.floquet_spin import (floquet_spin_observables,
                                                  spin_spectral_density)
from floquet_soc.observables.spin_texture import spin_texture
from floquet_soc.benchmarks.rashba import RashbaModel

CFG = os.path.join(HERE, "configs", "pbnv.json")
TMP = os.path.join(HERE, "outputs", "_tests_b")
OLD_ENGINE_REF = os.path.join(HERE, "outputs", "_tests_b",
                              "old_engine_ref.npz")
REF_IDX = [0, 40, 80, 120, 160]


def drive(hw, gamma, a_ref=1.0):
    return Drive(hw=hw, gamma=gamma, a_ref=a_ref,
                 polarization="circular", e1=(1, 0, 0), e2=(0, 1, 0),
                 sign=-1)


class TestFloquetSolvers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = RashbaModel(hbar2_over_2m=0.5, alpha=0.1)
        cls.ops = SPIN_HALF
        cls.k = np.array([0.10, 0.05, 0.0])

    def test_zero_drive_quasienergy(self):
        fl = floquet_modes(self.model, drive(1.0, 0.0), self.k, n_step=200)
        stat = np.linalg.eigvalsh(self.model.Hk(self.k))
        folded = (stat + 0.5) % 1.0 - 0.5
        self.assertLess(float(np.max(np.abs(np.sort(fl["eps"])
                                            - np.sort(folded)))), 1e-10)

    def test_zero_drive_spin_texture(self):
        obs = floquet_spin_observables(self.model, drive(1.0, 0.0), self.k,
                                       self.ops, n_step=200)
        stat = spin_texture(self.model.Hk(self.k), self.ops)["spins"]
        self.assertLess(float(np.max(np.abs(obs["S_bar"] - stat))), 1e-10)

    def test_photon_weights_normalized(self):
        fl = floquet_modes(self.model, drive(1.0, 0.05), self.k, n_step=300,
                           want_wspec=True, m_sp=12)
        self.assertLess(float(np.max(np.abs(
            fl["wspec"]["w_sum"] - 1.0))), 1e-8)

    def test_sambe_vs_evolution(self):
        d = drive(1.0, 0.05)
        res = sambe_vs_evolution(self.model, d, self.k, n_harm=3,
                                 n_step_sambe=400, n_step_evo=300)
        self.assertLess(res["max_mismatch"], 1e-6)

    def test_effective_convergence(self):
        kk = np.array([0.12, 0.04, 0.0])
        errs = {}
        for om in (1.0, 2.0, 4.0):
            d = drive(om, 0.01)
            eps = floquet_modes(self.model, d, kk, n_step=300)["eps"]
            r = effective_split_error(self.model, d, kk, eps, order=2,
                                      n_harm=3, check_validity=False)
            errs[om] = r["max_abs"]
        self.assertLess(errs[4.0], errs[1.0])
        self.assertLess(errs[4.0], 1e-4)

    def test_spin_resolved_spectral_identity(self):
        obs = floquet_spin_observables(self.model, drive(1.0, 0.05), self.k,
                                       self.ops, n_step=300, m_sp=12)
        ident = np.max(np.abs(obs["spin_weights"].sum(axis=0)
                              - obs["S_bar"]))
        self.assertLess(float(ident), 1e-12)
        E = np.linspace(-8, 8, 16001)
        spec = spin_spectral_density(obs, E, 2e-3)
        integ = spec["A_signed"].sum(axis=1) * (E[1] - E[0])
        target = obs["S_bar"].sum(axis=0)
        self.assertLess(float(np.max(np.abs(integ - target))), 2e-4)

    def test_micromotion_is_included(self):
        obs = floquet_spin_observables(self.model, drive(1.0, 0.10), self.k,
                                       self.ops, n_step=300, m_sp=12)
        diff = float(np.max(np.abs(obs["S_bar"] - obs["S_t0"])))
        self.assertTrue(np.isfinite(diff))
        self.assertGreater(diff, 1e-8,
                           "driven spin texture should differ from the t=0 "
                           "approximation (micromotion)")


class TestStageBBenchmarkReport(unittest.TestCase):
    def test_floquet_rashba_benchmark(self):
        from floquet_soc.benchmarks import floquet_rashba
        res = floquet_rashba.run(outdir=None)
        rep = res["report"]
        rep.print()
        self.assertEqual(rep.count(FAIL), 0)
        self.assertEqual(rep.count(BLOCKED), 0)
        self.assertGreaterEqual(rep.count(PASS), 5)


class TestPbNVRegression(unittest.TestCase):
    """The generalized stage-B engine must reproduce the validated PbNV
    engine at the same parameters (n_step=150, pmax=3, circ+ hw=0.4 gamma=0.1).
    """

    @classmethod
    def setUpClass(cls):
        if not os.path.exists(OLD_ENGINE_REF):
            raise unittest.SkipTest(
                "cross-engine reference missing; run gen_old_engine_ref.py")
        cls.cfg = load_config(CFG)
        cls.model = build_model(cls.cfg)
        cls.ref = np.load(OLD_ENGINE_REF)
        cls.path, _ = kpath_from_cfg(cls.cfg, nk_seg=40)
        cls.drive = Drive(hw=0.4, gamma=0.1, a_ref=cls.model.a_ref,
                          polarization="circular", sign=-1)

    def test_quasienergies_match_validated_engine(self):
        dev = 0.0
        for n, ik in enumerate(REF_IDX):
            fl = floquet_modes(self.model, self.drive, self.path[ik],
                               n_step=150, pmax=3)
            dev = max(dev, float(np.max(np.abs(
                np.sort(fl["eps"]) - np.sort(self.ref["eps"][n])))))
        self.assertLess(dev, 1e-8,
                        "stage-B engine deviates from the validated PbNV "
                        "engine (same data, should be machine precision)")


class TestStageBCLI(unittest.TestCase):
    def _small_analytic_cfg(self):
        cfg = json.load(open(os.path.join(HERE, "configs", "rashba.json"),
                             encoding="utf-8"))
        cfg["kmesh"] = {"n": 7, "kmax": 0.10}
        cfg["numerics"] = {"n_step": 80, "m_sp": 8}
        os.makedirs(TMP, exist_ok=True)
        fn = os.path.join(TMP, "rashba_small.json")
        json.dump(cfg, open(fn, "w", encoding="utf-8"), ensure_ascii=False)
        return fn

    def test_spectrum_analytic(self):
        rc = cli_main(["spectrum", "--config", self._small_analytic_cfg(),
                       "--nk", "2", "--kmax", "0.15",
                       "--outdir", os.path.join(TMP, "spec")])
        self.assertEqual(rc, 0)
        self.assertTrue(os.path.exists(os.path.join(TMP, "spec",
                                                    "spectrum.npz")))

    def test_spin_floquet_analytic(self):
        rc = cli_main(["spin_floquet", "--config",
                       self._small_analytic_cfg(),
                       "--outdir", os.path.join(TMP, "spin")])
        self.assertEqual(rc, 0)
        self.assertTrue(os.path.exists(os.path.join(TMP, "spin",
                                                    "spin_floquet.npz")))

    def test_spin_floquet_blocked_without_transform(self):
        import struct
        d = os.path.join(TMP, "fake_spn")
        os.makedirs(d, exist_ok=True)
        fake = os.path.join(d, "wannier90.spn")
        nb, nk = 3, 1
        tri = nb * (nb + 1) // 2
        with open(fake, "wb") as f:
            b = (("Created from WAVECAR at test" + " " * 60)[:60]).encode()
            f.write(struct.pack("<i", len(b)) + b + struct.pack("<i", len(b)))
            dims = struct.pack("<ii", nb, nk)
            f.write(struct.pack("<i", len(dims)) + dims
                    + struct.pack("<i", len(dims)))
            idx = [(i, j) for j in range(nb) for i in range(j + 1)]
            data = np.zeros((tri, 3), dtype=np.complex128)
            for kk, (i, j) in enumerate(idx):
                if i == j:
                    data[kk, 2] = 0.5
            raw = data.reshape(-1).tobytes()
            f.write(struct.pack("<i", len(raw)) + raw
                    + struct.pack("<i", len(raw)))
        cfg = json.load(open(os.path.join(HERE, "configs", "rashba.json"),
                             encoding="utf-8"))
        cfg["spin"] = {"source": "spn", "file": fake}
        fn = os.path.join(TMP, "blocked.json")
        json.dump(cfg, open(fn, "w", encoding="utf-8"), ensure_ascii=False)
        rc = cli_main(["spin_floquet", "--config", fn,
                       "--outdir", os.path.join(TMP, "blocked")])
        self.assertEqual(rc, 2)

    def test_status_lists_stage_b(self):
        self.assertEqual(cli_main(["status"]), 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
