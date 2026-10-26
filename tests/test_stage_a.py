# -*- coding: utf-8 -*-
"""Stage A tests: input parsing, spin-operator handling, analytic benchmarks.

Run with:  python -m unittest discover -s tests -v
"""
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.cli import load_config, build_model, main as cli_main
from floquet_soc.core.wannier import WannierModel, read_win, read_hr
from floquet_soc.core.spin_operator import (
    read_wannier_spn, validate_spin_matrices, load_spin_for_model,
    MissingTransformError, SPIN_HALF)
from floquet_soc.core.validation import PASS, FAIL, WARN, BLOCKED
from floquet_soc.observables.spin_texture import spin_texture
from floquet_soc.observables.pst import pst_collinearity, axis_angle
from floquet_soc.drives.peierls import Drive, peierls_Hk
from floquet_soc.benchmarks.rashba import RashbaModel
from floquet_soc.benchmarks.rashba_dresselhaus import RashbaDresselhausModel

CFG = os.path.join(HERE, "configs", "pbnv.json")
TMP = os.path.join(HERE, "outputs", "_tests")
PBHV = os.path.join(HERE, "..", "PbNV")


class TestWannierInput(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config(CFG)
        cls.model = build_model(cls.cfg)

    def test_parse_dimensions(self):
        self.assertEqual(self.model.nw, 156)
        self.assertEqual(self.model.nrpts, 49)
        self.assertEqual(len(self.model.win["kpoints"]), 36)
        self.assertEqual(self.model.win["num_bands"], 192)
        self.assertTrue(self.model.win["spinors"])
        self.assertEqual(list(self.model.win["mp_grid"]), [6, 6, 1])

    def test_hermiticity(self):
        err = self.model.hermiticity_errors()
        # the folded R set must be closed under R -> -R, otherwise the H(R)
        # residual is undefined (None) rather than a false 0.0
        self.assertEqual(err["n_unpaired_R"], 0)
        self.assertIsNotNone(err["HR_hermiticity_folded"])
        self.assertLess(err["HR_hermiticity_folded"], 1e-10)
        self.assertLess(err["Hk_hermiticity_expanded"], 1e-10)
        nd = self.model.ndegen_consistency()
        self.assertEqual(nd["ndegen_mismatch"], 0)

    def test_centres(self):
        self.assertEqual(self.model.tau.shape, (156, 3))
        self.assertIn("centres.xyz", self.model.tau_source)

    def test_static_gap_matches_pbnv_reference(self):
        kpts, labels = None, None
        from floquet_soc.cli import kpath_from_cfg
        kpts, labels = kpath_from_cfg(self.cfg)
        e = self.model.static_bands(kpts)
        nocc = self.cfg["num_electrons"]
        gap = float(e[:, nocc].min() - e[:, nocc - 1].max())
        self.assertAlmostEqual(gap, 1.4870, delta=0.01)


class TestSpinOperator(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config(CFG)
        cls.model = build_model(cls.cfg)
        spn_path = cls.cfg["spin"]["file"]
        if not os.path.exists(spn_path):
            raise unittest.SkipTest(f"spn file missing: {spn_path}")
        cls.spn = read_wannier_spn(spn_path)

    def test_spn_structure(self):
        self.assertEqual(self.spn["num_bands"], 192)
        self.assertEqual(self.spn["num_kpts"], 36)
        self.assertEqual(self.spn["matrices"].shape, (36, 3, 192, 192))
        self.assertIn(self.spn["format"], ("formatted-spn", "vasp-binary-spn"))
        self.assertIn("VASP", self.spn["header"])

    def test_spn_numerical_quality_metrics(self):
        metrics, checks = validate_spin_matrices(self.spn["matrices"])
        self.assertLess(metrics["hermiticity"], 1e-12)
        self.assertLess(metrics["diag_imag"], 1e-12)
        self.assertLessEqual(metrics["diag_pol_max"], 1.0 + 1e-3)
        # New (VASP LWRITE_SPN) data: Pauli bound satisfied, sum rule close to 3
        self.assertLessEqual(metrics["eig_max"], 1.0 + 2e-3)
        self.assertGreaterEqual(metrics["eig_min"], -1.0 - 2e-3)
        self.assertLessEqual(metrics["S2_max"], 3.0 + 2e-2)

    def test_missing_transform_error(self):
        """In a directory with only a spn file, rotation must be refused."""
        import struct

        d = os.path.join(TMP, "spn_only")
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
        cfg = {"spin": {"source": "spn", "file": fake}}

        class FakeModel:
            nw = 2

        with self.assertRaises(MissingTransformError) as ctx:
            load_spin_for_model(cfg, FakeModel())
        self.assertIn("wannier90.chk", str(ctx.exception))
        self.assertTrue(any("wannier90.chk" in p
                            for p in ctx.exception.missing_files))


class TestFormattedSpn(unittest.TestCase):
    def test_formatted_reader(self):
        from floquet_soc.core.spin_operator import read_formatted_spn
        d = os.path.join(TMP, "fmt_spn")
        os.makedirs(d, exist_ok=True)
        nb, nk = 3, 2
        path = os.path.join(d, "wannier90.spn")
        lines = ["Written by VASP: test system", f"{nb} {nk}"]
        expected = {}
        for ik in range(nk):
            for m in range(nb):
                for n in range(m + 1):
                    for c in range(3):
                        z = complex(0.1 * (ik + 1) + 0.01 * n,
                                    0.2 * m + 0.001 * c)
                        expected[(ik, c, n, m)] = z
                        lines.append(f"{z.real:.16E} {z.imag:.16E}")
        with open(path, "w", encoding="ascii") as f:
            f.write("\n".join(lines) + "\n")
        spn = read_formatted_spn(path)
        self.assertEqual(spn["num_bands"], nb)
        self.assertEqual(spn["num_kpts"], nk)
        self.assertEqual(spn["format"], "formatted-spn")
        M = spn["matrices"]
        for (ik, c, n, m), z in expected.items():
            self.assertLess(abs(M[ik, c, n, m] - z), 1e-12)
        for ik in range(nk):
            for c in range(3):
                for m in range(nb):
                    for n in range(m + 1):
                        if n != m:
                            self.assertLess(
                                abs(M[ik, c, n, m] - np.conj(M[ik, c, m, n])),
                                1e-12)


class TestAnalyticBenchmarks(unittest.TestCase):
    def test_rashba_analytic(self):
        m = RashbaModel()
        for k in np.linspace(0.02, 0.3, 21):
            self.assertAlmostEqual(m.splitting([k, 0, 0]), 2 * m.alpha * k,
                                   places=12)
            for ang in np.linspace(0, 2 * np.pi, 7, endpoint=False):
                kk = np.array([k * np.cos(ang), k * np.sin(ang), 0.0])
                res = m.spin_at(kk)
                s = res["spins"][0]
                s_an = 0.5 * np.array([kk[1], -kk[0], 0.0]) / k
                sign = 1.0 if np.dot(s, s_an) >= 0 else -1.0
                self.assertLess(np.linalg.norm(s - sign * s_an), 1e-10)
                self.assertLess(abs(s[2]), 1e-12)

    def test_rashba_benchmark_report(self):
        from floquet_soc.benchmarks import rashba
        res = rashba.run(outdir=None)
        rep = res["report"]
        self.assertEqual(rep.count(FAIL), 0)
        self.assertEqual(rep.count(BLOCKED), 0)
        self.assertGreaterEqual(rep.count(PASS), 4)

    def test_rashba_dresselhaus_pst(self):
        m = RashbaDresselhausModel(alpha=0.1, beta=0.1)
        axis = np.array([1.0, -1.0, 0.0]) / np.sqrt(2)
        circ = np.array([[0.2 * np.cos(a), 0.2 * np.sin(a), 0.0]
                         for a in np.linspace(0, 2 * np.pi, 25, endpoint=False)])
        spins = np.array([m.spin_at(k)["spins"][0] for k in circ])
        pst = pst_collinearity(spins)
        self.assertAlmostEqual(pst["collinearity"], 1.0, places=8)
        self.assertLess(axis_angle(pst["pst_axis"], axis), 1e-8)

    def test_rashba_dresselhaus_benchmark_report(self):
        from floquet_soc.benchmarks import rashba_dresselhaus
        res = rashba_dresselhaus.run(outdir=None)
        rep = res["report"]
        self.assertEqual(rep.count(FAIL), 0)
        self.assertGreaterEqual(rep.count(PASS), 4)


class TestDriveAndPeierls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = load_config(CFG)
        cls.model = build_model(cls.cfg)

    def test_drive_parameters(self):
        d = Drive(hw=0.4, gamma=0.1, a_ref=self.model.a_ref,
                  polarization="circular", sign=+1)
        self.assertAlmostEqual(d.period, 2 * np.pi / 0.4)
        A = d.A(np.array([0.0, d.period / 4]))
        self.assertEqual(A.shape, (2, 3))
        # circular: |A| constant
        self.assertAlmostEqual(np.linalg.norm(A[0]), np.linalg.norm(A[1]),
                               places=12)

    def test_peierls_reduces_to_hk_at_zero_field(self):
        k = np.array([0.13, 0.27, 0.0])
        H0 = self.model.Hk(k)
        Hp = peierls_Hk(self.model, k, np.zeros(3))
        self.assertLess(np.max(np.abs(H0 - Hp)), 1e-12)

    def test_peierls_hermiticity(self):
        k = np.array([0.1, 0.2, 0.0])
        A = np.array([0.01, -0.02, 0.0])
        H = peierls_Hk(self.model, k, A)
        self.assertLess(np.max(np.abs(H - H.conj().T)), 1e-10)


class TestCLI(unittest.TestCase):
    def test_validate_pbnv(self):
        rc = cli_main(["validate", "--config", CFG, "--outdir", TMP])
        self.assertEqual(rc, 0,
                         "chk rotation should make the PbNV spin operator usable")
        md = os.path.join(TMP, "validate.md")
        self.assertTrue(os.path.exists(md))
        with open(md, encoding="utf-8") as f:
            text = f.read()
        self.assertIn("spin_operator_wannier_basis", text)
        self.assertIn("PASS", text)
        self.assertIn("wannier90.spn", text)

    def test_benchmark_all(self):
        rc = cli_main(["benchmark", "--which", "all",
                       "--outdir", os.path.join(TMP, "bench")])
        self.assertEqual(rc, 0)

    def test_bands_and_status(self):
        self.assertEqual(cli_main(["bands", "--config", CFG,
                                   "--outdir", os.path.join(TMP, "bands")]), 0)
        self.assertTrue(os.path.exists(
            os.path.join(TMP, "bands", "bands.npz")))
        self.assertEqual(cli_main(["status"]), 0)

    def test_unknown_command(self):
        with self.assertRaises(SystemExit) as ctx:
            cli_main(["not-a-command"])
        self.assertEqual(ctx.exception.code, 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
