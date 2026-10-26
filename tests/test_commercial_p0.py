# -*- coding: utf-8 -*-
"""Commercial-hardening P0 regressions: read_hr truncation, effective order,
config validation."""
from __future__ import annotations

import os
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.core.wannier import read_hr, WannierInputError        # noqa: E402
from floquet_soc.core.model_factory import validate_cfg, load_config   # noqa: E402
from floquet_soc.solvers.effective import effective_hamiltonian        # noqa: E402

_HR_DIALECT_A = "\n".join([
    "2",            # nw
    "1",            # nrpts
    "1",            # ndegen
    "   0   0   0",                       # R block header
    "   1   1   0.5   0.0",
    "   1   2   0.0   0.1",
    "   2   1   0.0  -0.1",
    "   2   2   0.5   0.0",
]) + "\n"

_HR_DIALECT_B = "\n".join([
    "2", "1", "1",
    "   0   0   0   1   1   0.5   0.0",
    "   0   0   0   1   2   0.0   0.1",
    "   0   0   0   2   1   0.0  -0.1",
    "   0   0   0   2   2   0.5   0.0",
]) + "\n"


class TestReadHrTruncation(unittest.TestCase):

    def _parse(self, text):
        with tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False,
                                         encoding="utf-8") as f:
            f.write(text)
            path = f.name
        try:
            return read_hr(path)
        finally:
            os.remove(path)

    def test_dialect_a_truncated_block_raises(self):
        truncated = _HR_DIALECT_A.rsplit("   2   2   0.5   0.0\n", 1)[0] + "\n"
        with self.assertRaises(WannierInputError):
            self._parse(truncated)

    def test_dialect_b_truncated_block_raises(self):
        truncated = _HR_DIALECT_B.rsplit("   0   0   0   2   2   0.5   0.0\n",
                                         1)[0] + "\n"
        with self.assertRaises(WannierInputError):
            self._parse(truncated)

    def test_complete_files_parse(self):
        for text in (_HR_DIALECT_A, _HR_DIALECT_B):
            hr = self._parse(text)
            self.assertEqual(hr["num_wann"], 2)
            self.assertEqual(hr["nrpts"], 1)
            self.assertEqual(hr["H"].shape, (1, 2, 2))


class TestEffectiveOrder(unittest.TestCase):

    def test_order_above_two_is_rejected(self):
        with self.assertRaises(ValueError):
            effective_hamiltonian(None, None, None, order=5)


class TestConfigValidation(unittest.TestCase):

    def test_unknown_key_raises(self):
        with self.assertRaises(ValueError):
            validate_cfg({"kind": "analytic", "name": "x", "bogus_key": 1})

    def test_negative_gamma_raises(self):
        with self.assertRaises(ValueError):
            validate_cfg({"drive": {"gamma": -1.0}})

    def test_low_m_sp_raises(self):
        with self.assertRaises(ValueError):
            validate_cfg({"numerics": {"m_sp": 1}})

    def test_zero_jobs_raises(self):
        with self.assertRaises(ValueError):
            validate_cfg({"jobs": 0})

    def test_valid_config_passes(self):
        validate_cfg({"kind": "analytic", "name": "x",
                      "drive": {"gamma": 0.1},
                      "numerics": {"m_sp": 8, "n_step": 150}, "jobs": 2})

    def test_shipped_configs_pass(self):
        for name in ("rashba.json", "rashba_dresselhaus.json",
                     "pt.json", "pbnv.json"):
            load_config(os.path.join(HERE, "configs", name))

    def test_unknown_subkey_raises(self):
        with self.assertRaises(ValueError):
            validate_cfg({"kind": "analytic", "name": "x",
                          "spin": {"sorce": "pauli"}})   # typo: sorce
        with self.assertRaises(ValueError):
            validate_cfg({"kind": "analytic", "name": "x",
                          "numerics": {"m_sp": 8, "nstep": 150}})  # typo

    def test_known_subkeys_pass(self):
        validate_cfg({"kind": "analytic", "name": "x",
                      "spin": {"source": "pauli"},
                      "drive": {"hw": 1.0, "gamma": 0.05,
                                "polarization": "circular", "sign": -1},
                      "numerics": {"n_step": 150, "m_sp": 16, "pmax": 3}})


if __name__ == "__main__":
    unittest.main()
