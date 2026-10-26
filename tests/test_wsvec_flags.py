# -*- coding: utf-8 -*-
"""Batch D: wannier90 wsvec.dat flag semantics (no silent folded fallback)."""
from __future__ import annotations

import os
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.core.wannier import (WannierModel, WannierInputError,   # noqa: E402
                                     header_flag, read_wsvec)

BODY_ONE_BLOCK = "   0   0   0   1   2\n1\n   1   0   0\n"


def _write(tmp, header):
    path = os.path.join(tmp, "wannier90_wsvec.dat")
    with open(path, "w", encoding="utf-8") as f:
        f.write(header)
        f.write(BODY_ONE_BLOCK)
    return path


def _model(wsvec_path=None, override=None):
    win = dict(num_wann=2, num_bands=2, spinors=True,
               lattice=np.eye(3) * 4.0, mp_grid=np.array([2, 2, 1]),
               kpoints=np.zeros((1, 3)))
    h0 = np.array([[0.0, 0.25], [0.25, 0.3]], complex)
    h1 = np.array([[0.1, 0.2j], [-0.2j, -0.1]], complex)
    hr = dict(num_wann=2, nrpts=3, R=np.array([(0, 0, 0), (1, 0, 0),
                                               (-1, 0, 0)]),
              H=np.asarray([h0, h1, h1.conj().T], complex),
              ndegen=np.ones(3, int))
    return WannierModel(win, hr, tau=np.zeros((2, 3)), wsvec_path=wsvec_path,
                        wsvec_use_ws_distance=override)


class TestHeaderFlag(unittest.TestCase):

    def test_parses_true_false_and_absent(self):
        self.assertTrue(header_flag(" use_ws_distance=.true.", "use_ws_distance"))
        self.assertTrue(header_flag(" USE_WS_DISTANCE = .TRUE. ", "use_ws_distance"))
        self.assertFalse(header_flag(" use_ws_distance=.false.", "use_ws_distance"))
        self.assertFalse(header_flag("use_ws_distance = .false.", "use_ws_distance"))
        self.assertIsNone(header_flag(" just a comment line", "use_ws_distance"))
        self.assertIsNone(header_flag("", "use_ws_distance"))

    def test_absent_token_is_none_not_false(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, " wannier90 wsvec\n")
            wv = read_wsvec(path)
            self.assertIsNone(wv["use_ws_distance"])
            self.assertIsNone(wv["header_use_ws_distance"])
            self.assertEqual(wv["blocks"], 1)


class TestModelSemantics(unittest.TestCase):

    def test_token_free_header_infers_true_and_records_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _model(_write(tmp, " wannier90 wsvec\n"))
            self.assertIsNotNone(m.HR_full, "expansion must be applied")
            self.assertIs(m.wsvec_info["use_ws_distance"], True)
            self.assertEqual(m.wsvec_info["use_ws_distance_source"],
                             "inferred (wsvec.dat exists)")
            self.assertIsNone(m.wsvec_info["header_use_ws_distance"])

    def test_header_true_is_honoured(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _model(_write(tmp, " use_ws_distance=.true.\n"))
            self.assertIsNotNone(m.HR_full)
            self.assertEqual(m.wsvec_info["use_ws_distance_source"], "header")

    def test_explicit_false_with_shift_blocks_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write(tmp, " use_ws_distance=.false.\n")
            with self.assertRaises(WannierInputError):
                _model(path)

    def test_config_override_wins_and_disables_expansion(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _model(_write(tmp, " use_ws_distance=.true.\n"), override=False)
            self.assertIsNone(m.HR_full, "override must disable the expansion")
            self.assertFalse(m.wsvec_info["applied_to_model"])
            self.assertEqual(m.wsvec_info["use_ws_distance_source"],
                             "config override")


if __name__ == "__main__":
    unittest.main()
