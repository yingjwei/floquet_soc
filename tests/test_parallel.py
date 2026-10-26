# -*- coding: utf-8 -*-
"""Parallel subprocess path (python -m floquet_soc _kchunk) round-trip.

This is the only test that exercises ``run_chunks`` without needing real
Wannier90 data: an analytic 2-band model uses ``pauli`` spin operators, so the
worker runs end to end.  It is the regression net for the Windows-only
NpzFile-handle bug (open ``np.load`` results locked ``out_*.npz`` and broke the
temp-dir cleanup) and for the log-tail surfaced on worker failure.
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.workflows.parallel import run_chunks          # noqa: E402


def _rashba_cfg(path):
    cfg = {"name": "Rashba", "kind": "analytic",
           "model": {"type": "rashba", "hbar2_over_2m": 0.5, "alpha": 0.1},
           "spin": {"source": "pauli"},
           "drive": {"hw": 1.0, "gamma": 0.05, "polarization": "circular",
                     "a_ref": 1.0}}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(cfg, f)
    return path


class TestParallelRunChunks(unittest.TestCase):

    def test_spin_static_roundtrip_two_workers(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = _rashba_cfg(os.path.join(d, "cfg.json"))
            kpts = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.2, 0.0, 0.0]]
            out = run_chunks("spin_static", cfg, kpts, jobs=2)
            self.assertEqual(out["spins"].shape, (3, 2, 3))
            self.assertEqual(out["energies"].shape, (3, 2))
            self.assertTrue(np.all(np.isfinite(out["spins"])))

    def test_worker_failure_surfaces_log_tail(self):
        with tempfile.TemporaryDirectory() as d:
            # no "files" and no analytic model -> build_model raises in worker
            cfg = os.path.join(d, "bad.json")
            with open(cfg, "w", encoding="utf-8") as f:
                json.dump({"name": "x", "kind": "wannier90",
                           "files": {"win": os.path.join(d, "missing.win")}}, f)
            with self.assertRaises(RuntimeError) as cm:
                run_chunks("spin_static", cfg,
                           [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]], jobs=2)
            self.assertIn("worker rc=", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
