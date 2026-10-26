# -*- coding: utf-8 -*-
"""Provenance: version, config/input hashes, missing-file handling."""
import hashlib
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.core.provenance import collect, sha256_file  # noqa: E402
from floquet_soc.version import __version__              # noqa: E402


class TestProvenance(unittest.TestCase):
    def test_collect_hashes_and_missing(self):
        d = tempfile.mkdtemp(prefix="floquet_soc_prov_")
        self.addCleanup(shutil.rmtree, d, True)
        data = os.path.join(d, "hr.dat")
        with open(data, "wb") as f:
            f.write(b"dummy\n")
        cfg = dict(name="fake", kind="wannier90", num_electrons=2,
                   files=dict(hr=data, win=os.path.join(d, "missing.win")),
                   numerics=dict(n_step=150), spin=dict(source="pauli"),
                   drive=dict(hw=0.4, gamma=0.02))
        cpath = os.path.join(d, "cfg.json")
        with open(cpath, "w", encoding="utf-8") as f:
            json.dump(cfg, f)
        prov = collect(cpath, cfg)
        self.assertEqual(prov["version"], __version__)
        self.assertEqual(prov["files"]["hr"]["sha256"],
                         hashlib.sha256(b"dummy\n").hexdigest())
        self.assertFalse(prov["files"]["win"]["exists"])
        self.assertEqual(prov["numerics"]["n_step"], 150)
        self.assertEqual(prov["config"]["sha256"], sha256_file(cpath))


if __name__ == "__main__":
    unittest.main(verbosity=2)
