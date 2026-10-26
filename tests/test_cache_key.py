# -*- coding: utf-8 -*-
"""A2: the model-cache key must change when the *build inputs* change.

Before this batch the key missed ``spin.source`` (two completely different
code paths) and any boolean switch in ``files`` (e.g.
``wsvec_use_ws_distance``), and it fed booleans straight into
``os.path.exists`` (a ``TypeError``).  It also relied on size+mtime only, so a
content edit with a preserved timestamp went unnoticed.
"""
from __future__ import annotations

import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.core import cache as cache_mod      # noqa: E402


def _cfg(**overrides):
    cfg = {"kind": "analytic", "name": "x", "model": {"alpha": 0.1},
           "_config_dir": os.path.dirname(__file__)}
    cfg.update(overrides)
    return cfg


class TestCacheKeyCompleteness(unittest.TestCase):

    def test_spin_source_is_hashed(self):
        a = _cfg(spin={"source": "pauli"})
        b = _cfg(spin={"source": "spn"})
        self.assertNotEqual(cache_mod._cache_key(a), cache_mod._cache_key(b),
                            "switching spin.source must invalidate the cache")

    def test_boolean_switch_is_hashed_and_does_not_crash(self):
        t = _cfg(files={"wsvec_use_ws_distance": True})
        f = _cfg(files={"wsvec_use_ws_distance": False})
        self.assertNotEqual(cache_mod._cache_key(t), cache_mod._cache_key(f))
        # regression: the boolean used to flow into os.path.exists(True) ->
        # TypeError
        cache_mod._cache_key(t)

    def test_file_content_is_hashed_not_just_stat(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "hr.dat")
            with open(p, "w", encoding="utf-8") as fh:
                fh.write("1.0 0.0\n" * 16)
            st = os.stat(p)
            a = _cfg(files={"hr": p})
            key_a = cache_mod._cache_key(a)
            # same size, restore the same mtime, different content
            with open(p, "w", encoding="utf-8") as fh:
                fh.write("2.0 0.0\n" * 16)
            os.utime(p, (st.st_atime, st.st_mtime))
            key_b = cache_mod._cache_key(a)
            self.assertNotEqual(key_a, key_b,
                                "an edited file with preserved mtime must "
                                "produce a different key")


if __name__ == "__main__":
    unittest.main()
