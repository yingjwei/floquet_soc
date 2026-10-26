# -*- coding: utf-8 -*-
"""Batch G: model cache location, schema version, atomic write, eviction.

Before this batch the cache lived in ``<repo>/outputs/_cache`` (any run
recreated it after a manual delete), the key had no schema version (an old
pickle was silently unpickled into the current classes), writes were not
atomic, and nothing was ever evicted (two identical 182 MB entries were found
on disk).
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


class TestCacheLayout(unittest.TestCase):

    def test_dir_is_outside_the_package(self):
        d = os.path.abspath(cache_mod.cache_dir())
        self.assertFalse(d.startswith(os.path.abspath(HERE)),
                         f"cache must not live inside the source tree: {d}")

    def test_env_override_wins(self):
        with tempfile.TemporaryDirectory() as tmp:
            old = os.environ.get("FLOQUET_SOC_CACHE")
            os.environ["FLOQUET_SOC_CACHE"] = tmp
            try:
                self.assertEqual(os.path.abspath(cache_mod.cache_dir()),
                                 os.path.abspath(tmp))
            finally:
                if old is None:
                    os.environ.pop("FLOQUET_SOC_CACHE", None)
                else:
                    os.environ["FLOQUET_SOC_CACHE"] = old

    def test_key_carries_the_schema_version(self):
        cfg = {"kind": "analytic", "name": "x", "model": {"alpha": 0.1}}
        before = cache_mod._cache_key(cfg)
        old = cache_mod.CACHE_SCHEMA
        try:
            cache_mod.CACHE_SCHEMA = old + 1
            self.assertNotEqual(cache_mod._cache_key(cfg), before,
                                "a schema bump must produce a different key")
        finally:
            cache_mod.CACHE_SCHEMA = old


class TestAtomicWriteAndEviction(unittest.TestCase):

    def test_no_temp_files_left_behind(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "entry.pkl")
            cache_mod._atomic_write(path, {"a": 1})
            self.assertTrue(os.path.exists(path))
            self.assertEqual([f for f in os.listdir(tmp) if ".tmp" in f], [])

    def test_eviction_keeps_newest(self):
        with tempfile.TemporaryDirectory() as tmp:
            for i in range(4):
                p = os.path.join(tmp, f"e{i}.pkl")
                with open(p, "wb") as f:
                    f.write(b"x")
                os.utime(p, (1000 + i, 1000 + i))
            cache_mod._evict(tmp, keep=2)
            left = sorted(os.listdir(tmp))
            self.assertEqual(left, ["e2.pkl", "e3.pkl"],
                             "oldest entries must be dropped first")

    def test_eviction_survives_missing_dir(self):
        cache_mod._evict(os.path.join(tempfile.gettempdir(),
                                      "definitely-not-here-xyz"), 2)


if __name__ == "__main__":
    unittest.main()
