# -*- coding: utf-8 -*-
"""Run provenance: version, config hash, input-file hashes, numerical settings.

The framework must record enough to reproduce a run: version, config path and
SHA256, every existing input file with size/hash, precision/solver settings
and the spin-operator source.  This module is material agnostic.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from datetime import datetime, timezone

import numpy as np


def sha256_file(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def collect(config_path: str, cfg: dict | None = None,
            extra: dict | None = None) -> dict:
    from ..version import __version__

    cfg = cfg if cfg is not None else json.load(
        open(config_path, encoding="utf-8"))
    files = {}
    for key, path in (cfg.get("files") or {}).items():
        if not path:
            continue
        abspath = os.path.abspath(path)
        entry = dict(path=abspath, exists=os.path.exists(abspath))
        if entry["exists"]:
            entry.update(size=int(os.path.getsize(abspath)),
                         sha256=sha256_file(abspath))
        files[key] = entry
    prov = dict(
        framework="floquet_soc",
        version=__version__,
        created_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        python=sys.version.split()[0],
        numpy=np.__version__,
        platform=platform.platform(),
        config=dict(path=os.path.abspath(config_path),
                    sha256=sha256_file(os.path.abspath(config_path))),
        material=cfg.get("name"),
        kind=cfg.get("kind"),
        num_electrons=cfg.get("num_electrons"),
        files=files,
        spin=dict(cfg.get("spin") or {}),
        numerics=dict(cfg.get("numerics") or {}),
        drive=dict(cfg.get("drive") or {}),
    )
    if extra:
        prov["extra"] = extra
    return prov


def write(prov: dict, outdir: str, name: str = "provenance.json") -> str:
    os.makedirs(outdir, exist_ok=True)
    path = os.path.join(outdir, name)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=2, ensure_ascii=False)
    return path
