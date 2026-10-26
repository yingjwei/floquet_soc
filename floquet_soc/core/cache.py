# -*- coding: utf-8 -*-
"""Parsed-model cache (saves ~20 s of parsing per command / worker).

Design constraints learned the hard way:

* the cache must **not** live inside the package tree (``outputs/_cache`` in the
  repo): any run recreated it after a manual ``rm``, so "delete the cache" was
  never a durable action.  Location: ``$FLOQUET_SOC_CACHE`` or
  ``<system temp>/floquet_soc-cache``.
* the key must carry a **schema version**, otherwise a cache written by an
  older layout is unpickled into the current classes (``__setstate__`` fills
  the missing fields silently) and quietly feeds stale objects.
* writes must be **atomic** and the store **bounded**: worker processes hit the
  same key concurrently, and each entry of a 156-orbital model is ~180 MB.
"""
import hashlib
import json
import os
import pickle
import tempfile

# bump when the pickled payload layout changes (new/renamed model attributes,
# different spin-operator container, ...)
CACHE_SCHEMA = 3
DEFAULT_MAX_ENTRIES = 8
# content-hash files up to this size (key inputs are hr.dat/chk/eig/spn, which
# are small; the ~180 MB objects are the *value*, never a key input)
_CONTENT_HASH_MAX = 64 * 1024 * 1024


def _content_hash(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def cache_dir() -> str:
    """Cache root: explicit env var wins, otherwise the system temp dir."""
    env = os.environ.get("FLOQUET_SOC_CACHE")
    if env:
        return os.path.abspath(env)
    return os.path.join(tempfile.gettempdir(), "floquet_soc-cache")


def _max_entries() -> int:
    try:
        return max(1, int(os.environ.get("FLOQUET_SOC_CACHE_MAX",
                                         DEFAULT_MAX_ENTRIES)))
    except (TypeError, ValueError):
        return DEFAULT_MAX_ENTRIES


def _cache_key(cfg):
    from ..version import __version__
    h = hashlib.sha1()
    h.update(f"floquet_soc:{__version__}:schema{CACHE_SCHEMA}".encode())
    h.update(os.path.abspath(cfg.get("_config_dir", "")).encode())
    for name, p in sorted((cfg.get("files") or {}).items()):
        if p is None:
            continue
        if isinstance(p, (str, os.PathLike)):
            if not os.path.exists(p):
                continue
            st = os.stat(p)
            h.update(f"{name}:{st.st_size}:{st.st_mtime_ns}".encode())
            # content hash for small inputs: size+mtime is blind to cp -p /
            # rsync -a / git-preserve-timestamps copies of an *edited* file
            if st.st_size <= _CONTENT_HASH_MAX:
                h.update(f":sha256:{_content_hash(p)}".encode())
        else:
            # config-only switches (e.g. files.wsvec_use_ws_distance=True/False)
            # directly change how the model is built; they are part of the key.
            # NB: booleans must NOT go through os.path.exists above (TypeError).
            h.update(f"{name}:{type(p).__name__}:{p}".encode())
    # v0.1.5-rc1: analytic-model configs have no data files, so the file-stat
    # hash alone gave the SAME key to every analytic config in one directory
    # (e.g. rashba.json and rashba_dresselhaus.json) and a stale model could
    # be loaded from the cache.  Include kind/name/model parameters.
    h.update(f"kind:{cfg.get('kind')}".encode())
    h.update(f"name:{cfg.get('name')}".encode())
    h.update(json.dumps(cfg.get("model"), sort_keys=True,
                        default=str).encode())
    spin = cfg.get("spin") or {}
    # spin.source ("pauli" vs "spn") selects entirely different code paths;
    # switching it must invalidate the cache even if no file changed.
    h.update(f"spin_source:{spin.get('source')}".encode())
    for key in ("file", "chk", "eig"):
        p = spin.get(key)
        if p and os.path.exists(p):
            st = os.stat(p)
            h.update(f"spin:{key}:{st.st_size}:{st.st_mtime_ns}".encode())
            if st.st_size <= _CONTENT_HASH_MAX:
                h.update(f":sha256:{_content_hash(p)}".encode())
    return h.hexdigest()[:24]


def _evict(directory: str, keep: int) -> None:
    """Drop the oldest entries beyond ``keep`` (best effort, never raises)."""
    try:
        entries = [os.path.join(directory, f) for f in os.listdir(directory)
                   if f.endswith(".pkl")]
        entries.sort(key=lambda p: os.path.getmtime(p))
        for path in entries[:-keep] if keep < len(entries) else []:
            try:
                os.remove(path)
            except OSError:
                pass
    except OSError:
        pass


def _atomic_write(path: str, obj) -> None:
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "wb") as f:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def load_model_spin_cached(cfg):
    """Return (model, spin) using an on-disk pickle cache when possible."""
    from .model_factory import build_model
    from .spin_operator import load_spin_for_model

    directory = cache_dir()
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, _cache_key(cfg) + ".pkl")
    if os.path.exists(path):
        try:
            with open(path, "rb") as f:
                obj = pickle.load(f)
            return obj["model"], obj["spin"]
        except Exception:                                 # noqa: BLE001
            # stale/corrupt entry (half-written by an older version, or from
            # another schema): rebuild instead of failing the command
            try:
                os.remove(path)
            except OSError:
                pass
    model = build_model(cfg)
    spin = load_spin_for_model(cfg, model)
    try:
        _atomic_write(path, dict(model=model, spin=spin,
                                 schema=CACHE_SCHEMA))
        _evict(directory, _max_entries())
    except OSError:
        pass
    return model, spin
