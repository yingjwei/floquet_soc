# -*- coding: utf-8 -*-
"""Process-based k-point parallelism for the spin-texture commands.

Each worker process loads the material model once, computes its k-point chunk
and returns plain numpy arrays; the parent merges the chunks, analyses and
plots.  This scales nearly linearly with ``--jobs`` because the expensive work
(156x156 diagonalisations, 150-step Floquet evolution) is spread over cores.
"""
from __future__ import annotations

import os

import numpy as np


def _chunks(seq, n):
    n = max(1, int(n))
    k = int(np.ceil(len(seq) / n))
    return [seq[i:i + k] for i in range(0, len(seq), k)]


# ---------------------------------------------------------------- worker
def _worker(args):
    """Worker entry point (top-level, picklable).

    args = (job, cfg_path)
      job = ("spin_static", kpts)
      job = ("pst", kpts)
      job = ("spin_floquet", kpts)
    """
    job, cfg_path = args
    kind, kpts = job
    from ..core.model_factory import load_config
    from ..cli import drive_from_cfg, numerics_from_cfg
    from ..core.cache import load_model_spin_cached
    from ..observables.spin_texture import (spin_texture,
                                            spin_texture_wannier_batch)
    from ..observables.floquet_spin import floquet_spin_observables

    cfg = load_config(cfg_path)
    model, spin = load_model_spin_cached(cfg)
    spin_at_k = spin.get("spin_at_k")
    op = spin.get("spin_operator")
    kpts = np.atleast_2d(np.asarray(kpts, float))

    if kind in ("spin_static", "pst"):
        if spin_at_k is not None and op is not None and hasattr(op, "S_R"):
            tex = spin_texture_wannier_batch(model, op, kpts, chunk=64)
            return dict(kpts=kpts, energies=tex["energies"],
                        spins=tex["spins"])
        E = np.empty((len(kpts), model.nw))
        S = np.empty((len(kpts), model.nw, 3))
        for ik, k in enumerate(kpts):
            ops_k = spin_at_k(k) if spin_at_k is not None else spin.get("ops")
            r = spin_texture(model.Hk(k), ops_k)
            E[ik] = r["energies"]
            S[ik] = r["spins"]
        return dict(kpts=kpts, energies=E, spins=S)

    # ---- Floquet spin texture chunk
    drive = drive_from_cfg(cfg, model)
    num = numerics_from_cfg(cfg, model)
    n = len(kpts)
    Sb = np.empty((n, model.nw, 3))
    St = np.empty((n, model.nw, 3))
    Sb_m = np.empty((n, model.nw, 3))
    purity = np.empty((n, model.nw))
    eps = np.empty((n, model.nw))
    st = np.empty((n, model.nw, 3))
    for ik, k in enumerate(kpts):
        ops_k = spin_at_k(k) if spin_at_k is not None else spin.get("ops")
        obs = floquet_spin_observables(model, drive, k, ops_k,
                                       n_step=num["n_step"],
                                       m_sp=num["m_sp"])
        _w, V = np.linalg.eigh(model.Hk(k))
        ov = np.abs(V.conj().T @ obs["evec"]) ** 2
        best = np.argmax(ov, axis=1)
        Sb[ik] = obs["S_bar"]
        St[ik] = obs["S_t0"]
        Sb_m[ik] = obs["S_bar"][best]
        purity[ik] = ov[np.arange(len(best)), best]
        eps[ik] = np.sort(obs["eps"])
        r = spin_texture(model.Hk(k), ops_k)
        st[ik] = r["spins"]
    return dict(kpts=kpts, S_bar=Sb, S_t0=St, S_bar_matched=Sb_m,
                purity=purity, eps=eps, S_static=st)


def run_chunks(kind, cfg_path, kpts, jobs):
    """Split k points over ``jobs`` subprocesses (file-based, no pipes).

    Each worker is a plain ``python -m floquet_soc _kchunk`` subprocess that
    writes its chunk to an .npz file; the parent merges the chunks and removes
    the scratch directory.  Using ``-m`` (rather than a repo-relative
    ``floquet.py``) keeps the parallel path working after ``pip install``.
    """
    import shutil
    import subprocess
    import sys
    import tempfile

    kpts = np.atleast_2d(np.asarray(kpts, float))
    parts = _chunks(list(kpts), jobs)
    tmpdir = tempfile.mkdtemp(prefix="floquet_soc_kchunk_")
    procs = []
    try:
        for i, part in enumerate(parts):
            kf = os.path.join(tmpdir, f"k_{i}.npy")
            np.save(kf, np.asarray(part, float))
            outf = os.path.join(tmpdir, f"out_{i}.npz")
            logf_path = os.path.join(tmpdir, f"log_{i}.txt")
            logf = open(logf_path, "w", encoding="utf-8")
            cmd = [sys.executable, "-m", "floquet_soc", "_kchunk",
                   "--cfg", cfg_path, "--kind", kind,
                   "--kpts", kf, "--out", outf]
            procs.append((subprocess.Popen(cmd, stdout=logf, stderr=logf),
                          outf, logf, logf_path))
        results = []
        errors = []
        for proc, outf, logf, logf_path in procs:
            rc = proc.wait()
            logf.close()
            if rc != 0 or not os.path.exists(outf):
                tail = _tail(logf_path, 20)
                errors.append(f"worker rc={rc} ({logf_path}):\n{tail}")
                continue
            with np.load(outf) as data:
                # copy out immediately: the NpzFile stays open otherwise and
                # locks out_*.npz on Windows, making rmtree below fail
                results.append({k: data[k] for k in data.files})
        if errors:
            raise RuntimeError("parallel k-chunk failure:\n" + "\n".join(errors))
        out = {}
        keys = set().union(*[set(r) for r in results])
        for key in keys:
            out[key] = np.concatenate([r[key] for r in results], axis=0)
        return out
    finally:
        # Windows can briefly hold the scratch dir open right after the last
        # worker exits, so retry once and then drop an already-empty dir.
        shutil.rmtree(tmpdir, ignore_errors=True)
        try:
            os.rmdir(tmpdir)
        except OSError:
            pass


def _tail(path: str, n: int) -> str:
    """Last ``n`` lines of a worker log (for surfacing the real error)."""
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            return "".join(f.readlines()[-n:]).rstrip()
    except OSError:
        return "<log unavailable>"
