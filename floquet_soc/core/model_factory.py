# -*- coding: utf-8 -*-
"""Model construction from a material config (moved out of ``cli``).

``load_config`` and ``build_model`` were the reason ``core.cache`` and
``workflows.parallel`` had to reach into the CLI (``from ..cli import ...``).
Keeping them in ``core`` breaks that dependency: the numerical core no longer
imports the command-line layer.
"""
from __future__ import annotations

import json
import os

from .wannier import WannierModel, WannierInputError

# Top-level keys the framework understands.  Everything else is refused: a
# misspelled key is silently ignored by every ``cfg.get``, which is the single
# most common source of "the config looks right but does nothing".
_KNOWN_TOP_LEVEL_KEYS = {
    "name", "kind", "model", "root", "files", "spin", "drive", "numerics",
    "kmesh", "kmesh_pst", "kmesh_spin_floquet", "kplane", "kpath", "nk_seg",
    "num_electrons", "num_wann", "band_window", "system", "jobs",
    "soc", "pst", "rashba_fit", "scan", "converge", "spectrum", "notes",
    "floquet_match_tol", "spin_inplane_ratio", "spin_plot_band",
    "_root", "_config_dir",
}

# Sub-dict keys the framework reads.  A misspelled sub-key is silently ignored
# by every ``block.get("k")`` and wastes hours; refuse it like top-level keys.
# These sets are the authoritative allowlists -- extend them here when a new
# option is added, then document it in docs/CONFIG_REFERENCE.md.
_KNOWN_SUBKEYS = {
    "files": {
        "win", "hr", "centres", "chk", "eig", "spn", "amn", "mmn", "wsvec",
        "wout", "bands_ref", "bands_ref_kpt", "wsvec_use_ws_distance",
    },
    "spin": {"source", "file", "chk", "eig", "u_mat", "u_dis_mat", "notes"},
    "drive": {"hw", "gamma", "polarization", "sign", "a_ref", "ratio", "phi",
              "e1", "e2", "ellipticity", "amplitude"},
    "numerics": {"n_step", "m_sp", "pmax"},
    "kmesh": {"n", "kmax", "dk", "kind"},
    "kmesh_pst": {"n", "kmax", "dk", "kind"},
    "kmesh_spin_floquet": {"n", "kmax", "dk", "kind"},
    "scan": {"hw", "gamma", "kpt", "n_step", "m_sp", "window"},
    "converge": {"kpt", "n_step", "m_sp", "window", "tol_ebar", "tol_spin",
                 "tol_wsum"},
    "rashba_fit": {"anchor", "direction", "nk", "kmax", "kmin", "step",
                   "r2_thr", "alpha_ratio_thr", "split_r2_thr"},
    "system": {"type", "inversion", "time_reversal"},
    "soc": {"emin", "emax", "max_gap"},
    "pst": {"bands", "min_polarization", "min_splitting"},
    "spectrum": {"nE", "emin", "emax", "sigma"},
    "model": {"type", "hbar2_over_2m", "alpha", "beta"},
}


def validate_cfg(cfg: dict) -> None:
    """Fail fast on unknown top-level/sub keys and out-of-range numerics."""
    unknown = sorted(set(cfg) - _KNOWN_TOP_LEVEL_KEYS)
    if unknown:
        raise ValueError(
            f"unknown config key(s): {', '.join(unknown)}; known top-level "
            f"keys: {', '.join(sorted(_KNOWN_TOP_LEVEL_KEYS))}")

    for block, allowed in _KNOWN_SUBKEYS.items():
        b = cfg.get(block)
        if isinstance(b, dict):
            bad = sorted(set(b) - allowed)
            if bad:
                raise ValueError(
                    f"unknown key(s) in config[{block!r}]: {', '.join(bad)}; "
                    f"known: {', '.join(sorted(allowed))}")

    drive = cfg.get("drive") or {}
    gamma = drive.get("gamma")
    if gamma is not None:
        try:
            gamma = float(gamma)
        except (TypeError, ValueError):
            raise ValueError(f"drive.gamma is not a number: {gamma!r}")
        if gamma < 0:
            raise ValueError(f"drive.gamma must be >= 0, got {gamma}")

    numerics = cfg.get("numerics") or {}
    for key, lo in (("m_sp", 2), ("n_step", 1)):
        v = numerics.get(key)
        if v is not None and (not isinstance(v, (int, float)) or v < lo):
            raise ValueError(f"numerics.{key} must be >= {lo}, got {v!r}")

    jobs = cfg.get("jobs")
    if jobs is not None and (not isinstance(jobs, (int, float)) or jobs < 1):
        raise ValueError(f"jobs must be >= 1, got {jobs!r}")


def load_config(path: str) -> dict:
    # utf-8-sig: Windows editors (Notepad, PowerShell 5.1, Excel) happily add a
    # BOM, which plain utf-8 rejects with a JSONDecodeError on column 1.
    with open(path, encoding="utf-8-sig") as f:
        cfg = json.load(f)
    root = cfg.get("root")
    if root:
        # ${VAR} / ~ let the same config run on a different machine or OS
        root = os.path.expandvars(os.path.expanduser(str(root)))
    base = cfg.get("_config_dir", os.path.dirname(os.path.abspath(path)))
    if root:
        base = root if os.path.isabs(root) else os.path.join(base, root)
    cfg["_root"] = os.path.normpath(base)
    cfg["_config_dir"] = os.path.dirname(os.path.abspath(path))
    for key in ("files",):
        for name, p in (cfg.get(key) or {}).items():
            if p is None or not isinstance(p, (str, os.PathLike)):
                # booleans/numbers are switches (e.g. wsvec_use_ws_distance),
                # not paths -- never feed them to os.path.isabs
                continue
            p = os.path.expandvars(os.path.expanduser(str(p)))
            if not os.path.isabs(p):
                p = os.path.normpath(os.path.join(base, p))
            cfg[key][name] = p
    spin = cfg.get("spin", {})
    for key in ("file", "chk", "eig", "u_mat", "u_dis_mat"):
        if spin.get(key):
            p = os.path.expandvars(os.path.expanduser(str(spin[key])))
            if not os.path.isabs(p):
                p = os.path.normpath(os.path.join(base, p))
            spin[key] = p
    validate_cfg(cfg)
    return cfg


def build_model(cfg):
    kind = cfg.get("kind")
    if kind == "wannier90":
        files = cfg["files"]
        win_path = files.get("win")
        hr_path = files.get("hr")
        chk_path = files.get("chk")
        eig_path = files.get("eig")
        has_win = bool(win_path and os.path.exists(win_path))
        has_hr = bool(hr_path and os.path.exists(hr_path))
        has_chk = bool(chk_path and os.path.exists(chk_path))
        has_eig = bool(eig_path and os.path.exists(eig_path))
        wsvec_override = files.get("wsvec_use_ws_distance")
        if has_win and has_hr:
            return WannierModel.from_files(
                win_path, hr_path, files.get("centres"),
                wsvec_path=files.get("wsvec"),
                wsvec_use_ws_distance=wsvec_override)
        if has_chk and has_hr:
            # win file not available: reconstruct metadata from the checkpoint
            return WannierModel.from_chk(
                chk_path, hr_path, files.get("centres"),
                wsvec_path=files.get("wsvec"),
                wsvec_use_ws_distance=wsvec_override)
        if has_chk and has_eig:
            # no hr.dat: rebuild H(R) from the checkpoint on its k-mesh
            return WannierModel.from_chk_eig(chk_path, eig_path,
                                             win_path if has_win else None)
        raise WannierInputError(
            "Wannier model needs one of: win+hr, chk+hr, or chk+eig "
            "(win/hr missing)")
    if kind == "analytic":
        from ..benchmarks.rashba import RashbaModel
        from ..benchmarks.rashba_dresselhaus import RashbaDresselhausModel
        m = cfg["model"]
        typ = m["type"]
        if typ == "rashba":
            return RashbaModel(hbar2_over_2m=m.get("hbar2_over_2m", 0.5),
                               alpha=m.get("alpha", 0.1))
        if typ == "rashba_dresselhaus":
            return RashbaDresselhausModel(
                hbar2_over_2m=m.get("hbar2_over_2m", 0.5),
                alpha=m.get("alpha", 0.1), beta=m.get("beta", 0.1))
        raise ValueError(f"unknown analytic model type {typ!r}")
    raise ValueError(f"unknown material kind {kind!r}")
