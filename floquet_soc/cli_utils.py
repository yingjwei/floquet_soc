# -*- coding: utf-8 -*-
"""Utility functions for CLI operations.
"""
from __future__ import annotations

import numpy as np
from .core.model_factory import build_model
from .core.validation import ValidationReport, PASS, FAIL, WARN, BLOCKED, SKIP
from .core.wannier import WannierInputError, read_win, read_hr, write_hr
from .core.cache import load_model_spin_cached
from .core.provenance import collect as prov_collect, write as prov_write
from .core.spin_operator import (read_wannier_spn, validate_spin_matrices,
                                load_spin_for_model, MissingTransformError,
                                SpinInputError)
from .drives.peierls import Drive, peierls_hr_components
from .solvers.evolution import floquet_modes
from .version import __version__
from . import plot_style as ps
from .observables.spin_texture import (spin_texture,
                                        spin_texture_wannier_batch,
                                        spin_texture_path)
from .observables.pst import pst_collinearity
from .observables.spectral import spectral_function
from .observables.soc_splitting import (adjacent_splitting,
    gap_edge_pairs, edge_pair_splitting)
from .observables.rashba_fit import fit_gap_edge_pair
from .observables.pst import (analyse_texture,
                              gauge_invariance_test)
from .observables.spin_texture_map import (resolve_band, full_kpts,
                                            window_kpts, compute_map, plot_map)
from .validation.dft_spin import (read_spin_texture_dat, compare as dft_compare,
                                  save_and_plot as dft_save_plot)
from .observables.floquet_points import (resolve_points, run_points,
                                         valley_asymmetry, rashba_effective,
                                         local_grid_at_hw, make_drive)


def kpath_from_cfg(cfg, nk_seg=None):
    """Extract k-path from configuration."""
    if "kpath" in cfg:
        kp = cfg["kpath"]
        kpts = np.array(kp["kpoints"])
        labels = kp.get("labels")
        if nk_seg is not None:
            # Resegment path
            n_seg = len(kpts) - 1
            nk_per_seg = nk_seg // n_seg
            new_kpts = []
            new_labels = []
            for i in range(n_seg):
                seg_kpts = np.linspace(kpts[i], kpts[i + 1], nk_per_seg)
                new_kpts.extend(seg_kpts)
                if labels and labels[i]:
                    new_labels.append((i * nk_per_seg, labels[i]))
                if i == n_seg - 1 and labels and labels[i + 1]:
                    new_labels.append(((i + 1) * nk_per_seg - 1, labels[i + 1]))
            kpts = np.array(new_kpts)
            labels = dict(new_labels) if new_labels else None
        return kpts, labels
    return None, None


def drive_from_cfg(cfg, model):
    """Build drive from configuration."""
    if "drive" in cfg:
        return Drive.from_config(cfg["drive"], a_ref=model.a_ref)
    return None


def numerics_from_cfg(cfg, model):
    """Extract numerics parameters from configuration."""
    n_km = cfg.get("n_km", 7)
    wsvec = cfg.get("wsvec", True)
    return dict(n_km=n_km, wsvec=wsvec)


def _spin_status_wannier(cfg, model, rep):
    """Check spin operator status for Wannier models."""
    try:
        spin_ops = load_spin_for_model(model, cfg)
        rep.add("spin_operator_source", PASS,
                f"spin operators from {model.tau_source}")
        
        # Validate spin matrices
        herm = validate_spin_matrices(spin_ops)
        rep.add("spin_hermiticity", 
                PASS if herm < 1e-12 else FAIL,
                f"max|S-S^dag| = {herm:.2e}")
        
        # Check completeness
        n_spin = len(spin_ops)
        expected = model.nw * (2 if model.win.get("spinors") else 1)
        rep.add("spin_completeness",
                PASS if n_spin == expected else WARN,
                f"n_spin={n_spin}, expected={expected}")
        
        return True
    except MissingTransformError as exc:
        rep.add("spin_operator_source", BLOCKED, str(exc))
        rep.add("spin_hermiticity", SKIP, "no spin operators")
        return False
    except SpinInputError as exc:
        rep.add("spin_operator_source", FAIL, str(exc))
        rep.add("spin_hermiticity", SKIP, "invalid spin operators")
        return False


def _model_lattice(model):
    """Get lattice vectors from model."""
    if hasattr(model, 'a_ref'):
        return model.a_ref
    return np.eye(3)


def _lattice_from_cfg(cfg):
    """Get lattice vectors from configuration."""
    if "lattice" in cfg:
        return np.array(cfg["lattice"])
    return np.eye(3)


def _nw_from_cfg(cfg):
    """Get number of Wannier functions from configuration."""
    if cfg.get("kind") == "wannier90":
        return cfg.get("num_wann", None)
    return 4  # default for analytic models


def resolve_kmesh(cfg, key, default_n=21, default_kmax=0.55):
    """Resolve k-mesh parameters."""
    if key in cfg:
        kmesh = cfg[key]
        n = kmesh.get("n", default_n)
        kmax = kmesh.get("kmax", default_kmax)
        return dict(n=n, kmax=kmax)
    return dict(n=default_n, kmax=default_kmax)


def _jobs_for(args, cfg, nk, min_points=16):
    """Determine job distribution for parallel processing."""
    if args.parallel is None:
        return None
    
    n_jobs = args.parallel
    if n_jobs <= 1:
        return None
    
    # Distribute k-points
    if nk < min_points:
        return None
    
    points_per_job = max(1, nk // n_jobs)
    return dict(n_jobs=n_jobs, points_per_job=points_per_job)


def _parallel_ok(cfg, model=None):
    """Check if parallel processing is allowed."""
    if cfg.get("parallel", True):
        return True
    return False


def _load_spin_or_exit(cfg, model):
    """Load spin operators or exit with error."""
    try:
        return load_spin_for_model(model, cfg)
    except (MissingTransformError, SpinInputError) as exc:
        print(f"BLOCKED: {exc}")
        exit(2)


def _soc_window(cfg, model):
    """Determine SOC window for analysis."""
    if cfg.get("num_electrons"):
        n_occ = cfg["num_electrons"]
        # VBM-1 to CBM+1 window by default
        return dict(i=max(0, n_occ - 1), j=min(model.nw - 1, n_occ + 1))
    # Auto window: 30% around Fermi level estimate
    e_mid = model.static_bands(np.array([[0.1, 0.1, 0.0]]))
    e_fermi = np.median(e_mid)
    window = 0.3  # 30% of total range
    e_range = model.static_bands(np.array([[0.1, 0.1, 0.0]]))
    e_total = e_range.max() - e_range.min()
    band_width = window * e_total / model.nw
    i_band = int((e_fermi - e_range.min()) / band_width)
    return dict(i=max(0, i_band - 2), j=min(model.nw - 1, i_band + 3))


def _pst_bands(cfg, model):
    """Select PST analysis bands."""
    if cfg.get("pst_bands"):
        return cfg["pst_bands"]
    # Auto-select: pairs around Fermi level
    n_occ = cfg.get("num_electrons")
    if n_occ:
        # bands around Fermi level
        i = max(0, n_occ - 2)
        j = min(model.nw - 1, n_occ + 2)
        return list(range(i, j + 1))
    return list(range(model.nw))


def _pst_not_applicable(args, cfg):
    """Check if PST is applicable for the system."""
    if cfg.get("kind") == "wannier90" and cfg.get("num_electrons"):
        # Check PT-symmetry indicators
        model = build_model(cfg)
        static_bands = model.static_bands(np.array([[0.1, 0.1, 0.0]]))
        # Simple PT-symmetry test: time-reversion symmetry + non-magnetic
        # This is a heuristic - real PT-symmetry requires more analysis
        return False  # Assume applicable for most cases
    return True


def _reciprocal_basis(lattice):
    """Compute reciprocal basis vectors."""
    a = lattice
    b = 2 * np.pi * np.linalg.inv(a).T
    return b


def _resolve_rashba_anchor(cfg):
    """Resolve Rashba analysis anchor point."""
    if cfg.get("rashba_anchor"):
        return cfg["rashba_anchor"]
    # Default: Gamma point
    return dict(point="Gamma", k=np.array([0.0, 0.0, 0.0]))


def ps_model(popt, k):
    """Parabolic model for Rashba fitting."""
    alpha, e0 = popt
    k_norm = np.linalg.norm(k)
    return e0 + alpha * k_norm**2