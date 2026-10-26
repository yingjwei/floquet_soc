# -*- coding: utf-8 -*-
"""Basic CLI commands: validate, bands, spin_static.
"""
from __future__ import annotations

import argparse
import json
import numpy as np
from pathlib import Path

from .cli_utils import (kpath_from_cfg, drive_from_cfg, numerics_from_cfg,
                        _spin_status_wannier, _model_lattice, _lattice_from_cfg,
                        _nw_from_cfg, resolve_kmesh, _jobs_for, _parallel_ok,
                        _load_spin_or_exit, _soc_window, _pst_bands,
                        _pst_not_applicable, _reciprocal_basis,
                        _resolve_rashba_anchor, ps_model)
from .core.model_factory import build_model
from .core.validation import ValidationReport, PASS, FAIL, WARN, BLOCKED, SKIP
from .core.spin_operator import MissingTransformError, SpinInputError
from .version import __version__


def cmd_validate(args):
    """Validate model configuration and inputs."""
    cfg = args.cfg
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    rep = ValidationReport(title=f"validate: {cfg.get('name', 'config')}")
    
    # Parse configuration
    try:
        model = build_model(cfg)
        rep.add("model_parse", PASS,
                f"kind={cfg.get('kind')}, num_wann={getattr(model, 'nw', None)}")
    except Exception as exc:                      # noqa: BLE001
        rep.add("model_parse", FAIL, str(exc))
        rep.print()
        rep.save(outdir, "validate")
        return 1

    if cfg.get("kind") == "wannier90":
        win = model.win
        nkpts_win = (len(win["kpoints"]) if win["kpoints"] is not None
                     else None)
        win_src = win.get("source", "wannier90.win")
        rep.add("win_metadata",
                PASS if "reconstructed" not in str(win_src) else WARN,
                f"num_wann={win['num_wann']}, num_bands={win['num_bands']}, "
                f"spinors={win['spinors']}, kpoints={nkpts_win}, "
                f"mp_grid={win['mp_grid']}; source={win_src}")
        
        herr = model.hermiticity_errors()
        folded = herr["HR_hermiticity_folded"]
        if folded is None:
            rep.add("H_R_hermiticity", WARN,
                    f"-R pairing incomplete on the folded R set "
                    f"(n_unpaired={herr['n_unpaired_R']}); H(R) hermiticity "
                    f"is undefined, judged through H(k) instead")
        else:
            rep.add("H_R_hermiticity",
                    PASS if folded < 1e-10 else FAIL,
                    f"max|H(-R)-H(R)^dag| = {folded:.2e} (tol 1e-10)")
        rep.add("Hk_hermiticity",
                PASS if herr["Hk_hermiticity_expanded"] < 1e-10 else FAIL,
                f"max|H(k)-H(k)^dag| = {herr['Hk_hermiticity_expanded']:.2e} "
                f"(tol 1e-10, {len(herr['ktest'])} k probes incl. off-mesh, "
                f"wsvec={'applied' if herr['wsvec_expanded'] else 'no'})")
        
        nd = model.ndegen_consistency()
        rep.add("ndegen_consistency",
                PASS if nd["ndegen_mismatch"] == 0 else FAIL,
                f"ndegen in [{nd['ndegen_min']},{nd['ndegen_max']}], "
                f"mismatched R/-R pairs: {nd['ndegen_mismatch']}")
        rep.add("wannier_centres",
                PASS if "centres.xyz" in model.tau_source else WARN,
                f"tau source = {model.tau_source}")
        
        kpts, labels = kpath_from_cfg(cfg)
        if kpts is not None:
            e = model.static_bands(kpts)
            nocc = cfg.get("num_electrons")
            if nocc:
                gap = float(e[:, nocc].min() - e[:, nocc - 1].max())
                if not np.isfinite(gap):
                    rep.add("static_gap", FAIL, "non-finite gap")
                elif gap <= 0.0:
                    rep.add("static_gap", WARN,
                            f"no gap on the chosen k-path (gap {gap:.4f} eV): "
                            f"metallic / semimetallic (num_electrons={nocc})")
                else:
                    rep.add("static_gap", PASS,
                            f"indirect gap on chosen k-path = {gap:.6f} eV "
                            f"(num_electrons={nocc})")
            else:
                rep.add("static_gap", SKIP, "num_electrons not configured")
        
        _spin_status_wannier(cfg, model, rep)
        
        if cfg.get("drive"):
            try:
                Drive.from_config(cfg["drive"], a_ref=model.a_ref)
                rep.add("drive_config", PASS, "drive parameters parsed")
            except Exception as exc:                      # noqa: BLE001
                rep.add("drive_config", FAIL, str(exc))
    else:
        ops = model.spin_ops()
        rep.add("spin_operator_source", PASS,
                "analytic Pauli spin operators S=sigma/2")
        herm = max(float(np.max(np.abs(model.Hk(k) - model.Hk(k).conj().T)))
                   for k in np.array([[0.1, 0, 0], [0, 0.2, 0], [0.1, 0.1, 0]]))
        rep.add("Hk_hermiticity", PASS if herm < 1e-14 else FAIL,
                f"max|H-H^dag| = {herm:.2e}")

        km = cfg.get("kmesh", dict(n=31, kmax=0.2))
        kk = np.array([[x, y, 0.0] for x in np.linspace(-km["kmax"], km["kmax"], 15)
                       for y in np.linspace(-km["kmax"], km["kmax"], 15)])
        spins = np.array([spin_texture(model.Hk(k), ops)["spins"][0] for k in kk])
        pst = pst_collinearity(spins)
        rep.add("pst_metric", PASS if np.isfinite(pst["collinearity"]) else FAIL,
                f"M eigenvalues = {np.round(pst['eigenvalues'], 4)}")

    rep.print()
    md, js = rep.save(outdir, "validate")
    print(f"\nreport: {md}\n        {js}")
    if rep.count(FAIL):
        return 1
    if rep.count(BLOCKED):
        return 2
    return 0


def cmd_bands(args):
    """Generate static band structure."""
    cfg = args.cfg
    kpts, labels = kpath_from_cfg(cfg, nk_seg=args.nk)
    if kpts is None:
        kpts = np.array([[k, 0.0, 0.0] for k in np.linspace(0, 0.3, 121)])
        labels = None
    
    model = build_model(cfg)
    e = model.static_bands(kpts)
    npz = args.outdir / "bands.npz"
    np.savez_compressed(npz, kpts=kpts, e=e, labels=labels if labels else [])
    print(f"bands: {len(kpts)} k-points x {e.shape[1]} bands -> {npz}")
    print(f"  energy range: [{e.min():.3f}, {e.max():.3f}] eV")
    
    if args.plot or args.plot_only:
        path = args.outdir / "bands.png"
        lattice = _model_lattice(model)
        window = args.window if args.window else None
        _plot_bands(kpts, e, labels, path, window=window, lattice=lattice)
        if not args.plot_only:
            print(f"figure: {path}")


def _plot_bands(kpts, e, labels, path, window=None, lattice=None):
    """Plot band structure."""
    from matplotlib import pyplot as plt
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot bands
    for i in range(e.shape[1]):
        ax.plot(kpts[:, 0], e[:, i], 'k-', linewidth=0.8)
    
    # Add labels
    if labels:
        for (idx, label) in labels.items():
            ax.axvline(kpts[idx, 0], color='gray', linestyle='--', alpha=0.5)
            ax.text(kpts[idx, 0], ax.get_ylim()[1], label, 
                    ha='center', va='bottom')
    
    # Energy window
    if window:
        ax.set_ylim(window)
    
    ax.set_xlabel('k-path')
    ax.set_ylabel('Energy (eV)')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def spin_plot_band(cfg):
    """Determine which band to plot for spin texture."""
    if cfg.get("spin_plot_band"):
        return cfg["spin_plot_band"]
    n_occ = cfg.get("num_electrons")
    if n_occ:
        return n_occ - 1  # VBM
    return 0  # Lowest band


def cmd_spin_static(args):
    """Generate static spin texture."""
    cfg = args.cfg
    outdir = Path(args.outdir)
    model = build_model(cfg)
    
    km = resolve_kmesh(cfg, "kmesh", default_n=21, default_kmax=0.55)
    kk = full_kpts(km["n"], km["kmax"])
    
    try:
        spin_ops = _load_spin_or_exit(cfg, model)
    except (MissingTransformError, SpinInputError) as exc:
        print("BLOCKED: spin texture cannot be computed gauge-consistently.")
        print(str(exc))
        print("missing files:", ", ".join(exc.missing_files))
        return 2
    
    band = spin_plot_band(cfg)
    spin_batch = spin_texture_wannier_batch(model.Hk, spin_ops, kk, band)
    spins = spin_batch["spins"]
    kpts = spin_batch["kpts"]
    
    npz = outdir / "spin_static.npz"
    np.savez_compressed(npz, kpts=kpts, spins=spins, band=band)
    print(f"spin texture: {len(kk)} k-points -> {npz}")
    
    if args.plot or args.plot_only:
        path = outdir / "spin_static.png"
        lattice = _model_lattice(model)
        _plot_spin_texture(kpts, spins, path, title=f"Static spin texture (band {band})", 
                           iband=band, lattice=lattice)
        if not args.plot_only:
            print(f"figure: {path}")


def _plot_spin_texture(kpts, spins, path, title="", iband=0, lattice=None):
    """Plot spin texture."""
    from matplotlib import pyplot as plt
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Plot spin vectors
    ax.quiver(kpts[:, 0], kpts[:, 1], spins[:, 0], spins[:, 1], 
              spins[:, 2], cmap='coolwarm', scale=5)
    
    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


load_config = None  # Will be imported from main CLI module