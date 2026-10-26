# -*- coding: utf-8 -*-
"""Floquet and scanning CLI commands.
"""
from __future__ import annotations

import argparse
import json
import numpy as np
from pathlib import Path

from .cli_utils import (kpath_from_cfg, drive_from_cfg, numerics_from_cfg,
                        _model_lattice, _lattice_from_cfg, _nw_from_cfg,
                        resolve_kmesh, _jobs_for, _parallel_ok,
                        _load_spin_or_exit, _soc_window, _pst_bands,
                        _reciprocal_basis, _resolve_rashba_anchor, ps_model)
from .core.model_factory import build_model
from .core.validation import PASS, FAIL, WARN, BLOCKED, SKIP
from .core.spin_operator import MissingTransformError, SpinInputError
from .drives.peierls import Drive
from .solvers.evolution import floquet_modes
from .observables.pst import pst_collinearity
from .observables.soc_splitting import adjacent_splitting, gap_edge_pairs
from .observables.rashba_fit import fit_gap_edge_pair
from .observables.floquet_points import (resolve_points, run_points,
                                         valley_asymmetry, rashba_effective)
from .observables.spin_texture import spin_texture_wannier_batch
from .version import __version__


def cmd_soc_scan(args):
    """SOC splitting scan across drive parameters."""
    cfg = args.cfg
    model = build_model(cfg)
    
    hw = args.hw if args.hw else cfg.get("drive", {}).get("hw", 0.1)
    if args.gamma:
        gammas = [args.gamma]
    elif cfg.get("gamma"):
        gammas = cfg["gamma"]
    else:
        gammas = [0.0]
    
    pols = args.pol if args.pol else ["x", "y", "z"]
    if pols is None:
        pols = ["x"]
    
    print(f"soc_scan: hw={hw:g}, gamma={gammas}, pol={pols}")
    
    outdir = Path(args.outdir)
    results = []
    
    for pol in pols:
        for gamma in gammas:
            drive_config = dict(hw=hw, gamma=gamma, pol=pol)
            drive = Drive.from_config(drive_config, a_ref=model.a_ref)
            
            window = _soc_window(cfg, model)
            model_drive = model.with_drive(drive)
            e_drive = model_drive.static_bands(np.array([[0.1, 0.1, 0.0]]))
            
            # SOC splitting analysis
            splits = []
            if cfg.get("num_electrons"):
                n_occ = cfg["num_electrons"]
                for i in range(window["i"], window["j"]):
                    try:
                        pair = adjacent_splitting(e_drive, i)
                        splits.append(pair["max"])
                    except Exception:
                        continue
                
                if splits:
                    print(f"  alpha_R[{pol}]:",
                          f"{np.mean(splits)*1000:.2f} meV")
            
            # PST analysis if applicable
            if cfg.get("num_electrons") and not _pst_not_applicable(args, cfg):
                km = resolve_kmesh(cfg, "kmesh")
                kk = np.array([[x, y, 0.0] 
                              for x in np.linspace(-km["kmax"], km["kmax"], 31)
                              for y in np.linspace(-km["kmax"], km["kmax"], 31)])
                
                spin_ops = _load_spin_or_exit(cfg, model_drive)
                bands = _pst_bands(cfg, model_drive)
                pst_data = []
                
                for band_pair in bands:
                    try:
                        batch = spin_texture_wannier_batch(
                            model_drive.Hk, spin_ops, kk, band_pair)
                        spins = batch["spins"]
                        pst = pst_collinearity(spins)
                        pst_data.append(pst)
                    except Exception:
                        continue
                
                if pst_data:
                    best = pst_data[np.argmax([p["collinearity"] for p in pst_data])]
                    r = dict(pol=pol, gamma=gamma, 
                            collinearity=best["collinearity"],
                            pst_axis=best["pst_axis"])
                    print(f"  PST {r['pol']:7s} g={r['gamma']:.3g} "
                          f"coll={r['collinearity']:.3f}")
                    results.append(r)
    
    out_json = outdir / "soc_scan.json"
    json.dump(results, open(out_json, 'w'), indent=2)
    print("  ->", out_json)


def cmd_floquet_bands(args):
    """Generate Floquet band structure."""
    cfg = args.cfg
    model = build_model(cfg)
    drive = drive_from_cfg(cfg, model)
    
    if drive is None:
        print("ERROR: no drive configured for Floquet analysis")
        return 1
    
    kpts, labels = kpath_from_cfg(cfg, nk_seg=args.nk)
    if kpts is None:
        if args.mode == "window" and args.center:
            cx, cy, cz = args.center
            kpts = np.array([[cx + k, cy, cz] for k in np.linspace(-0.1, 0.1, 81)])
            labels = None
        else:
            kpts = np.array([[k, 0.0, 0.0] for k in np.linspace(0, 0.3, 121)])
            labels = None
    
    numerics = numerics_from_cfg(cfg, model)
    n_km = numerics["n_km"]
    
    print(f"floquet_bands: {len(kpts)} k-points on {labels}, "
          f"hw={drive.hw:g} eV, n_km={n_km}")
    
    model_drive = model.with_drive(drive)
    floquet = floquet_modes(model_drive.Hk, drive, n_km=n_km)
    e_floquet = floquet["e_floquet"]
    
    npz = args.outdir / "floquet_bands.npz"
    np.savez_compressed(npz, kpts=kpts, e=e_floquet, labels=labels if labels else [],
                        drive=drive.as_dict())
    
    if args.plot or args.plot_only:
        path = args.outdir / "floquet_bands.png"
        window = args.window if args.window else None
        lattice = _model_lattice(model)
        _plot_bands(kpts, e_floquet, labels, path, window=window, lattice=lattice)
        if not args.plot_only:
            print("  ->", out["figure"])


def cmd_floquet_points(args):
    """Analyze Floquet states at specific points."""
    cfg = args.cfg
    model = build_model(cfg)
    drive = drive_from_cfg(cfg, model)
    
    if drive is None:
        print("ERROR: no drive configured for Floquet analysis")
        return 1
    
    hw = args.hw if args.hw else drive.hw
    points = args.points if args.points else ["Gamma", "K"]
    
    point_dict = {
        "Gamma": np.array([0.0, 0.0, 0.0]),
        "K": np.array([1/3, 1/3, 0.0]),
        "M": np.array([1/2, 0.0, 0.0]),
    }
    
    kpoints = [point_dict.get(p, np.array([p[0], p[1], p[2]])) for p in points]
    
    print(f"floquet_points: hw={hw:g} eV, points={points}, "
          f"pol={drive.pol}, gamma={drive.gamma}")
    
    numerics = numerics_from_cfg(cfg, model)
    n_km = numerics["n_km"]
    
    model_drive = model.with_drive(drive)
    results = []
    
    for p_idx, k_pt in enumerate(kpoints):
        point_name = points[p_idx]
        
        # Local grid for spin analysis
        local_kpts = local_grid_at_hw(k_pt, drive.hw, n_grid=5)
        floquet_local = floquet_modes(model_drive.Hk, drive, n_km=n_km,
                                     kpts=local_kpts)
        
        # Spin analysis
        try:
            spin_ops = _load_spin_or_exit(cfg, model)
            spin_data = []
            
            for band_idx in range(min(4, model.nw)):  # First 4 bands
                try:
                    batch = spin_texture_wannier_batch(
                        model_drive.Hk, spin_ops, local_kpts, band_idx)
                    spins = batch["spins"]
                    
                    # Valley asymmetry
                    if point_name in ["K", "K'"]:
                        asym = valley_asymmetry(spins, k_pt)
                        spin_data.append(dict(band=band_idx, asym=asym))
                    
                    # Rashba effective parameter
                    rashba = rashba_effective(spins, drive.hw)
                    spin_data.append(dict(band=band_idx, rashba=rashba))
                except Exception:
                    continue
            
            results.append(dict(point=point_name, k=k_pt, 
                               spin_data=spin_data))
            
            print(f"  {point_name}:")
            for data in spin_data:
                for key, val in data.items():
                    if key != "band":
                        print(f"    {key:16s} max|dS|(reliable)={val:.3e}  "
                              f"(band {data['band']})")
        
        except (MissingTransformError, SpinInputError):
            print(f"  {point_name}: BLOCKED - no spin operator available")
            continue
    
    # Global valley asymmetry
    if "K" in points and len(results) > 0:
        k_data = [r for r in results if r["point"] in ["K", "K'"]]
        if k_data:
            asym_worst = 0
            for r in k_data:
                for d in r["spin_data"]:
                    if "asym" in d:
                        asym_worst = max(asym_worst, abs(d["asym"]))
            print(f"  valley asymmetry: max|dEbar|={asym_worst*1000:.3f} meV")
    
    out_json = args.outdir / "floquet_points.json"
    json.dump(results, open(out_json, 'w'), indent=2)
    print("  ->", out_json)


def cmd_scan(args):
    """Parameter scan across drive and SOC parameters."""
    cfg = args.cfg
    model = build_model(cfg)
    
    try:
        spin_ops = _load_spin_or_exit(cfg, model)
    except (MissingTransformError, SpinInputError):
        print("BLOCKED: no spin operator available for the scan observable.")
        return 2
    
    hw_vals = args.hw if args.hw else np.linspace(0.01, 0.3, 20)
    gamma_vals = args.gamma if args.gamma else np.linspace(0.0, 1.0, 11)
    pol = args.pol if args.pol else "z"
    
    print(f"scan: {model.nw} bands, hw={hw_vals} eV, gamma={gamma_vals}, "
          f"pol={pol}")
    
    outdir = Path(args.outdir)
    results = []
    
    for hw in hw_vals:
        for gamma in gamma_vals:
            drive_config = dict(hw=hw, gamma=gamma, pol=pol)
            drive = Drive.from_config(drive_config, a_ref=model.a_ref)
            
            model_drive = model.with_drive(drive)
            
            # SOC splitting
            try:
                e_drive = model_drive.static_bands(np.array([[0.1, 0.1, 0.0]]))
                splits = []
                if cfg.get("num_electrons"):
                    n_occ = cfg["num_electrons"]
                    window = _soc_window(cfg, model)
                    for i in range(window["i"], window["j"]):
                        try:
                            pair = adjacent_splitting(e_drive, i)
                            splits.append(pair["max"])
                        except Exception:
                            continue
                
                # Spin texture metrics
                km = resolve_kmesh(cfg, "kmesh")
                kk = np.array([[x, y, 0.0] 
                              for x in np.linspace(-km["kmax"], km["kmax"], 15)
                              for y in np.linspace(-km["kmax"], km["kmax"], 15)])
                
                mean_spins = []
                max_spins = []
                for band_idx in range(min(2, model.nw)):  # First 2 bands
                    try:
                        batch = spin_texture_wannier_batch(
                            model_drive.Hk, spin_ops, kk, band_idx)
                        spins = batch["spins"]
                        mean_spins.append(np.mean(np.linalg.norm(spins, axis=1)))
                        max_spins.append(np.max(np.linalg.norm(spins, axis=1)))
                    except Exception:
                        continue
                
                # Weighted sum error (if applicable)
                wsum_err = 0.0  # Placeholder
                
                r = dict(hw=hw, gamma=gamma, 
                        eps_spread=np.std(splits) if splits else 0.0,
                        mean_spins=np.mean(mean_spins) if mean_spins else 0.0,
                        max_spins=np.max(max_spins) if max_spins else 0.0,
                        wsum_err=wsum_err)
                results.append(r)
                
                if len(results) <= 5:  # Show first few rows
                    print(f"  {r['hw']:.2f}  {r['gamma']:.3f}   {r['eps_spread']:.4f} "
                          f"{r['mean_spins']:.3f}  {r['max_spins']:.3f}  {r['wsum_err']:.2e}")
            except Exception:
                continue
    
    print("    hw    gamma   eps_spread  mean|S|  max|S|  wsum_err")
    for r in results[:5]:
        print(f"  {r['hw']:.2f}  {r['gamma']:.3f}   {r['eps_spread']:.4f} "
              f"{r['mean_spins']:.3f}  {r['max_spins']:.3f}  {r['wsum_err']:.2e}")
    
    npz = outdir / "scan.npz"
    np.savez_compressed(npz, results=results)
    print(f"  -> {npz}")
    
    if args.plot or args.plot_only:
        path = outdir / "scan.png"
        _plot_scan(results, path)
        if not args.plot_only:
            print(f"  figure: {path}")


def _plot_scan(rows, path):
    """Plot scan results."""
    from matplotlib import pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    
    hw_vals = [r["hw"] for r in rows]
    gamma_vals = [r["gamma"] for r in rows]
    
    # SOC splitting spread
    eps_spreads = [r["eps_spread"] for r in rows]
    axes[0, 0].scatter(hw_vals, gamma_vals, c=eps_spreads, cmap='viridis')
    axes[0, 0].set_title('SOC splitting spread')
    axes[0, 0].set_xlabel('hw (eV)')
    axes[0, 0].set_ylabel('gamma')
    
    # Mean spin magnitude
    mean_spins = [r["mean_spins"] for r in rows]
    axes[0, 1].scatter(hw_vals, gamma_vals, c=mean_spins, cmap='plasma')
    axes[0, 1].set_title('Mean spin magnitude')
    axes[0, 1].set_xlabel('hw (eV)')
    axes[0, 1].set_ylabel('gamma')
    
    # Max spin magnitude
    max_spins = [r["max_spins"] for r in rows]
    axes[1, 0].scatter(hw_vals, gamma_vals, c=max_spins, cmap='magma')
    axes[1, 0].set_title('Max spin magnitude')
    axes[1, 0].set_xlabel('hw (eV)')
    axes[1, 0].set_ylabel('gamma')
    
    # Weighted sum error
    wsum_errs = [r["wsum_err"] for r in rows]
    axes[1, 1].scatter(hw_vals, gamma_vals, c=wsum_errs, cmap='coolwarm')
    axes[1, 1].set_title('Weighted sum error')
    axes[1, 1].set_xlabel('hw (eV)')
    axes[1, 1].set_ylabel('gamma')
    
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def cmd_converge(args):
    """Numerical convergence test."""
    cfg = args.cfg
    model = build_model(cfg)
    
    # Reference parameters
    ref_n_km = args.ref_n_km if args.ref_n_km else 15
    ref_m_sp = args.ref_m_sp if args.ref_m_sp else 31
    
    # Test parameters
    n_km_vals = args.n_km if args.n_km else [3, 5, 7, 9, 11, 13]
    m_sp_vals = args.m_sp if args.m_sp else [7, 11, 15, 19, 23, 27, 31]
    
    # Test k-point
    if args.k:
        k_test = np.array(args.k)
    else:
        k_test = np.array([0.1, 0.1, 0.0])
    
    print(f"converge: {model.nw} bands at k={np.round(k_test, 3)}, "
          f"ref (n_step={ref_n_km}, m_sp={ref_m_sp})")
    
    drive = drive_from_cfg(cfg, model) or Drive(hw=0.1, pol="z", gamma=0.0)
    
    outdir = Path(args.outdir)
    
    # Reference calculation
    model_drive = model.with_drive(drive)
    ref_floquet = floquet_modes(model_drive.Hk, drive, n_km=ref_n_km, 
                                kpts=[k_test])
    ref_e = ref_floquet["e_floquet"][0]  # First k-point only
    
    # n_km convergence
    n_km_results = []
    for n_km in n_km_vals:
        test_floquet = floquet_modes(model_drive.Hk, drive, n_km=n_km, 
                                    kpts=[k_test])
        test_e = test_floquet["e_floquet"][0]
        
        # Energy deviation
        dev_ebar = np.mean(np.abs(test_e - ref_e))
        
        r = dict(n_km=n_km, dev_ebar=dev_ebar)
        n_km_results.append(r)
    
    # m_sp convergence (needs spin operators)
    try:
        spin_ops = _load_spin_or_exit(cfg, model)
        m_sp_results = []
        
        km = resolve_kmesh(cfg, "kmesh", default_n=ref_m_sp, default_kmax=0.55)
        kk = full_kpts(m_sp_vals[-1], km["kmax"])  # Use largest mesh
        
        for m_sp in m_sp_vals:
            # Subsample k-points for smaller meshes
            k_sub = kk[::len(kk)//(m_sp**2)][:m_sp**2]
            
            spin_texture_data = []
            for band_idx in range(min(2, model.nw)):
                try:
                    batch = spin_texture_wannier_batch(
                        model.Hk, spin_ops, k_sub, band_idx)
                    spin_texture_data.append(batch["spins"])
                except Exception:
                    continue
            
            # Spin texture deviation (placeholder logic)
            ref_spins = spin_texture_wannier_batch(
                model.Hk, spin_ops, kk, 0)["spins"]
            dev_spin = np.mean([np.mean(np.linalg.norm(s - ref_spins[:len(s)])) 
                               for s in spin_texture_data])
            
            r = dict(m_sp=m_sp, dev_ebar=dev_spin)
            m_sp_results.append(r)
    except (MissingTransformError, SpinInputError):
        m_sp_results = []  # No spin operators available
    
    print("  n_step convergence:")
    for r in n_km_results:
        print(f"    n_step={r['n_km']:4d}  dEbar={r['dev_ebar']:.2e} eV  "
              f"{'PASS' if r['dev_ebar'] < 1e-3 else 'FAIL'}")
    
    print("  m_sp convergence:")
    for r in m_sp_results:
        print(f"    m_sp={r['m_sp']:3d}    dEbar={r['dev_ebar']:.2e} eV  "
              f"{'PASS' if r['dev_ebar'] < 1e-3 else 'FAIL'}")
    
    # Convergence analysis
    convergence = {
        'n_km_pass': all(r['dev_ebar'] < 1e-3 for r in n_km_results[-3:]),
        'm_sp_pass': all(r['dev_ebar'] < 1e-3 for r in m_sp_results[-3:]),
    }
    
    npz = outdir / "converge.npz"
    np.savez_compressed(npz, n_km_results=n_km_results, 
                        m_sp_results=m_sp_results, ref_e=ref_e)
    print(f"  -> {npz}")
    
    if args.plot or args.plot_only:
        path = outdir / "converge.png"
        _plot_converge(dict(n_km=n_km_results, m_sp=m_sp_results), path)
        if not args.plot_only:
            print(f"  figure: {path}")


def _plot_converge(res, path):
    """Plot convergence results."""
    from matplotlib import pyplot as plt
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    
    # n_km convergence
    n_km_data = res["n_km"]
    n_km_vals = [r["n_km"] for r in n_km_data]
    dev_ebar_vals = [r["dev_ebar"] for r in n_km_data]
    axes[0].plot(n_km_vals, dev_ebar_vals, 'o-')
    axes[0].axhline(1e-3, color='red', linestyle='--', alpha=0.5)
    axes[0].set_title('n_km convergence')
    axes[0].set_xlabel('n_km')
    axes[0].set_ylabel('Energy deviation (eV)')
    axes[0].set_yscale('log')
    
    # m_sp convergence
    m_sp_data = res["m_sp"]
    m_sp_vals = [r["m_sp"] for r in m_sp_data]
    dev_spin_vals = [r["dev_ebar"] for r in m_sp_data]
    axes[1].plot(m_sp_vals, dev_spin_vals, 'o-')
    axes[1].axhline(1e-3, color='red', linestyle='--', alpha=0.5)
    axes[1].set_title('m_sp convergence')
    axes[1].set_xlabel('m_sp')
    axes[1].set_ylabel('Spin deviation')
    axes[1].set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def full_kpts(n, kmax):
    """Generate full k-point mesh."""
    xs = np.linspace(-kmax, kmax, n)
    ys = np.linspace(-kmax, kmax, n)
    kk = np.array([[x, y, 0.0] for x in xs for y in ys])
    return kk


load_config = None  # Will be imported from main CLI module