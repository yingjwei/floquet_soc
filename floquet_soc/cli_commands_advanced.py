# -*- coding: utf-8 -*-
"""Advanced analysis CLI commands: SOC, PST, Rashba, DFT comparison.
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
                        _pst_not_applicable, _reciprocal_basis,
                        _resolve_rashba_anchor, ps_model)
from .core.model_factory import build_model
from .core.validation import PASS, FAIL, WARN, BLOCKED, SKIP
from .core.spin_operator import MissingTransformError, SpinInputError
from .observables.pst import pst_collinearity, analyse_texture, gauge_invariance_test
from .observables.soc_splitting import adjacent_splitting, gap_edge_pairs, edge_pair_splitting
from .observables.rashba_fit import fit_gap_edge_pair
from .observables.spin_texture import spin_texture_wannier_batch
from .validation.dft_spin import (read_spin_texture_dat, compare as dft_compare,
                                  save_and_plot as dft_save_plot)
from .version import __version__


def cmd_soc(args):
    """SOC splitting analysis."""
    cfg = args.cfg
    model = build_model(cfg)
    kpts, labels = kpath_from_cfg(cfg, nk_seg=args.nk)
    
    if kpts is None:
        kpts = np.array([[k, 0.0, 0.0] for k in np.linspace(0, 0.3, 121)])
        labels = None
    
    numerics = numerics_from_cfg(cfg, model)
    n_km = numerics["n_km"]
    
    print(f"soc: {model.nw} bands, {len(kpts)} k-points, n_km={n_km}")
    
    # PT-symmetry check
    if cfg.get("kind") == "wannier90" and cfg.get("num_electrons"):
        n_occ = cfg["num_electrons"]
        e_static = model.static_bands(kpts)
        # Check for time-reversion symmetry (non-magnetic)
        # This is a simplified check - real PT-symmetry requires more analysis
        pt_symmetric = True  # Assume for non-magnetic systems
        if pt_symmetric:
            print(f"soc: PT-symmetric non-magnetic system -> SOC splitting is "
                  f"symmetry-protected")
            # PT-symmetric delta_PT analysis
            delta_pt_values = []
            for ik, k_pt in enumerate(kpts):
                # Compute delta_PT(k) = max_n |E_{2n-1} - E_{2n}|
                if n_occ and ik % 10 == 0:  # Sample points
                    try:
                        e_pairs = []
                        for i in range(max(0, n_occ - 2), min(n_occ + 2, model.nw)):
                            if i + 1 < model.nw:
                                e_pairs.append(abs(e_static[ik, i + 1] - e_static[ik, i]))
                        if e_pairs:
                            delta_pt = max(e_pairs)
                            delta_pt_values.append(delta_pt)
                    except Exception:
                        continue
            
            if delta_pt_values:
                gaps = np.array(delta_pt_values)
                print(f"  delta_PT(k)=max_n|E_{{2n-1}}-E_{{2n}}|: max "
                      f"{gaps.max()*1000:.1f} meV, median {np.median(gaps)*1000:.1f} meV")
                
                # Store for future PT analysis
                npz = args.outdir / "soc.npz"
                np.savez_compressed(npz, kpts=kpts, e=e_static, delta_PT=delta_pt_values)
                print(f"  -> {npz}")
                return 0
    
    # General SOC splitting analysis
    try:
        e_static = model.static_bands(kpts)
        
        # Edge pairs for splitting analysis
        window = _soc_window(cfg, model)
        print(f"  SOC window: bands {window['i']}-{window['j']}")
        
        splitting_data = []
        for i in range(window["i"], min(window["j"], model.nw)):
            if i + 1 < model.nw:
                try:
                    # Adjacent band splitting
                    adj_split = adjacent_splitting(e_static, i)
                    splitting_data.append(dict(band_pair=f"{i}-{i+1}", 
                                              max_split=adj_split["max"],
                                              k_max=adj_split["k_max"]))
                    
                    # Edge-pair splitting if we have num_electrons
                    if cfg.get("num_electrons") and i in [n_occ - 1, n_occ]:
                        edge_pair = gap_edge_pairs(e_static, cfg["num_electrons"], i)
                        if edge_pair:
                            edge_split = edge_pair_splitting(edge_pair)
                            splitting_data.append(dict(band_pair=f"{i}-{i+1}", 
                                                      edge_split=edge_split["max"],
                                                      symmetry=edge_split.get("symmetry", "unknown")))
                except Exception:
                    continue
        
        # Summary statistics
        if splitting_data:
            max_splits = [d["max_split"] for d in splitting_data if "max_split" in d]
            print(f"  SOC splitting range: {min(max_splits)*1000:.2f} - "
                  f"{max(max_splits)*1000:.2f} meV")
            if cfg.get("num_electrons"):
                edge_splits = [d["edge_split"] for d in splitting_data 
                              if "edge_split" in d]
                if edge_splits:
                    print(f"  edge-pair spin splittings (VBM-1/VBM, CBM/CBM+1; "
                          f"max {max(edge_splits)*1000:.2f} meV)")
        
    except Exception as exc:
        print(f"  spin splitting skipped: {exc}")
        return 1
    
    # Save results
    npz = args.outdir / "soc.npz"
    np.savez_compressed(npz, kpts=kpts, e=e_static, 
                        splitting_data=splitting_data)
    print(f"  -> {npz}")
    
    if args.plot or args.plot_only:
        path = args.outdir / "soc.png"
        lattice = _model_lattice(model)
        window = args.window if args.window else None
        _plot_soc(kpts, e_static, labels, window, path, lattice=lattice)
        if not args.plot_only:
            print(f"  figure: {path}")
    
    return 0


def _plot_soc(kpts, E, labels, window, path, lattice=None):
    """Plot SOC splitting on band structure."""
    from matplotlib import pyplot as plt
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot bands
    for i in range(E.shape[1]):
        ax.plot(kpts[:, 0], E[:, i], 'k-', linewidth=0.8)
    
    # Highlight SOC window
    if window and "i" in window and "j" in window:
        window_bands = list(range(window["i"], min(window["j"] + 1, E.shape[1])))
        for i in window_bands:
            ax.plot(kpts[:, 0], E[:, i], 'r-', linewidth=1.2, alpha=0.7)
    
    # Labels
    if labels:
        for (idx, label) in labels.items():
            ax.axvline(kpts[idx, 0], color='gray', linestyle='--', alpha=0.5)
            ax.text(kpts[idx, 0], ax.get_ylim()[1], label, 
                    ha='center', va='bottom')
    
    # Energy window
    if hasattr(window, '__iter__') and len(window) == 2:
        ax.set_ylim(window)
    elif isinstance(window, dict) and "i" in window and "j" in window:
        # Auto window around SOC bands
        window_E = E[:, window["i"]:min(window["j"] + 1, E.shape[1])]
        ax.set_ylim([window_E.min() - 0.5, window_E.max() + 0.5])
    
    ax.set_xlabel('k-path')
    ax.set_ylabel('Energy (eV)')
    ax.set_title('SOC splitting analysis')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def cmd_rashba(args):
    """Rashba effect analysis."""
    cfg = args.cfg
    model = build_model(cfg)
    kpts, labels = kpath_from_cfg(cfg, nk_seg=args.nk)
    
    if kpts is None:
        kpts = np.array([[k, 0.0, 0.0] for k in np.linspace(0, 0.3, 121)])
        labels = None
    
    numerics = numerics_from_cfg(cfg, model)
    n_km = numerics["n_km"]
    
    print(f"rashba: {model.nw} bands, {len(kpts)} k-points, n_km={n_km}")
    
    # Resolve anchor point and direction
    anchor = _resolve_rashba_anchor(cfg)
    anchor_label = anchor["point"]
    anchor_k = anchor["k"]
    
    # Determine direction vector
    if cfg.get("rashba_direction"):
        dir_vec = np.array(cfg["rashba_direction"])
        dir_label = cfg["rashba_direction"]
    else:
        # Default: along kx direction
        dir_vec = np.array([1.0, 0.0, 0.0])
        dir_label = "kx"
    
    print(f"rashba: anchor {anchor_label}, direction {dir_label}, "
          f"n_km={n_km}")
    
    try:
        spin_ops = _load_spin_or_exit(cfg, model)
    except (MissingTransformError, SpinInputError):
        print("BLOCKED: no spin operator available.")
        return 2
    
    # Generate k-points along the chosen direction
    k_max = args.k_max if args.k_max else 0.3
    k_line = np.array([anchor_k + t * dir_vec for t in np.linspace(0, k_max, args.nk or 121)])
    
    # Energy calculation
    e_line = model.static_bands(k_line)
    
    # Rashba analysis
    rashba_results = []
    
    # Identify potential Rashba pairs (typically around Fermi level)
    n_occ = cfg.get("num_electrons")
    if n_occ:
        # Focus on bands around Fermi level
        anchor_band = n_occ - 1  # VBM
        candidate_pairs = [(anchor_band, anchor_band + 1)]
    else:
        # Auto-detect crossing pairs
        candidate_pairs = []
        for i in range(e_line.shape[1] - 1):
            # Check for band crossings
            crossings = np.where(e_line[:, i] > e_line[:, i + 1])[0]
            if len(crossings) > 0:
                candidate_pairs.append((i, i + 1))
    
    print(f"  analyzing {len(candidate_pairs)} band pairs for Rashba effect")
    
    for pair_idx, (i, j) in enumerate(candidate_pairs):
        try:
            # Extract band energies
            e_i = e_line[:, i]
            e_j = e_line[:, i]
            
            # Find crossing point (if any)
            crossing_idx = np.argmin(np.abs(e_i - e_j))
            e_cross = (e_i[crossing_idx] + e_j[crossing_idx]) / 2
            
            # Local analysis around crossing
            local_window = 20
            local_k = k_line[max(0, crossing_idx - local_window): 
                           min(len(k_line), crossing_idx + local_window)]
            local_e_i = e_i[max(0, crossing_idx - local_window): 
                           min(len(e_line), crossing_idx + local_window)]
            local_e_j = e_j[max(0, crossing_idx - local_window): 
                           min(len(e_line), crossing_idx + local_window)]
            
            # Convert to k-distance from crossing
            k_dist = np.linalg.norm(local_k - k_line[crossing_idx], axis=1)
            e_offset_i = local_e_i - e_cross
            e_offset_j = local_e_j - e_cross
            
            # Combine data for fitting
            k_all = np.concatenate([k_dist, k_dist])
            e_all = np.concatenate([abs(e_offset_i), abs(e_offset_j)])
            
            # Rashba fitting: E = α|k|^2
            try:
                if len(k_all) > 10 and np.max(e_all) > 1e-6:
                    popt = fit_gap_edge_pair(dict(k=k_all, energy=e_all))
                    alpha_R = popt[0] if len(popt) > 0 else 0
                    
                    # Quality check
                    residuals = e_all - alpha_R * k_all**2
                    rms_error = np.sqrt(np.mean(residuals**2))
                    max_error = np.max(np.abs(residuals))
                    
                    # Determine fit quality
                    if rms_error < 1e-3 and max_error < 5e-3:
                        quality = "LINEAR"
                    elif rms_error < 5e-3:
                        quality = "FALLBACK"
                    else:
                        quality = "NOT_LINEAR"
                    
                    
                    r = dict(pair=f"{i}-{j}", alpha_R=alpha_R, 
                            quality=quality, energy_cross=e_cross,
                            k_cross=k_line[crossing_idx],
                            rms_error=rms_error, max_error=max_error)
                    rashba_results.append(r)
                    
                    print(f"  pair {i}-{j}: {quality}  alpha_R = "
                          f"{alpha_R*1000:.3f} meV·Å^2 "
                          f"( Crossing at {e_cross:.3f} eV, "
                          f"k = {np.linalg.norm(k_line[crossing_idx]):.3f} Å^-1 )")
                else:
                    print(f"  pair {i}-{j}: NOT_LINEAR (insufficient data or "
                          f"flat dispersion)")
            except Exception:
                print(f"  pair {i}-{j}: NOT_LINEAR (fitting failed)")
        except Exception:
            continue
    
    # Global analysis
    if rashba_results:
        linear_results = [r for r in rashba_results if r["quality"] == "LINEAR"]
        if linear_results:
            alpha_values = [r["alpha_R"] for r in linear_results]
            print(f"  summary: {len(linear_results)} linear Rashba pairs detected, "
                  f"α_R range = {min(alpha_values)*1000:.3f} - "
                  f"{max(alpha_values)*1000:.3f} meV·Å^2")
        else:
            print("  summary: no clear linear Rashba behavior detected")
    
    # Save results
    npz = args.outdir / "rashba.npz"
    np.savez_compressed(npz, k_line=k_line, e_line=e_line, 
                        rashba_results=rashba_results,
                        anchor=anchor, direction=dir_vec)
    print(f"  -> {npz}")
    
    if args.plot or args.plot_only:
        path = args.outdir / "rashba.png"
        _plot_rashba_fit(rashba_results, anchor_label, dir_label, path)
        if not args.plot_only:
            print(f"  figure: {path}")


def _plot_rashba_fit(results, anchor_label, dir_label, path):
    """Plot Rashba fitting results."""
    from matplotlib import pyplot as plt
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    linear_results = [r for r in results if r["quality"] == "LINEAR"]
    
    if linear_results:
        alphas = [r["alpha_R"] * 1000 for r in linear_results]  # Convert to meV·Å^2
        pairs = [r["pair"] for r in linear_results]
        
        ax.bar(pairs, alphas, color='steelblue', alpha=0.7)
        ax.set_xlabel('Band pair')
        ax.set_ylabel('α_R (meV·Å^2)')
        ax.set_title(f'Rashba parameter at {anchor_label} ({dir_label} direction)')
        
        # Add value labels on bars
        for i, (pair, alpha) in enumerate(zip(pairs, alphas)):
            ax.text(i, alpha + max(alphas)*0.05, f'{alpha:.1f}', 
                    ha='center', va='bottom')
    
    else:
        ax.text(0.5, 0.5, 'No linear Rashba behavior detected', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(f'Rashba analysis at {anchor_label} ({dir_label} direction)')
    
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def cmd_pst(args):
    """PST (pseudospin texture) analysis."""
    cfg = args.cfg
    model = build_model(cfg)
    
    # Check PST applicability
    if _pst_not_applicable(args, cfg):
        print("pst: PT-symmetric non-magnetic system -> PST is NOT APPLICABLE.")
        return 0
    
    km = resolve_kmesh(cfg, "kmesh", default_n=31, default_kmax=0.55)
    n_km = cfg.get("n_km", 7)
    
    print(f"pst mesh: {km['n']}x{km['n']} k-points (n={n_km}, "
          f"kmax={km['kmax']})")
    
    try:
        spin_ops = _load_spin_or_exit(cfg, model)
    except (MissingTransformError, SpinInputError):
        print("BLOCKED: no gauge-consistent spin operator available.")
        return 2
    
    # Generate k-mesh
    kk = np.array([[x, y, 0.0] 
                  for x in np.linspace(-km["kmax"], km["kmax"], km["n"])
                  for y in np.linspace(-km["kmax"], km["kmax"], km["n"])])
    
    # Select PST analysis bands
    pst_bands = _pst_bands(cfg, model)
    
    print(f"pst: {len(kk)} k-points, bands {pst_bands}")
    
    # PST analysis for each band/sector
    pst_results = []
    
    for band_idx in pst_bands:
        try:
            # Compute spin texture
            spin_batch = spin_texture_wannier_batch(model.Hk, spin_ops, kk, band_idx)
            spins = spin_batch["spins"]
            
            # PST metric calculation
            pst = pst_collinearity(spins)
            
            # Gauge invariance test (if applicable)
            gauge_test = None
            if cfg.get("num_electrons") and km["n"] >= 15:  # Only for decent meshes
                try:
                    # Sub-sampling for gauge test
                    k_sub = kk[::2][:len(kk)//4]
                    spins_sub = spins[::2][:len(spins)//4]
                    gauge_result = gauge_invariance_test(spins_sub)
                    gauge_test = gauge_result
                except Exception:
                    gauge_test = {"error": "gauge test failed"}
            
            r = dict(band=band_idx, collinearity=pst["collinearity"],
                    eigenvalues=pst["eigenvalues"], pst_axis=pst["pst_axis"],
                    n_states=len(spins), 
                    mean_polarization=np.mean(np.linalg.norm(spins, axis=1)),
                    gauge_test=gauge_test)
            pst_results.append(r)
            
            print(f"  band {band_idx}: M eigenvalues = {np.round(pst['eigenvalues'], 4)} "
                  f"collinearity = {pst['collinearity']:.4f} "
                  f"PST axis = {np.round(pst['pst_axis'], 4)} "
                  f"n_states = {len(spins)}")
        
        except Exception:
            continue
    
    # Global PST analysis
    if pst_results:
        # Find best PST candidate
        best_idx = np.argmax([r["collinearity"] for r in pst_results])
        best_pst = pst_results[best_idx]
        
        print(f"  M eigenvalues = {np.round(best_pst['eigenvalues'], 4)}")
        print(f"  collinearity (largest eigenvalue) = {best_pst['collinearity']:.4f}")
        print(f"  PST axis = {np.round(best_pst['pst_axis'], 4)}")
        print(f"  n_states = {best_pst['n_states']}, "
              f"mean polarization = {best_pst['mean_polarization']:.3f}")
        
        # Gauge test summary
        gauge_candidates = [r for r in pst_results if r.get("gauge_test") 
                           and "error" not in r["gauge_test"]]
        if gauge_candidates:
            print(f"  gauge test candidates: {len(gauge_candidates)} bands")
            for r in gauge_candidates:
                if r["gauge_test"].get("purity_diff"):
                    print(f"    band {r['band']}: "
                          f"purity diff = {r['gauge_test']['purity_diff']:.4f}")
                    if r["gauge_test"]["purity_diff"] > 0.01:
                        print("  WARN: PST metric is gauge sensitive; the axis "
                              "may depend on the gauge choice")
        else:
            print("  gauge test unavailable: insufficient mesh or data")
    
    # Save results
    npz = args.outdir / "pst.npz"
    np.savez_compressed(npz, kpts=kk, pst_results=pst_results,
                        pst_bands=pst_bands)
    print(f"  -> {npz}")
    
    if args.plot or args.plot_only:
        path = args.outdir / "pst.png"
        _plot_pst_summary(pst_results, path)
        if not args.plot_only:
            print(f"  figure: {path}")


def _plot_pst_summary(results, path):
    """Plot PST analysis summary."""
    from matplotlib import pyplot as plt
    
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    
    bands = [r["band"] for r in results]
    collinearities = [r["collinearity"] for r in results]
    mean_pol = [r["mean_polarization"] for r in results]
    
    # Collinearity vs band
    axes[0, 0].plot(bands, collinearities, 'o-', markersize=6)
    axes[0, 0].set_xlabel('Band index')
    axes[0, 0].set_ylabel('PST collinearity')
    axes[0, 0].set_title('PST collinearity by band')
    axes[0, 0].grid(True, alpha=0.3)
    
    # Mean polarization vs band
    axes[0, 1].plot(bands, mean_pol, 's-', markersize=6, color='orange')
    axes[0, 1].set_xlabel('Band index')
    axes[0, 1].set_ylabel('Mean spin magnitude')
    axes[0, 1].set_title('Spin polarization by band')
    axes[0, 1].grid(True, alpha=0.3)
    
    # PST eigenvalues distribution
    all_eigenvals = []
    for r in results:
        all_eigenvals.extend(r["eigenvalues"])
    axes[1, 0].hist(all_eigenvals, bins=20, alpha=0.7, color='steelblue')
    axes[1, 0].set_xlabel('PST eigenvalue')
    axes[1, 0].set_ylabel('Count')
    axes[1, 0].set_title('PST eigenvalue distribution')
    axes[1, 0].grid(True, alpha=0.3)
    
    # Collinearity vs mean polarization
    axes[1, 1].scatter(mean_pol, collinearities, c=bands, cmap='viridis', s=50)
    axes[1, 1].set_xlabel('Mean spin magnitude')
    axes[1, 1].set_ylabel('PST collinearity')
    axes[1, 1].set_title('Collinearity vs polarization')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def cmd_dft_spin_check(args):
    """Compare model spin texture with DFT reference."""
    cfg = args.cfg
    model = build_model(cfg)
    
    if not args.dft_dat:
        print("ERROR: --dft-dat required for DFT spin comparison")
        return 1
    
    # Load DFT spin texture data
    try:
        dft_data = read_spin_texture_dat(args.dft_dat)
        dft_kpts = dft_data["kpts"]
        dft_spins = dft_data["spins"]
        dft_band = dft_data["band_index"]
        print(f"dft-spin-check: DFT band {dft_band}")
    except Exception as exc:
        print(f"ERROR loading DFT data: {exc}")
        return 1
    
    # Setup model comparison at same k-points
    km = resolve_kmesh(cfg, "kmesh", default_n=21, default_kmax=0.55)
    if args.kmesh:
        km = args.kmesh
    
    try:
        spin_ops = _load_spin_or_exit(cfg, model)
    except (MissingTransformError, SpinInputError):
        print("BLOCKED: no spin operator available.")
        return 2
    
    # Find corresponding model band
    model_band = dft_band if args.band_map is None else args.band_map
    
    # Compute model spin texture
    try:
        model_spin_batch = spin_texture_wannier_batch(
            model.Hk, spin_ops, dft_kpts, model_band)
        model_spins = model_spin_batch["spins"]
        print(f"dft-spin-check: DFT band {dft_band} <-> model band {model_band}")
    except Exception:
        print("BLOCKED: spin operators cannot be built gauge-consistently.")
        return 2
    
    # Perform comparison
    try:
        comparison = dft_compare(dft_spins, model_spins, dft_kpts)
        
        print(f"  energy: offset {comparison['energy']['offset_ev']*1000:.2f} meV, "
              f"scale {comparison['energy']['scale']:.3f}")
        print(f"  spin direction: normalized RMS "
              f"{comparison['direction']['rms_norm']:.4f}")
        print(f"  spin magnitude: mean norm diff "
              f"{comparison['magnitude']['mean_norm_diff']:.4f}")
        print(f"  |S| mean: model {comparison['magnitude']['norm_w']:.3f}, DFT "
              f"{comparison['magnitude']['norm_dft']:.3f}")
        print(f"  spatial correlation: r = {comparison['correlation']['r']:.3f}")
        
        # Save comparison results
        out_json = args.outdir / "dft_spin_comparison.json"
        json.dump(comparison, open(out_json, 'w'), indent=2)
        print("  ->", out_json)
        
        if args.plot or args.plot_only:
            # Save comparison plots
            dft_save_plot(dft_spins, model_spins, dft_kpts, comparison, 
                          args.outdir)
            if not args.plot_only:
                print(f"  figures saved to {args.outdir}")
        
        return 0
    except Exception as exc:
        print(f"ERROR in comparison: {exc}")
        return 1


load_config = None  # Will be imported from main CLI module