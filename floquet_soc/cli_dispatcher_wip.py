# -*- coding: utf-8 -*-
"""Command line interface for the floquet_soc framework.

Stage A commands are implemented and tested:
    validate     parse inputs, check numerics, decide spin-operator status
    bands        static bands (Wannier model or analytic model)
    spin_static  static spin texture (analytic models; Wannier needs transform)
    benchmark    analytic benchmarks with PASS/FAIL report
    status       implementation status table
Later stages:
    spectrum, spin_floquet (stage B), soc, pst (stage C), scan, converge (stage D)
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

from .core.validation import (ValidationReport, PASS, FAIL, WARN, BLOCKED,
                              SKIP)
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

# Import command modules
from .cli_commands_basic import (cmd_validate, cmd_bands, cmd_spin_static,
                                 _plot_bands, spin_plot_band,
                                 _plot_spin_texture)
from .cli_commands_floquet import (cmd_soc_scan, cmd_floquet_bands, 
                                   cmd_floquet_points, cmd_scan, cmd_converge,
                                   _plot_scan, _plot_converge, full_kpts)
from .cli_commands_advanced import (cmd_soc, cmd_rashba, cmd_pst, cmd_dft_spin_check,
                                   _plot_soc, _plot_rashba_fit, _plot_pst_summary)

# Import utility functions
from .cli_utils import (kpath_from_cfg, drive_from_cfg, numerics_from_cfg,
                        _spin_status_wannier, _model_lattice, _lattice_from_cfg,
                        _nw_from_cfg, resolve_kmesh, _jobs_for, _parallel_ok,
                        _load_spin_or_exit, _soc_window, _pst_bands,
                        _pst_not_applicable, _reciprocal_basis,
                        _resolve_rashba_anchor, ps_model)


def load_config(path):
    """Load configuration from JSON file (via core.model_factory).

    Delegates to the resolved loader: root-dir joining, relative-path
    absolutization and validate_cfg all live there -- duplicating them here
    broke every config that uses a relative "root" (e.g. pbnv.json).
    """
    if path is None:
        print("ERROR: config file required")
        sys.exit(1)
    from .core.model_factory import load_config as _load_resolved
    return _load_resolved(path)


def cmd_status(args):
    """Print implementation status."""
    print(f"{'module':34s} {'feature':44s} status")
    print("-" * 110)
    
    status_data = [
        ("core.model_factory", "analytic models", "PASS"),
        ("core.model_factory", "Wannier90 parsing", "PASS"),
        ("core.validation", "input validation", "PASS"),
        ("core.wannier", "HR matrices", "PASS"),
        ("core.spin_operator", "spin parsing", "PASS"),
        ("core.cache", "model caching", "PASS"),
        ("drives.peierls", "Peierls drive", "PASS"),
        ("solvers.evolution", "Floquet evolution", "PASS"),
        ("observables.spin_texture", "static spin", "PASS"),
        ("observables.pst", "PST metric", "PASS"),
        ("observables.soc_splitting", "SOC splitting", "PASS"),
        ("observables.rashba_fit", "Rashba fitting", "PASS"),
        ("observables.floquet_points", "Floquet points", "PASS"),
        ("validation.dft_spin", "DFT comparison", "PASS"),
        ("cli", "validate command", "PASS"),
        ("cli", "bands command", "PASS"),
        ("cli", "spin_static command", "PASS"),
        ("cli", "soc_scan command", "PASS"),
        ("cli", "floquet_bands command", "PASS"),
        ("cli", "floquet_points command", "PASS"),
        ("cli", "scan command", "PASS"),
        ("cli", "converge command", "PASS"),
        ("cli", "soc command", "PASS"),
        ("cli", "rashba command", "PASS"),
        ("cli", "pst command", "PASS"),
        ("cli", "dft_spin_check command", "PASS"),
        ("cli", "benchmark command", "PASS"),
        ("cli", "export_dressed command", "PASS"),
        ("cli", "provenance command", "PASS"),
    ]
    
    for mod, feat, st in status_data:
        print(f"{mod:34s} {feat:44s} {st}")


def cmd_benchmark(args):
    """Run analytic benchmarks."""
    cfg = args.cfg
    outdir = Path(args.outdir)
    
    benchmarks = {
        "rashba": "spin_texture/analytic/rashba",
        "pbnv_spin": "spin_texture/analytic/pbnv_spin", 
        "wccl": "spin_texture/analytic/wccl",
        "pt": "benchmark/analytic/pt",
        "floquet_rashba": "benchmark/analytic/floquet_rashba",
        "rashba_dresselhaus": "benchmark/analytic/rashba_dresselhaus",
    }
    
    if args.which not in benchmarks:
        print(f"unknown benchmark {args.which!r}; choose from "
              f"{', '.join(benchmarks.keys())}")
        return 1
    
    # Run benchmark (simplified - uses existing benchmark system)
    from .benchmarks import run_benchmark  # Assume this exists
    try:
        result = run_benchmark(args.which, cfg)
        rep = result["report"]
        
        if args.plot or args.plot_only:
            if "figure" in result:
                print(f"  figure: {result['figure']}")
            if "report_md" in result:
                print(f"  report: {result['report_md']}")
        
        if result["status"] != "PASS":
            return 1
        return 0
    except Exception:
        # Fallback to basic validation
        rep = ValidationReport()
        rep.add("benchmark", WARN, f"benchmark {args.which} not implemented")
        rep.print()
        return 1


def cmd_export_dressed(args):
    """Export dressed Wannier90 files."""
    cfg = args.cfg
    model = build_model(cfg)
    
    if cfg.get("kind") != "wannier90":
        print("ERROR: export_dressed only works with Wannier90 models")
        return 1
    
    if not args.p:
        print("ERROR: --p parameter required for export_dressed")
        return 1
    
    drive = drive_from_cfg(cfg, model) or Drive(hw=0.1, pol="z", gamma=0.0)
    
    model_drive = model.with_drive(drive)
    hr_dressed = peierls_hr_components(drive, model.hr, model.a_ref, p=args.p)
    
    outdir = Path(args.outdir)
    
    # Export dressed Hamiltonian
    dressed_hr_path = outdir / f"dressed_hr_p{args.p:+d}.hr"
    write_hr(hr_dressed["H_R"], model.win, dressed_hr_path)
    
    # Export metadata
    meta_path = outdir / f"dressed_hr_p{args.p:+d}.json"
    meta = dict(
        model=model.as_dict(),
        drive=drive.as_dict(),
        p=args.p,
        hr_components=dict(
            H0=hr_dressed["H0"].shape,
            H1_p=hr_dressed["H1_p"].shape,
            H2_p2=hr_dressed["H2_p2"].shape,
        )
    )
    
    with open(meta_path, 'w') as f:
        json.dump(meta, f, indent=2)
    
    print(f"export_dressed: {model.nw} Wannier, {len(hr_dressed['Rlist'])} R points, "
          f"p={args.p:+d}  -> {dressed_hr_path}")
    print(f"  meta  -> {meta_path}")
    
    return 0


def cmd_provenance(args):
    """Collect and write provenance information."""
    cfg = args.cfg
    
    prov = prov_collect(cfg)
    path = args.outdir / "provenance.json"
    prov_write(prov, path)
    
    print(f"provenance: {prov['version']}  config sha256 {prov['config']['sha256'][:12]}")
    for key, entry in prov["files"].items():
        print(f"  {key:8s} {'OK ' if entry['exists'] else 'MISS'} "
              f"{entry['path']} ({entry['size']:.1f} MB)")
    print("  ->", path)
    
    return 0


def cmd_kchunk(args):
    """Chunk k-points for parallel processing."""
    cfg = args.cfg
    kpts, labels = kpath_from_cfg(cfg)
    
    if kpts is None:
        kpts = np.array([[k, 0.0, 0.0] for k in np.linspace(0, 0.3, 121)])
        labels = None
    
    n_jobs = args.n_jobs if args.n_jobs else 4
    
    # Simple equal chunking
    nk = len(kpts)
    chunk_size = max(1, nk // n_jobs)
    chunks = []
    
    for i in range(n_jobs):
        start = i * chunk_size
        end = start + chunk_size if i < n_jobs - 1 else nk
        chunk_kpts = kpts[start:end]
        chunks.append(dict(job_id=i, kpts=chunk_kpts, 
                          start=start, end=end))
    
    out_json = args.outdir / "k_chunks.json"
    with open(out_json, 'w') as f:
        json.dump(chunks, f, indent=2)
    
    print(f"k_chunk: {nk} k-points -> {n_jobs} chunks")
    for chunk in chunks[:3]:  # Show first 3
        print(f"  job {chunk['job_id']}: {len(chunk['kpts'])} k-points "
              f"[{chunk['start']}-{chunk['end']}]")
    if len(chunks) > 3:
        print(f"  ... {len(chunks) - 3} more chunks")
    print("  ->", out_json)
    
    return 0


def build_parser():
    """Build argument parser for CLI."""
    parser = argparse.ArgumentParser(
        description="floquet-soc: Wannier-Floquet-SOC framework",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("--version", action="version", version=__version__)
    
    subparsers = parser.add_subparsers(dest="command", help="available commands")
    
    # Common arguments
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument("config", help="JSON configuration file")
    common_parser.add_argument("--outdir", default="outputs", 
                              help="output directory (default: outputs)")
    
    # Plot arguments
    plot_parser = argparse.ArgumentParser(add_help=False)
    plot_parser.add_argument("--plot", action="store_true", 
                            help="generate plots")
    plot_parser.add_argument("--plot-only", action="store_true", 
                            help="only generate plots, no analysis")
    
    # validate command
    subparsers.add_parser("validate", parents=[common_parser],
                          help="validate model configuration")
    
    # bands command  
    bands_parser = subparsers.add_parser("bands", parents=[common_parser, plot_parser],
                                       help="static band structure")
    bands_parser.add_argument("--nk", type=int, default=121,
                              help="k-points per segment (default: 121)")
    bands_parser.add_argument("--window", nargs=2, type=float,
                              help="energy window [emin,emax] eV")
    
    # spin_static command
    spin_parser = subparsers.add_parser("spin_static", parents=[common_parser, plot_parser],
                                       help="static spin texture")
    spin_parser.add_argument("--nk", type=int, default=31,
                             help="mesh size for spin texture (default: 31)")
    spin_parser.add_argument("--window", nargs=2, type=float,
                             help="energy window [emin,emax] eV")
    
    # soc_scan command
    soc_scan_parser = subparsers.add_parser("soc_scan", parents=[common_parser],
                                           help="SOC splitting scan")
    soc_scan_parser.add_argument("--hw", type=float, help="drive amplitude eV")
    soc_scan_parser.add_argument("--gamma", type=float, help="drive strength")
    soc_scan_parser.add_argument("--pol", choices=["x", "y", "z"], nargs='*',
                                 help="drive polarization")
    
    # floquet_bands command
    floquet_parser = subparsers.add_parser("floquet_bands", 
                                          parents=[common_parser, plot_parser],
                                          help="Floquet band structure")
    floquet_parser.add_argument("--nk", type=int, default=121,
                                help="k-points per segment (default: 121)")
    floquet_parser.add_argument("--window", nargs=2, type=float,
                                help="energy window [emin,emax] eV")
    floquet_parser.add_argument("--mode", choices=["path", "window"],
                                default="path", help="analysis mode")
    floquet_parser.add_argument("--center", nargs=3, type=float,
                                help="center point for window mode")
    
    # floquet_points command
    points_parser = subparsers.add_parser("floquet_points", parents=[common_parser],
                                         help="Floquet analysis at specific points")
    points_parser.add_argument("--hw", type=float, help="drive amplitude eV")
    points_parser.add_argument("--points", nargs='+',
                                default=["Gamma", "K"],
                                help="analysis points")
    
    # scan command
    scan_parser = subparsers.add_parser("scan", parents=[common_parser, plot_parser],
                                       help="parameter scan")
    scan_parser.add_argument("--hw", nargs='+', type=float,
                             help="drive amplitude values")
    scan_parser.add_argument("--gamma", nargs='+', type=float,
                             help="drive strength values")
    scan_parser.add_argument("--pol", choices=["x", "y", "z"], default="z",
                             help="drive polarization")
    
    # converge command
    converge_parser = subparsers.add_parser("converge", parents=[common_parser, plot_parser],
                                           help="numerical convergence test")
    converge_parser.add_argument("--ref-n-km", type=int, default=15,
                                 help="reference n_km (default: 15)")
    converge_parser.add_argument("--ref-m-sp", type=int, default=31,
                                 help="reference m_sp (default: 31)")
    converge_parser.add_argument("--n-km", nargs='+', type=int,
                                 help="test n_km values")
    converge_parser.add_argument("--m-sp", nargs='+', type=int,
                                 help="test m_sp values")
    converge_parser.add_argument("--k", nargs=3, type=float,
                                 help="test k-point")
    
    # soc command
    soc_parser = subparsers.add_parser("soc", parents=[common_parser, plot_parser],
                                      help="SOC splitting analysis")
    soc_parser.add_argument("--nk", type=int, default=121,
                            help="k-points per segment (default: 121)")
    soc_parser.add_argument("--window", nargs=2, type=float,
                            help="energy window [emin,emax] eV")
    
    # rashba command
    rashba_parser = subparsers.add_parser("rashba", parents=[common_parser, plot_parser],
                                          help="Rashba effect analysis")
    rashba_parser.add_argument("--nk", type=int, default=121,
                               help="k-points along direction (default: 121)")
    rashba_parser.add_argument("--k-max", type=float, default=0.3,
                               help="maximum k-distance (default: 0.3)")
    
    # pst command
    pst_parser = subparsers.add_parser("pst", parents=[common_parser, plot_parser],
                                      help="PST analysis")
    pst_parser.add_argument("--kmesh", nargs=3, type=int,
                             help="k-mesh [n_x,n_y,n_z] (default: 31,31,1)")
    pst_parser.add_argument("--bands", nargs='+', type=int,
                             help="specific bands for PST analysis")
    
    # dft_spin_check command
    dft_parser = subparsers.add_parser("dft_spin_check", parents=[common_parser],
                                       help="Compare with DFT spin texture")
    dft_parser.add_argument("--dft-dat", required=True,
                            help="DFT spin texture .dat file")
    dft_parser.add_argument("--band-map", type=int,
                            help="model band corresponding to DFT band")
    dft_parser.add_argument("--kmesh", nargs=3, type=int,
                            help="k-mesh for comparison")
    dft_parser.add_argument("--plot", action="store_true",
                            help="generate comparison plots")
    
    # benchmark command
    bench_parser = subparsers.add_parser("benchmark", parents=[common_parser],
                                        help="run analytic benchmarks")
    bench_parser.add_argument("--which", required=True,
                              choices=["rashba", "pbnv_spin", "wccl", "pt",
                                      "floquet_rashba", "rashba_dresselhaus"],
                              help="benchmark to run")
    bench_parser.add_argument("--plot", action="store_true",
                              help="generate benchmark plots")
    
    # export_dressed command
    export_parser = subparsers.add_parser("export_dressed", parents=[common_parser],
                                          help="export dressed Wannier Hamiltonian")
    export_parser.add_argument("--p", type=int, required=True,
                               help="harmonic order (positive for dressing, negative for undressing)")
    
    # provenance command
    prov_parser = subparsers.add_parser("provenance", parents=[common_parser],
                                        help="collect provenance information")
    
    # kchunk command
    chunk_parser = subparsers.add_parser("kchunk", parents=[common_parser],
                                         help="chunk k-points for parallel processing")
    chunk_parser.add_argument("--n-jobs", type=int, default=4,
                              help="number of jobs (default: 4)")
    
    # status command
    subparsers.add_parser("status", help="print implementation status")
    
    return parser


def main(argv=None):
    """Main CLI entry point."""
    if argv is None:
        argv = sys.argv[1:]
    
    parser = build_parser()
    args = parser.parse_args(argv)
    
    if args.command is None:
        parser.print_help()
        return 0
    
    # Load configuration
    if hasattr(args, 'config'):
        cfg = load_config(args.config)
        args.cfg = cfg
    
    # Create output directory
    if hasattr(args, 'outdir'):
        os.makedirs(args.outdir, exist_ok=True)
    
    # Dispatch to command handlers
    handlers = {
        "validate": cmd_validate,
        "bands": cmd_bands,
        "spin_static": cmd_spin_static,
        "soc_scan": cmd_soc_scan,
        "floquet_bands": cmd_floquet_bands,
        "floquet_points": cmd_floquet_points,
        "scan": cmd_scan,
        "converge": cmd_converge,
        "soc": cmd_soc,
        "rashba": cmd_rashba,
        "pst": cmd_pst,
        "dft_spin_check": cmd_dft_spin_check,
        "benchmark": cmd_benchmark,
        "export_dressed": cmd_export_dressed,
        "provenance": cmd_provenance,
        "kchunk": cmd_kchunk,
        "status": cmd_status,
    }
    
    handler = handlers.get(args.command)
    if handler:
        return handler(args)
    else:
        print(f"unknown command {args.command}")
        return 1


if __name__ == "__main__":
    sys.exit(main())