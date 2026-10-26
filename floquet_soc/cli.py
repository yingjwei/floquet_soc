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
                                         local_grid_at_hw, make_drive,
                                         save_and_plot as fp_save)
from .observables.floquet_bands import (compute_path as fb_compute,
                                        save_and_plot as fb_save)
from .observables.soc_scan import (run_scan as soc_scan_run,
                                   save as soc_scan_save)
from .observables.spin_plot import (plot_spin_texture_2d,
                                    plot_spin_texture_panels,
                                    plane_kpts)
from .workflows.scan import run_scan, rows_to_arrays, observable_map
from .workflows.parallel import run_chunks
from .workflows.converge import run_convergence
from .observables.floquet_spin import floquet_spin_observables

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Output root: explicit env var, otherwise ./runs in the current directory.
# Never derive it from __file__: that made the output land next to the install
# (site-packages/runs) and depended on where the package happened to live.
DEFAULT_OUT = os.environ.get("FLOQUET_SOC_OUT") or os.path.join(
    os.getcwd(), "runs")


# ------------------------------------------------------------------ config
# load_config / build_model live in core.model_factory (no core->cli
# dependency); re-exported here for the existing callers.
from .core.model_factory import load_config, build_model  # noqa: F401


def kpath_from_cfg(cfg, nk_seg=None):
    kp = cfg.get("kpath")
    if not kp:
        return None, None
    labels = [p[0] for p in kp]
    pts = np.array([p[1] for p in kp], float)
    nseg = cfg.get("nk_seg", 20) if nk_seg is None else nk_seg
    kpts = []
    for i in range(len(pts) - 1):
        seg = np.linspace(pts[i], pts[i + 1], nseg, endpoint=False)
        kpts.extend(seg.tolist())
    kpts.append(pts[-1].tolist())
    return np.array(kpts), labels


# ------------------------------------------------------------------ drive
def drive_from_cfg(cfg, model):
    """Build a Drive from the material config (a_ref from lattice for Wannier)."""
    drv = dict(cfg.get("drive") or {})
    if not drv:
        raise ValueError("config has no 'drive' block")
    a_ref = getattr(model, "a_ref", None) or float(drv.get("a_ref", 1.0))
    return Drive.from_config(drv, a_ref=a_ref)


def numerics_from_cfg(cfg, model):
    num = cfg.get("numerics", {})
    return dict(n_step=int(num.get("n_step", 150)),
                m_sp=int(num.get("m_sp", 16)))


# ------------------------------------------------------------------ validate
def _spin_status_wannier(cfg, model, rep):
    spin_cfg = cfg.get("spin", {})
    spn_path = spin_cfg.get("file")
    if spin_cfg.get("source") != "spn":
        rep.add("spin_operator", SKIP,
                "no spin source configured (only Hamiltonian checks run)")
        return
    if not spn_path or not os.path.exists(spn_path):
        rep.add("spin_operator", BLOCKED,
                f"spin source 'spn' requested but file missing: {spn_path}")
        return
    spn = read_wannier_spn(spn_path)
    metrics, checks = validate_spin_matrices(spn["matrices"])
    fmt = spn.get("format", "unknown")
    rep.add("spn_format", PASS,
            f"{os.path.basename(spn_path)}: {fmt}, header={spn['header']!r}, "
            f"num_bands={spn['num_bands']}, num_kpts={spn['num_kpts']}")
    rep.add("spn_hermiticity", PASS if metrics["hermiticity"] < 1e-10 else FAIL,
            f"max|S-S^dag| = {metrics['hermiticity']:.2e}")
    rep.add("spn_real_diagonal", PASS if metrics["diag_imag"] < 1e-10 else FAIL,
            f"max|Im diag| = {metrics['diag_imag']:.2e}")
    # quality: Pauli bound and spin sum rule
    bound_ok = metrics["eig_min"] >= -1.0 - 1e-3 and metrics["eig_max"] <= 1.0 + 1e-3
    rep.add("spn_pauli_eigenvalue_bound",
            PASS if bound_ok else WARN,
            f"eigenvalues in [{metrics['eig_min']:+.3f},"
            f"{metrics['eig_max']:+.3f}]; Pauli bound [-1,1]"
            + ("" if bound_ok else
               " -> violated: the .spn matrices are not a set of orthonormal "
               "spin matrices (pseudo-wavefunction error or wrong file "
               "layout); spin results are approximate"))
    s2 = (metrics["S2_min"], metrics["S2_max"])
    s2_ok = abs(s2[0] - 3.0) < 0.05 and abs(s2[1] - 3.0) < 0.05
    rep.add("spn_spin_sum_rule", PASS if s2_ok else WARN,
            f"eig(sum_i S_i^2) in [{s2[0]:.3f},{s2[1]:.3f}] "
            f"(3 for a spin-complete basis; <3 means spin weight outside the "
            f"Wannier band window)")
    diag_ok = metrics["diag_pol_max"] <= 1.0 + 1e-3
    rep.add("spn_band_polarization_bound", PASS if diag_ok else WARN,
            f"max per-band |<sigma>| = {metrics['diag_pol_max']:.3f}")

    if spn["num_bands"] != model.nw:
        try:
            spin = load_spin_for_model(cfg, model)
        except MissingTransformError as exc:
            rep.add("spin_operator_wannier_basis", BLOCKED, str(exc))
            rep.add("spin_texture_wannier", BLOCKED,
                    "static/Floquet spin texture for this material cannot be "
                    "computed gauge-consistently until the transformation "
                    "files are provided; the framework refuses to guess band "
                    "ordering")
            return
        except SpinInputError as exc:
            rep.add("spin_operator_wannier_basis", FAIL, str(exc))
            return
        quality = ("spin_wannier_eigenvalue_bound", "spin_wannier_sum_rule",
                   "spin_wannier_trs",
                   "spn_pauli_eigenvalue_bound", "spn_spin_sum_rule",
                   "spin_eigenvalue_bound", "spin_trace_zero",
                   "spin_diagonal_polarization_bound")
        for name, passed, detail in spin["checks"]:
            if name in quality and not passed:
                rep.add(name, WARN, detail)
            else:
                rep.add(name, PASS if passed else FAIL, detail)
        hc = spin["metadata"].get("h_consistency_ev")
        rep.add("spin_operator_wannier_basis", PASS,
                f"chk rotation applied ({spin['metadata']['transform']}); "
                f"max|H_W-H_hr| = {hc:.2e} eV"
                if hc is not None else "chk rotation applied")
        rep.add("spin_texture_wannier", PASS,
                "gauge-consistent spin operators available; static spin "
                "texture can be computed (see `spin_static`)")
        return
    rep.add("spin_operator_wannier_basis", PASS,
            "spn matrices have num_bands == num_wann; gauge check required")


def cmd_validate(args):
    cfg = load_config(args.config)
    name = cfg.get("name", os.path.basename(args.config))
    rep = ValidationReport(f"validate: {name}")
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "static")
    try:
        prov_write(prov_collect(args.config, cfg), outdir)
    except Exception:                                     # noqa: BLE001
        pass
    try:
        model = build_model(cfg)
    except (WannierInputError, ValueError) as exc:
        rep.add("model_parse", FAIL, str(exc))
        rep.print()
        rep.save(outdir, "validate")
        return 1

    rep.add("model_parse", PASS,
            f"kind={cfg.get('kind')}, num_wann={getattr(model, 'nw', None)}")
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


# ------------------------------------------------------------------ bands
def cmd_bands(args):
    cfg = load_config(args.config)
    model = build_model(cfg)
    kpts, labels = kpath_from_cfg(cfg, nk_seg=args.nk)
    if kpts is None:
        kpts = np.array([[k, 0.0, 0.0] for k in np.linspace(0, 0.3, 121)])
        labels = None
    e = model.static_bands(kpts)
    outdir = args.outdir or os.path.join(DEFAULT_OUT, cfg.get("name", "material"), "static")
    os.makedirs(outdir, exist_ok=True)
    npz = os.path.join(outdir, "bands.npz")
    np.savez(npz, kpts=kpts, energies=e,
             labels=np.array(labels if labels else []))
    print(f"bands: {len(kpts)} k-points x {e.shape[1]} bands -> {npz}")
    print(f"  energy range: [{e.min():.3f}, {e.max():.3f}] eV")
    ne = cfg.get("num_electrons")
    e_zero = float(e[:, int(ne) - 1].max()) if ne else 0.0
    ylabel = ("Energy $-$ $E_{\\rm VBM}$ (eV)" if ne else "Energy (eV)")
    _plot_bands(kpts, e, labels, os.path.join(outdir, "bands.png"),
                window=cfg.get("band_window"), lattice=_model_lattice(model),
                e_zero=e_zero, ylabel=ylabel, ne=ne)
    return 0


def _plot_bands(kpts, e, labels, path, window=None, lattice=None,
                e_zero=0.0, ylabel="Energy (eV)", ne=None):
    """Static band structure in the phonkit PRL style (v0.1.2-rc1+).

    Bands are drawn as continuous lines over the full k path (no NaN masking
    inside the path, which used to make the curves look broken); the visible
    energy window only sets the y limits.
    """
    import matplotlib.pyplot as plt
    kd = ps.kpath_distance(kpts, lattice)
    lo, hi = (window if window else (e.min(), e.max()))
    ylo, yhi = lo - e_zero, hi - e_zero
    pad = 0.25 * (yhi - ylo)
    with ps.apply_style():
        fig, ax = plt.subplots(figsize=ps.FIG_SINGLE, constrained_layout=True)
        pair_lo = pair_hi = None
        for ib in range(e.shape[1]):
            band = e[:, ib] - e_zero
            if not ((band >= ylo - pad) & (band <= yhi + pad)).any():
                continue
            ax.plot(kd, band, color=ps.COLOR_SINGLE, lw=0.6, alpha=0.75)
        if ne:
            pair_lo = [int(ne) - 2, int(ne) - 1]
            pair_hi = [int(ne), int(ne) + 1]
            for ib in pair_lo:
                if 0 <= ib < e.shape[1]:
                    ax.plot(kd, e[:, ib] - e_zero, color="#0072B2", lw=1.4,
                            label="valence pair" if ib == pair_lo[0] else None)
            for ib in pair_hi:
                if 0 <= ib < e.shape[1]:
                    ax.plot(kd, e[:, ib] - e_zero, color="#D55E00", lw=1.4,
                            label="conduction pair" if ib == pair_hi[0] else None)
        ps.add_reference_line(ax, 0.0)
        if labels:
            ps.add_sym_labels(ax, kd, labels)
        else:
            ax.set_xlim(kd[0], kd[-1])
            ax.set_xlabel("k distance (1/$\\AA$)")
        ax.set_ylabel(ylabel, fontsize=ps.FONTS["label"])
        ax.set_ylim(ylo, yhi)
        if ne:
            ax.legend(frameon=False, fontsize=ps.FONTS["legend"],
                      loc="best")
        ps.save_figure(fig, path)
    print(f"figure: {path}")


def spin_plot_band(cfg):
    """Band index used in spin-texture plots (default: VBM)."""
    ne = cfg.get("num_electrons")
    if ne:
        return int(cfg.get("spin_plot_band", ne - 1))
    return int(cfg.get("spin_plot_band", 0))


# ------------------------------------------------------------------ spin_static
def cmd_spin_static(args):
    cfg = load_config(args.config)
    name = cfg.get("name", "material")
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "static")
    os.makedirs(outdir, exist_ok=True)
    n_km, kmax_km = resolve_kmesh(cfg, "kmesh", default_n=21,
                                  default_kmax=0.55)
    if getattr(args, "kmesh_n", None):
        n_km = int(args.kmesh_n)
    if getattr(args, "kmax", None) is not None:
        kmax_km = float(args.kmax)
    kk = plane_kpts(cfg.get("kplane"), n_km, kmax_km)
    jobs = _jobs_for(args, cfg, len(kk))
    if jobs > 1 and _parallel_ok(cfg):
        lat = _lattice_from_cfg(cfg)
        print(f"spin_static: {len(kk)} k-points (n={n_km}, "
              f"|k_i| <= {kmax_km}), jobs={jobs}")
        load_model_spin_cached(cfg)   # warm cache for workers
        tex = run_chunks("spin_static", os.path.abspath(args.config), kk, jobs)
        npz = os.path.join(outdir, "spin_static.npz")
        np.savez(npz, kpts=tex["kpts"], energies=tex["energies"],
                 spins=tex["spins"])
        print(f"spin texture: {len(kk)} k-points -> {npz}")
        ib = min(spin_plot_band(cfg), tex["spins"].shape[1] - 1)
        _plot_spin_texture(tex["kpts"], tex["spins"],
                           os.path.join(outdir, "spin_static.png"),
                           title=f"{name}: static spin texture, band {ib}",
                           iband=ib, lattice=lat, kplane=cfg.get("kplane"))
        return 0
    try:
        model, spin = load_model_spin_cached(cfg)
    except MissingTransformError as exc:
        print("BLOCKED: spin texture cannot be computed gauge-consistently.")
        print(str(exc))
        print("missing files:", ", ".join(exc.missing_files))
        return 2
    except NotImplementedError as exc:
        print("BLOCKED:", exc)
        return 2
    spin_at_k = spin.get("spin_at_k")
    lat = _model_lattice(model)
    kplane = cfg.get("kplane")
    if spin_at_k is not None:
        op = spin.get("spin_operator")
        if op is not None and hasattr(op, "S_R"):
            tex = spin_texture_wannier_batch(model, op, kk, chunk=64)
        else:
            nw = model.nw
            energies = np.empty((len(kk), nw))
            spins = np.empty((len(kk), nw, 3))
            for ik, k in enumerate(kk):
                ops_k = spin_at_k(k)
                Hk = model.Hk(k)
                r = spin_texture(Hk, ops_k)
                energies[ik] = r["energies"]
                spins[ik] = r["spins"]
            tex = dict(energies=energies, spins=spins)
    else:
        ops = spin["ops"]
        if ops is None:
            print("BLOCKED: spin operators are not available in the Wannier "
                  "basis.")
            return 2
        tex = spin_texture_path(model, ops, kk)
    npz = os.path.join(outdir, "spin_static.npz")
    np.savez(npz, kpts=kk, energies=tex["energies"], spins=tex["spins"])
    print(f"spin texture: {len(kk)} k-points -> {npz}")
    ib = min(spin_plot_band(cfg), tex["spins"].shape[1] - 1)
    _plot_spin_texture(kk, tex["spins"], os.path.join(outdir, "spin_static.png"),
                       title=f"{name}: static spin texture, band {ib}",
                       iband=ib, lattice=lat, kplane=kplane)
    return 0


def _model_lattice(model):
    """2D-lattice rows for BZ plotting (None for continuum/analytic models)."""
    lat = getattr(model, "lattice", None)
    return np.asarray(lat, float) if lat is not None else None


def _lattice_from_cfg(cfg):
    """Lattice rows without building the model (cheap read of win/chk)."""
    files = cfg.get("files") or {}
    win = files.get("win")
    if win and os.path.exists(win):
        try:
            return np.asarray(read_win(win)["lattice"], float)
        except Exception:                                 # noqa: BLE001
            pass
    return None


def _nw_from_cfg(cfg):
    files = cfg.get("files") or {}
    win = files.get("win")
    if win and os.path.exists(win):
        return int(read_win(win)["num_wann"])
    hr = files.get("hr")
    if hr and os.path.exists(hr):
        return int(read_hr(hr)["num_wann"])
    return int(cfg.get("num_wann") or 0)


def resolve_kmesh(cfg, key, default_n=21, default_kmax=0.55):
    """Resolve a k-mesh config, supporting physical spacing ``dk``.

    ``dk`` is in units of 2*pi/Angstrom (same convention as the VASP KPOINTS
    generator: K-spacing 0.01 -> n_i = round(BZ_span_i/(0.01*2*pi))).
    Returns (n_ints, kmax_values) where each may be a tuple (x, y).
    """
    km = dict(cfg.get(key) or cfg.get("kmesh") or {})
    if km.get("dk") and cfg.get("kind") == "wannier90":
        lat = _lattice_from_cfg(cfg)
        if lat is not None:
            from .observables.spin_plot import (bz_polygon_plane,
                                                plane_matrices)
            kplane = cfg.get("kplane")
            poly = bz_polygon_plane(lat, kplane)
            if len(poly) >= 3:
                E, M = plane_matrices(lat, kplane)
                cart = poly @ M
                span = np.ptp(cart[:, :2], axis=0)
                step = float(km["dk"]) * 2.0 * np.pi
                n1 = max(3, int(round(span[0] / step)))
                n2 = max(3, int(round(span[1] / step)))
                kmax = (float(np.max(np.abs(poly[:, 0]))) * 1.02,
                        float(np.max(np.abs(poly[:, 1]))) * 1.02)
                return (n1, n2), kmax
    return (int(km.get("n", default_n)),
            float(km.get("kmax", default_kmax)))


def _jobs_for(args, cfg, nk, min_points=16):
    env = os.environ.get("FLOQUET_SOC_JOBS")
    jobs = getattr(args, "jobs", None)
    if jobs is None and env:
        jobs = env                      # validated by int(jobs) below
    if jobs is None:
        jobs = cfg.get("jobs", 1)
    try:
        jobs = int(jobs)
    except (TypeError, ValueError):
        jobs = 1
    if jobs <= 1 or nk < min_points:
        return 1
    return max(1, min(jobs, nk))


def _parallel_ok(cfg, model=None):
    if cfg.get("kind") != "wannier90":
        return False
    if cfg.get("spin", {}).get("source") != "spn":
        return False
    files = cfg.get("files") or {}
    spn = cfg.get("spin", {}).get("file")
    chk = cfg.get("spin", {}).get("chk") or files.get("chk")
    eig = cfg.get("spin", {}).get("eig") or files.get("eig")
    return all(x and os.path.exists(x) for x in (spn, chk, eig))


def _load_spin_or_exit(cfg, model):
    """``load_spin_for_model`` with the framework-wide BLOCKED exit.

    Returns ``None`` after printing a BLOCKED diagnostic when the spin
    operators cannot be built (missing ``.chk``/``.eig`` transform, unknown
    source, ...).  Callers do ``if spin is None: return 2``.  Commands that
    used to call ``load_spin_for_model`` directly leaked a traceback on these
    inputs -- the same condition that ``spin_static``/``pst`` report cleanly.
    """
    try:
        return load_spin_for_model(cfg, model)
    except MissingTransformError as exc:
        print("BLOCKED: spin texture cannot be computed gauge-consistently.")
        print(str(exc))
        print("missing files:", ", ".join(exc.missing_files))
    except SpinInputError as exc:
        print("BLOCKED: spin operators cannot be built for this material.")
        print(str(exc))
    return None


def _plot_spin_texture(kpts, spins, path, title="", iband=0, lattice=None,
                       kplane=None):
    plot_spin_texture_2d(kpts, np.asarray(spins)[:, iband, :], path,
                         title=title, lattice=lattice, kplane=kplane)
    print(f"figure: {path}")


# ------------------------------------------------------------------ soc scan
def cmd_soc_scan(args):
    cfg = load_config(args.config)
    name = cfg.get("name", "material")
    model = build_model(cfg)
    spin = _load_spin_or_exit(cfg, model)
    if spin is None:
        return 2
    spin_at = spin.get("spin_at_k")
    if spin_at is None:
        print("BLOCKED: no spin operator available.")
        return 2
    drv = cfg.get("drive") or {}
    hw = float(args.hw if args.hw is not None else drv.get("hw", 0.4))
    gammas = args.gamma if args.gamma else [0.0, 0.02, 0.05]
    pols = args.pol if args.pol else ["circ+", "circ-", "linear"]
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "floquet",
                                         "soc_scan")
    print(f"soc_scan: hw={hw:g}, gamma={gammas}, pol={pols}")
    res = soc_scan_run(model, spin_at, cfg, hw, gammas, pols,
                       points=args.points, n_step=args.n_step,
                       m_sp=args.m_sp, pst_n=args.pst_n,
                       pst_kmax=args.pst_kmax)
    out = soc_scan_save(res, outdir)
    for pol, rows in res["alpha"].items():
        print(f"  alpha_R[{pol}]:",
              {g: round(v["alpha_eff_ev_a"] * 1000, 2)
               for g, v in rows.items()})
    for r in res["pst"]:
        print(f"  PST {r['pol']:7s} g={r['gamma']:.3g} "
              f"coll={r['collinearity']:.3f} rot={r['axis_rotation_deg']:.1f} deg")
    print("  ->", out["json"])
    return 0


# ------------------------------------------------------------------ floquet bands
def cmd_floquet_bands(args):
    cfg = load_config(args.config)
    name = cfg.get("name", "material")
    model = build_model(cfg)
    kpts, labels = kpath_from_cfg(cfg, nk_seg=args.nk_seg)
    if kpts is None:
        print("ERROR: config has no kpath; provide kpath in the JSON.")
        return 2
    drv = cfg.get("drive") or {}
    hw = float(args.hw if args.hw is not None else drv.get("hw", 0.4))
    gamma = float(args.gamma if args.gamma is not None
                  else drv.get("gamma", 0.02))
    pol = args.pol or ("linear" if drv.get("polarization") == "linear"
                       else "circ+")
    drive = make_drive(hw, gamma, pol, getattr(model, "a_ref", 1.0))
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "floquet")
    print(f"floquet_bands: {len(kpts)} k-points on {labels}, "
          f"hw={hw:g}, gamma={gamma:g}, pol={pol}")
    res = fb_compute(model, drive, kpts, n_step=args.n_step,
                     m_sp=args.m_sp, ne=cfg.get("num_electrons"))
    res["ne"] = cfg.get("num_electrons")
    out = fb_save(res, labels, outdir)
    print("  ->", out["figure"])
    return 0


# ------------------------------------------------------------------ floquet points
def cmd_floquet_points(args):
    """Band-resolved Floquet spin scan over selected k points."""
    cfg = load_config(args.config)
    name = cfg.get("name", "material")
    model = build_model(cfg)
    spin = _load_spin_or_exit(cfg, model)
    if spin is None:
        return 2
    spin_at = spin.get("spin_at_k")
    if spin_at is None:
        print("BLOCKED: no spin operator available.")
        return 2
    points = resolve_points(cfg, args.points)
    hw = float(args.hw if args.hw is not None
               else (cfg.get("drive") or {}).get("hw", 0.4))
    gammas = args.gamma if args.gamma else [0.0, 0.02, 0.05]
    pols = args.pol if args.pol else ["circ+", "circ-", "linear"]
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "floquet")
    print(f"floquet_points: hw={hw:g} eV, points={[p[0] for p in points]}, "
          f"gamma={gammas}, pol={pols}")
    res = run_points(model, spin_at, cfg, points, hw, gammas, pols,
                     n_step=args.n_step, m_sp=args.m_sp,
                     subspace=args.subspace,
                     subspace_tol=args.subspace_tol)
    asym = valley_asymmetry(res)
    extra = {}
    if args.rashba:
        kp = [k for lab, k in points if lab in ("K", "K′")]
        direction = kp[0] if kp else np.array([1 / 3, 1 / 3, 0.0])
        extra["rashba"] = rashba_effective(
            model, spin_at, cfg, hw, gammas, direction_frac=direction,
            tmax=args.rashba_tmax, npts=args.rashba_n,
            n_step=args.n_step, m_sp=args.m_sp)
    if args.grid:
        grid = {}
        for lab, k in resolve_points(cfg, args.grid_center):
            grid[lab] = local_grid_at_hw(
                model, spin_at, cfg, k, args.grid_half, args.grid_n, hw,
                args.grid_gamma, args.grid_pol,
                n_step=args.n_step, m_sp=args.m_sp)
        extra["grid"] = grid
    if extra:
        res["extra"] = extra
    res["provenance"] = prov_collect(args.config, cfg)
    out = fp_save(res, asym, outdir)
    names = [res["band_names"][b] for b in res["bands"]]
    for lab in res["points"]:
        print(f"  {lab}:")
        for key, entry in res["points"][lab]["cases"].items():
            rel = [nm for nm in names
                   if entry["targets"][nm].get("reliable", True)]
            dS = max(entry["targets"][nm]["dS_norm"]
                     for nm in rel) if rel else float("nan")
            print(f"    {key:16s} max|dS|(reliable)={dS:.3e}  "
                  f"n_unreliable={len(names)-len(rel)}  "
                  f"w_main_min={min(entry['targets'][nm]['w_main'] for nm in names):.3f}")
    for key, rows in asym.items():
        worst = max(max(rows[k][nm]["dEbar"] for nm in names) for k in rows)
        print(f"  valley asymmetry {key}: max|dEbar|={worst*1000:.3f} meV")
    print("  ->", out["json"])
    return 0


# ------------------------------------------------------------------ dressed export
def cmd_export_dressed(args):
    """Write the light-dressed Wannier Hamiltonian H_p(R) as _hr.dat files.

    Each file is a standard Wannier90 ``_hr.dat`` whose Fourier sum is the
    photon component H_p(k) of the driven Hamiltonian; a downstream code can
    read them and build the Floquet/Sambe problem (or take the p-sum at any
    time) without reimplementing the Peierls factor.
    """
    cfg = load_config(args.config)
    if cfg.get("kind") != "wannier90":
        print("ERROR: export_dressed needs a wannier90 config "
              "(kind='wannier90').")
        return 2
    model = build_model(cfg)
    drive = drive_from_cfg(cfg, model)
    num = cfg.get("numerics") or {}
    pmax = args.pmax if args.pmax is not None else int(num.get("pmax", 3))
    hr = peierls_hr_components(model, drive, pmax=pmax)
    name = cfg.get("name", "material")
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "dressed")
    os.makedirs(outdir, exist_ok=True)
    nw = int(model.nw)
    nd = hr["ndegen"]
    if nd is None:                    # expanded wsvec representation
        nd = np.ones(len(hr["Rlist"]), int)
    files = []
    for p in sorted(hr["Hp_R"]):
        fname = f"wannier90_hr_p{p:+d}.dat"
        write_hr(os.path.join(outdir, fname), nw, hr["Rlist"], nd,
                 hr["Hp_R"][p])
        files.append(fname)
    meta = dict(name=name, hw=float(drive.hw), gamma=float(drive.gamma),
                polarization=drive.polarization, sign=int(drive.sign),
                pmax=int(pmax), gmax=float(hr["gmax"]),
                bessel_truncation_residual=float(hr["truncation_residual"]),
                num_wann=nw, nrpts=int(len(hr["Rlist"])),
                wsvec_expanded=bool(getattr(model, "HR_full", None) is not None),
                convention=("H_p(k) = sum_R H_p(R)/ndegen exp(2*pi*i k.R); "
                            "H(k,t) = sum_p H_p(k) exp(-i p omega t)"),
                files=files)
    meta_path = os.path.join(outdir, "dressed_meta.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    print(f"export_dressed: {nw} Wannier, {len(hr['Rlist'])} R points, "
          f"pmax={pmax}, g_max={hr['gmax']:.4f}, "
          f"Bessel tail={hr['truncation_residual']:.2e}")
    for p in sorted(hr["Hp_R"]):
        print(f"  p={p:+d}  -> "
              f"{os.path.join(outdir, f'wannier90_hr_p{p:+d}.dat')}")
    print(f"  meta  -> {meta_path}")
    return 0


# ------------------------------------------------------------------ provenance
def cmd_provenance(args):
    cfg = load_config(args.config)
    name = cfg.get("name", "material")
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "static")
    prov = prov_collect(args.config, cfg)
    path = prov_write(prov, outdir)
    print(f"provenance: {prov['version']}  config sha256 {prov['config']['sha256'][:12]}")
    for key, e in prov["files"].items():
        print(f"  {key:8s} {'OK ' if e['exists'] else 'MISS'} "
              f"{e['path']}" + (f"  sha256 {e['sha256'][:12]}" if e["exists"] else ""))
    print("  ->", path)
    return 0


# ------------------------------------------------------------------ dft spin check
def cmd_dft_spin_check(args):
    """Compare a model spin texture with a DFT spin-texture table."""
    cfg = load_config(args.config)
    name = cfg.get("name", "material")
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "static")
    dft = read_spin_texture_dat(args.dft)
    try:
        model, spin = load_model_spin_cached(cfg)
    except MissingTransformError as exc:
        print("BLOCKED: spin operator cannot be built gauge-consistently.")
        print(str(exc))
        return 2
    spin_at = spin.get("spin_at_k")
    if spin_at is None:
        print("BLOCKED: no spin operator available.")
        return 2
    band = resolve_band(cfg, model, args.band) if args.band else None
    res = dft_compare(model, spin_at, dft, band=band,
                      mask_split=args.mask_split)
    res["provenance"] = prov_collect(args.config, cfg)
    out = dft_save_plot(res, model, spin_at, outdir, prefix=args.prefix)
    print(f"dft-spin-check: DFT band {dft['band_index']} <-> model band "
          f"{res['band_wannier']}, {res['n_points']} k-points")
    print(f"  energy: offset {res['energy']['offset_ev']*1000:.2f} meV, "
          f"residual RMS {res['energy']['rms_resid_ev']*1000:.2f} meV")
    print(f"  spin direction: normalized RMS "
          f"{res['direction']['normalized_rms_all']:.4f} (all), "
          f"{res['direction']['normalized_rms_non_degenerate']:.4f} "
          f"(non-degenerate)")
    print(f"  |S| mean: model {res['scale']['norm_w']:.3f}, DFT "
          f"{res['scale']['norm_dft']:.3f} (global scale "
          f"{res['scale']['global_fit']:.3f})")
    print("  ->", out["json"])
    return 0


# ------------------------------------------------------------------ spin_texture
def cmd_spin_texture(args):
    """Generic static spin-texture map (full BZ or a k window)."""
    cfg = load_config(args.config)
    name = cfg.get("name", "material")
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "static")
    os.makedirs(outdir, exist_ok=True)
    try:
        model, spin = load_model_spin_cached(cfg)
    except MissingTransformError as exc:
        print("BLOCKED: spin texture cannot be computed gauge-consistently.")
        print(str(exc))
        print("missing files:", ", ".join(exc.missing_files))
        return 2
    except NotImplementedError as exc:
        print("BLOCKED:", exc)
        return 2
    spin_at = spin.get("spin_at_k")
    if spin_at is None:
        print("BLOCKED: no spin operator available for this model.")
        return 2
    band = resolve_band(cfg, model, args.band)
    lat = _model_lattice(model)
    mode = args.mode
    if mode == "window":
        if args.center is None:
            print("ERROR: --mode window needs --center CX CY CZ "
                  "(fractional coordinates, e.g. K = 1/3 1/3 0).")
            return 2
        center = np.asarray(args.center, float)
        kpts, xy = window_kpts(lat, center, args.half, args.n)
        lattice = None
        kplane = None
        tag = "window"
    else:
        kpts = full_kpts(cfg, args.n, args.kmax)
        xy = kpts
        lattice = lat
        kplane = cfg.get("kplane")
        tag = "full"
    ratio = (args.inplane_ratio if args.inplane_ratio is not None
             else float(cfg.get("spin_inplane_ratio", 0.05)))
    tex = compute_map(model, spin_at, band, kpts,
                      mask_split=args.mask_split,
                      inplane_ratio=ratio, chunk=args.chunk)
    fig = os.path.join(outdir, f"spin_texture_{tag}.png")
    npz = os.path.join(outdir, f"spin_texture_{tag}.npz")
    np.savez(npz, kpts=tex["kpts"], energies=tex["energies"],
             spins_band=tex["spins_band"], split=tex["split"],
             band=tex["band"])
    plot_map(xy, tex["spins_band"], fig, lattice=lattice, kplane=kplane,
             figsize=tuple(args.figsize) if args.figsize else (3.2, 3.0),
             dot_ratio=ratio)
    meta = dict(mode=tag, band=tex["band"], n=int(args.n),
                kmax=float(args.kmax), half=float(args.half),
                center=(center.tolist() if mode == "window" else None),
                mask_split=float(args.mask_split),
                inplane_ratio=float(ratio),
                n_out_of_plane=int(tex["n_out_of_plane"]),
                n_points=int(len(kpts)), npz=npz, figure=fig)
    with open(os.path.join(outdir, f"spin_texture_{tag}.json"), "w",
              encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)
    print(f"spin_texture[{tag}]: band {band}, {len(kpts)} k-points "
          f"-> {npz}")
    return 0


# ------------------------------------------------------------------ stage C helpers
def _soc_window(cfg, model):
    """Energy window for SOC-splitting analysis (config or around the gap)."""
    soc = cfg.get("soc", {})
    if soc.get("emin") is not None and soc.get("emax") is not None:
        return float(soc["emin"]), float(soc["emax"])
    ne = cfg.get("num_electrons")
    if ne:
        kpts, _ = kpath_from_cfg(cfg, nk_seg=10)
        if kpts is None:
            kpts = np.array([[k, 0.0, 0.0] for k in np.linspace(0, 0.3, 21)])
        E = model.static_bands(kpts)
        vbm = float(E[:, ne - 1].max())
        cbm = float(E[:, ne].min())
        if not (np.isfinite(vbm) and np.isfinite(cbm)) or cbm <= vbm:
            # metal / semimetal: use a window around the Fermi level
            ef = 0.5 * (vbm + cbm)
            return ef - 1.5, ef + 1.5
        return vbm - 0.30, cbm + 0.30
    return -1.0, 1.0


def _pst_bands(cfg, model):
    ne = cfg.get("num_electrons")
    pst = cfg.get("pst", {})
    if pst.get("bands"):
        return [int(b) for b in pst["bands"]]
    if ne:
        return [int(x) for x in range(max(0, ne - 2), min(model.nw, ne + 2))]
    return list(range(min(4, model.nw)))


# ------------------------------------------------------------------ spectrum
def cmd_spectrum(args):
    cfg = load_config(args.config)
    model = build_model(cfg)
    drive = drive_from_cfg(cfg, model)
    kpts, labels = kpath_from_cfg(cfg, nk_seg=args.nk)
    if kpts is None:
        kmax = args.kmax
        kpts = np.array([[k, 0.0, 0.0] for k in np.linspace(0, kmax, 41)])
        labels = None
    spec_cfg = cfg.get("spectrum", {})
    # Wannier90 eigenvalues carry an arbitrary energy zero (Pt sits at
    # +7.6..+40.8 eV, PbNV at -64.8..+6.6 eV), so a hard-coded [-1.5, 1.5] eV
    # window can miss the whole band manifold and leave A(k,E) ~ 1e-29.
    # Use the same reference as bands / floquet_bands: E - E_ref.
    E_static = model.static_bands(kpts)
    ne = cfg.get("num_electrons")
    ref = (float(E_static[:, int(ne) - 1].max()) if ne
           else float(E_static.max()))
    emin = (args.emin if args.emin is not None
            else ref + float(spec_cfg.get("emin", -1.5)))
    emax = (args.emax if args.emax is not None
            else ref + float(spec_cfg.get("emax", 1.5)))
    nE = int(spec_cfg.get("nE", 1501))
    sigma = float(spec_cfg.get("sigma", 4e-3))
    num = numerics_from_cfg(cfg, model)
    n_step = args.n_step or num["n_step"]
    # Photon replicas sit at E = eps - m*hw with eps in [-hw/2, hw/2], so the
    # order needed to reach the plotted window is |m| ~ |E|/hw.  With the
    # default m_sp=16 a material whose bands sit at +19 eV (Pt) has *no* weight
    # inside the window -> A(k,E) comes out as ~1e-29.
    m_needed = int(np.ceil((max(abs(emin), abs(emax)) + drive.hw)
                           / drive.hw)) + 2
    m_sp = max(int(num["m_sp"]), m_needed)
    E_grid = np.linspace(emin, emax, nE)
    name = cfg.get("name", "material")
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "floquet")
    os.makedirs(outdir, exist_ok=True)
    print(f"spectrum: {model.nw} bands, {len(kpts)} k-points, "
          f"hw={drive.hw:g}, gamma={drive.gamma:g}, n_step={n_step}")
    data, spec = [], []
    for ik, k in enumerate(kpts):
        fl = floquet_modes(model, drive, k, n_step=n_step, want_wspec=True,
                           m_sp=m_sp)
        data.append(fl)
        spec.append(spectral_function(fl["eps_u"], fl["wspec"]["ms"],
                                      fl["wspec"]["wm"], E_grid, sigma,
                                      drive.hw))
        if (ik + 1) % max(1, len(kpts) // 5) == 0:
            print(f"  k {ik + 1}/{len(kpts)}")
    spec = np.array(spec)
    npz = os.path.join(outdir, "spectrum.npz")
    np.savez(npz, kpts=kpts, E_grid=E_grid, spec=spec, ref=ref,
             labels=np.array(labels if labels else []),
             eps=np.array([np.sort(f["eps"]) for f in data]),
             drive=np.array([json.dumps(drive.as_dict())]))
    print(f"  -> {npz}")
    _plot_spectrum(kpts, E_grid, spec, labels, ref,
                   os.path.join(outdir, "spectrum.png"))
    return 0


def _plot_spectrum(kpts, E_grid, spec, labels, ref, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ps.enforce_style()
    fig, ax = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)
    im = ax.imshow(spec.T, origin="lower", aspect="auto",
                   extent=[-0.5, len(kpts) - 0.5,
                           E_grid[0] - ref, E_grid[-1] - ref],
                   cmap="inferno", vmin=0,
                   vmax=np.percentile(spec, 99.5),
                   interpolation="nearest")
    if labels:
        nseg = (len(kpts) - 1) // (len(labels) - 1)
        ax.set_xticks([i * nseg for i in range(len(labels))])
        ax.set_xticklabels(labels)
    else:
        ax.set_xlabel("k index", fontsize=ps.FONTS["label"])
    ax.set_ylabel(r"$E - E_{\rm ref}$ (eV)", fontsize=ps.FONTS["label"])
    ps.add_colorbar(fig, im, ax, label="A(k,E) [1/eV]")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  figure: {path}")


# ------------------------------------------------------------------ spin_floquet
def cmd_spin_floquet(args):
    cfg = load_config(args.config)
    name = cfg.get("name", "material")
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "floquet")
    os.makedirs(outdir, exist_ok=True)
    n_km, kmax_km = resolve_kmesh(cfg, "kmesh_spin_floquet", default_n=11,
                                  default_kmax=0.55)
    if args.kmesh_n:
        n_km = int(args.kmesh_n)
    if args.kmax is not None:
        kmax_km = float(args.kmax)
    kpts = plane_kpts(cfg.get("kplane"), n_km, kmax_km)
    jobs = _jobs_for(args, cfg, len(kpts))
    if jobs > 1 and _parallel_ok(cfg):
        lat = _lattice_from_cfg(cfg)
        drv = cfg.get("drive") or {}
        print(f"spin_floquet: {len(kpts)} k-points (n={n_km}), "
              f"hw={drv.get('hw', 0.4):g}, gamma={drv.get('gamma', 0.1):g}, "
              f"jobs={jobs}")
        load_model_spin_cached(cfg)   # warm cache for workers
        tex = run_chunks("spin_floquet", os.path.abspath(args.config), kpts,
                         jobs)
        npz = os.path.join(outdir, "spin_floquet.npz")
        np.savez(npz, kpts=tex["kpts"], S_bar=tex["S_bar"], S_t0=tex["S_t0"],
                 S_static=tex["S_static"],
                 S_bar_matched=tex["S_bar_matched"], purity=tex["purity"],
                 eps=tex["eps"],
                 drive=np.array([json.dumps(drv)]))
        print(f"  -> {npz}")
        ib = min(spin_plot_band(cfg), tex["S_bar_matched"].shape[1] - 1)
        _plot_spin_floquet(tex["kpts"], tex["S_bar_matched"],
                           tex["S_static"],
                           os.path.join(outdir, "spin_floquet.png"),
                           iband=ib, lattice=lat,
                           kplane=cfg.get("kplane"))
        pur = tex["purity"]
        print(f"  overlap matching purity: median={np.median(pur):.3f}, "
              f"fraction>0.8={np.mean(pur > 0.8):.3f}")
        return 0
    try:
        model, spin = load_model_spin_cached(cfg)
    except MissingTransformError as exc:
        print("BLOCKED: Floquet spin texture requires a gauge-consistent spin "
              "operator.")
        print(str(exc))
        print("missing files:", ", ".join(exc.missing_files))
        return 2
    except NotImplementedError as exc:
        print("BLOCKED:", exc)
        return 2
    ops = spin.get("ops")
    spin_at_k = spin.get("spin_at_k")
    if ops is None and spin_at_k is None:
        print("BLOCKED: spin operators are not available in the Wannier basis.")
        return 2
    drive = drive_from_cfg(cfg, model)
    num = numerics_from_cfg(cfg, model)
    lat = _model_lattice(model)
    kplane = cfg.get("kplane")
    print(f"spin_floquet: {len(kpts)} k-points, hw={drive.hw:g}, "
          f"gamma={drive.gamma:g}, n_step={num['n_step']}")
    if spin_at_k is not None:
        print("  using chk-rotated, k-interpolated spin operators")
    Sb, St, Sb_m, purity, eps = [], [], [], [], []
    for ik, k in enumerate(kpts):
        ops_k = spin_at_k(k) if spin_at_k is not None else ops
        obs = floquet_spin_observables(model, drive, k, ops_k,
                                       n_step=num["n_step"],
                                       m_sp=num["m_sp"])
        # match Floquet states to static bands by max |<psi_static|u(0)>|^2
        ev_s, V = np.linalg.eigh(model.Hk(k))
        ov = np.abs(V.conj().T @ obs["evec"]) ** 2       # (nstatic, nfloq)
        best = np.argmax(ov, axis=1)
        Sb.append(obs["S_bar"])
        St.append(obs["S_t0"])
        Sb_m.append(obs["S_bar"][best])
        purity.append(ov[np.arange(len(best)), best])
        eps.append(np.sort(obs["eps"]))
        if (ik + 1) % max(1, len(kpts) // 4) == 0:
            print(f"  k {ik + 1}/{len(kpts)}")
    Sb, St, eps = np.array(Sb), np.array(St), np.array(eps)
    Sb_m = np.array(Sb_m)
    purity = np.array(purity)
    stat = np.array([spin_texture(model.Hk(k),
                                  spin_at_k(k) if spin_at_k is not None else ops
                                  )["spins"] for k in kpts])
    npz = os.path.join(outdir, "spin_floquet.npz")
    np.savez(npz, kpts=kpts, S_bar=Sb, S_t0=St, S_static=stat,
             S_bar_matched=Sb_m, purity=purity, eps=eps,
             drive=np.array([json.dumps(drive.as_dict())]))
    print(f"  -> {npz}")
    ib = min(spin_plot_band(cfg), Sb.shape[1] - 1)
    _plot_spin_floquet(kpts, Sb_m, stat,
                       os.path.join(outdir, "spin_floquet.png"),
                       iband=ib, lattice=lat, kplane=kplane)
    print(f"  overlap matching purity: median={np.median(purity):.3f}, "
          f"fraction>0.8={np.mean(purity > 0.8):.3f}")
    return 0


def _plot_spin_floquet(kpts, Sb, stat, path, iband=0,
                       lattice=None, kplane=None):
    panels = [
        dict(kpts=kpts, spins=stat[:, iband, :], lattice=lattice,
             kplane=kplane,
             label=f"(a) static spin texture (band {iband})",
             colorbar_label=r"$S_z$"),
        dict(kpts=kpts, spins=Sb[:, iband, :], lattice=lattice,
             kplane=kplane,
             label=f"(b) Floquet periodic-average spin (band {iband})",
             colorbar_label=r"$S_z$"),
    ]
    d = np.linalg.norm(Sb[:, iband, :] - stat[:, iband, :], axis=1)
    panels.append(dict(kpts=kpts,
                       spins=np.stack([np.zeros_like(d), np.zeros_like(d), d],
                                      axis=1),
                       lattice=lattice, kplane=kplane, cmap="bwr",
                       vmin=0.0,
                       vmax=max(float(d.max()), 1e-6),
                       label="(c) drive-induced |S_bar - S_static|",
                       colorbar_label=r"$|\Delta S|$"))
    plot_spin_texture_panels(panels, path)
    print(f"  figure: {path}")


# ------------------------------------------------------------------ soc
def cmd_soc(args):
    cfg = load_config(args.config)
    model = build_model(cfg)
    kpts, labels = kpath_from_cfg(cfg, nk_seg=args.nk)
    if kpts is None:
        kpts = np.array([[k, 0.0, 0.0] for k in np.linspace(0, args.kmax, 41)])
        labels = None
    E = model.static_bands(kpts)
    emin, emax = _soc_window(cfg, model)
    system = cfg.get("system") or {}
    stype = system.get("type", "generic")
    if stype == "pt_nonmagnetic":
        from .observables.kramers import (kramers_delta, kramers_summary,
                                          inter_doublet_gaps)
        outdir = args.outdir or os.path.join(
            DEFAULT_OUT, cfg.get("name", "material"), "static")
        os.makedirs(outdir, exist_ok=True)
        d = kramers_delta(E)
        gaps = inter_doublet_gaps(E)
        summ = kramers_summary(E, kpts)
        print(f"soc: PT-symmetric non-magnetic system -> SOC splitting is "
              f"NOT APPLICABLE. Reporting Kramers residual and "
              f"inter-doublet gaps only.")
        print(f"  delta_PT(k)=max_n|E_{{2n-1}}-E_{{2n}}|: max "
              f"{summ['delta_max']*1000:.3f} meV, median "
              f"{summ['delta_median']*1000:.3f} meV, worst k="
              f"{np.round(summ['worst_k'], 4).tolist()}")
        if gaps.size:
            print(f"  inter-doublet gaps: max {gaps.max()*1000:.1f} meV, "
                  f"median {np.median(gaps)*1000:.1f} meV "
                  f"(these are NOT SOC spin splittings)")
        npz = os.path.join(outdir, "soc_splitting.npz")
        np.savez(npz, kpts=kpts, energies=E, split=adjacent_splitting(E),
                 kramers_delta=d, inter_doublet_gaps=gaps,
                 window=np.array([emin, emax]), labels=labels,
                 system_type=np.array([stype]),
                 note=np.array(["Kramers residual / inter-doublet gaps; "
                                "SOC spin splitting is not defined for a "
                                "PT-symmetric non-magnetic system"]))
        print(f"  -> {npz}")
        return 0
    split = adjacent_splitting(E)
    outdir = args.outdir or os.path.join(DEFAULT_OUT, cfg.get("name", "material"), "static")
    os.makedirs(outdir, exist_ok=True)
    ne = cfg.get("num_electrons")
    edge_pairs = gap_edge_pairs(ne, model.nw) if ne else []
    print(f"soc: {model.nw} bands, {len(kpts)} k-points, "
          f"window [{emin:.2f},{emax:.2f}] eV; "
          f"edge pairs (0-based) = {edge_pairs} "
          f"[VBM-1/VBM, CBM/CBM+1]")
    if not edge_pairs:
        # no num_electrons (e.g. analytic models): keep a general fallback over
        # the adjacent pairs inside the SOC window
        win = (E[:, :-1] >= emin) & (E[:, :-1] <= emax)
        edge_pairs = [(i, i + 1) for i in range(E.shape[1] - 1)
                      if win[:, i].any()][:6]
        print("  num_electrons missing -> fallback to adjacent pairs in window: "
              f"{edge_pairs}")
    edge_adj = edge_pair_splitting(E, edge_pairs) if edge_pairs else []
    for r in edge_adj:
        print(f"    bands {r['i']}-{r['j']}: adjacent splitting "
              f"max {r['max']*1000:.2f} meV, median {r['median']*1000:.2f} meV")
    # ---- spin splitting of the two gap-edge pairs only (VBM-1/VBM, CBM/CBM+1)
    edge_spin = []
    try:
        spin = load_spin_for_model(cfg, model)
        spin_at_k = spin.get("spin_at_k")
        ops = spin.get("ops")
        if (spin_at_k is not None or ops is not None) and edge_pairs:
            S = np.empty((len(kpts), model.nw, 3))
            for ik, k in enumerate(kpts):
                ops_k = spin_at_k(k) if spin_at_k is not None else ops
                S[ik] = spin_texture(model.Hk(k), ops_k)["spins"]
            edge_spin = edge_pair_splitting(E, edge_pairs, spins=S)
            print("  edge-pair spin splittings (VBM-1/VBM, CBM/CBM+1; "
                  "spectral splitting, not a calibrated SOC constant):")
            for r in edge_spin:
                print(f"    bands {r['i']}-{r['j']}: max {r['max']*1000:.2f} meV, "
                      f"median {r['median']*1000:.2f} meV, "
                      f"spin-antiparallel median "
                      f"{r['spin_antiparallel_median']:.3f}")
            _plot_spin_pairs(kpts, edge_spin, labels,
                             os.path.join(outdir, "soc_spin_splitting.png"),
                             lattice=_model_lattice(model))
    except (MissingTransformError, SpinInputError) as exc:
        print(f"  spin splitting skipped: {exc}")

    npz = os.path.join(outdir, "soc_splitting.npz")
    save = dict(kpts=kpts, energies=E, split=split,
                labels=np.array(labels if labels else []),
                window=np.array([emin, emax]))
    if edge_pairs:
        save["edge_pairs"] = np.array(edge_pairs)
        save["edge_pair_split"] = np.stack(
            [r["gap"] for r in edge_adj], axis=1)
        save["edge_pair_rows"] = np.array(
            [(r["i"], r["j"], r["max"], r["median"]) for r in edge_adj])
    if edge_spin:
        save["spin_pair_rows"] = np.array(
            [(r["i"], r["j"], r["max"], r["median"], len(kpts))
             for r in edge_spin])
    np.savez(npz, **save)
    print(f"  -> {npz}")
    if edge_adj:
        _plot_soc(kpts, E, labels, (emin, emax),
                  os.path.join(outdir, "soc_splitting.png"),
                  lattice=_model_lattice(model),
                  edge_pairs=edge_pairs, edge_rows=edge_adj)
    return 0


def _plot_spin_pairs(kpts, rows, labels, path, lattice=None):
    """Spin splitting of the two gap-edge pairs (VBM-1/VBM, CBM/CBM+1)."""
    import matplotlib.pyplot as plt
    kd = ps.kpath_distance(kpts, lattice)
    F = ps.FONTS
    with ps.apply_style():
        fig, ax = plt.subplots(figsize=(10.0 / 2.54, 5.8 / 2.54),
                               constrained_layout=True)
        for n, r in enumerate(rows):
            lab = (f"{r['i']}-{r['j']}  (max {r['max'] * 1000:.1f}, "
                   f"med {r['median'] * 1000:.2f})")
            ax.plot(kd, r["gap"] * 1000.0, "-o", ms=2.5, lw=0.7, alpha=0.9,
                    color=ps.PALETTE[n % len(ps.PALETTE)], label=lab)
        ps.add_reference_line(ax, 0.0)
        ax.set_ylabel("edge-pair splitting (meV)", fontsize=F["label"])
        if labels:
            ps.add_sym_labels(ax, kd, labels)
        else:
            ax.set_xlim(kd[0], kd[-1])
            ax.set_xlabel("k distance (1/$\\AA$)", fontsize=F["label"])
        ax.legend(frameon=True, edgecolor="black", loc="upper left",
                  fontsize=F["legend"])
        ps.save_figure(fig, path)
    print(f"  figure: {path}")


def _plot_soc(kpts, E, labels, window, path, lattice=None,
              edge_pairs=(), edge_rows=()):
    """Gap-edge SOC figures in the phonkit PRL style.

    Only the two gap-edge pairs are shown.  Bands are drawn as CONTINUOUS
    lines over the whole k path: the visible energy window only sets the y
    limits (no NaN masking, which used to make the curves look broken).
    """
    import matplotlib.pyplot as plt
    kd = ps.kpath_distance(kpts, lattice)
    colors = [ps.COLOR_REF, ps.COLOR_CMP]
    pair_names = ["valence pair", "conduction pair"]
    F = ps.FONTS
    vbm = (float(E[:, edge_pairs[0][1]].max()) if edge_pairs else 0.0)
    # y limits: cover the SOC window AND the full gap-edge bands
    ylo = window[0] - vbm
    yhi = window[1] - vbm
    if edge_pairs:
        idx = [ib for p in edge_pairs for ib in p]
        vals = E[:, idx] - vbm
        ylo = min(ylo, float(vals.min()))
        yhi = max(yhi, float(vals.max()))
    pad = 0.05 * max(yhi - ylo, 1e-6)
    ylo, yhi = ylo - pad, yhi + pad
    with ps.apply_style():
        fig, axes = plt.subplots(1, 2, figsize=ps.FIG_2PANEL,
                                 constrained_layout=True)
        ax = axes[0]
        for ib in range(E.shape[1]):
            band = E[:, ib] - vbm            # full continuous curve
            if not ((band >= ylo) & (band <= yhi)).any():
                continue
            ax.plot(kd, band, color=ps.COLOR_SINGLE, lw=0.5, alpha=0.45)
        for n, (i, j) in enumerate(edge_pairs):
            color = colors[n % len(colors)]
            for ib in (i, j):
                ax.plot(kd, E[:, ib] - vbm, color=color, lw=0.9,
                        label=(pair_names[n] if ib == i and n < len(pair_names)
                               else None))
        if len(edge_pairs) >= 2:
            gap_curve = (E[:, edge_pairs[1][0]] - E[:, edge_pairs[0][1]])
            ax.annotate(f"min direct gap {gap_curve.min():.3f} eV",
                        xy=(0.03, 0.52), xycoords="axes fraction",
                        fontsize=F["legend"])
        ax.set_ylim(ylo, yhi)
        ax.set_ylabel("Energy $-$ $E_{\\rm VBM}$ (eV)", fontsize=F["label"])
        if edge_pairs:
            ax.legend(fontsize=F["legend"], frameon=True, edgecolor="black",
                      loc="lower right")
        ax = axes[1]
        for n, r in enumerate(edge_rows):
            lab = (f"{r['i']}-{r['j']}  (max {r['max'] * 1000:.1f}, "
                   f"med {r['median'] * 1000:.2f})")
            ax.plot(kd, r["gap"] * 1000.0, "-", lw=0.8,
                    color=colors[n % len(colors)], label=lab)
        ps.add_reference_line(ax, 0.0)
        ax.set_ylabel("$\\Delta E$ (meV)", fontsize=F["label"])
        if edge_rows:
            ax.legend(fontsize=F["legend"], frameon=True, edgecolor="black",
                      loc="upper left")
        for ax in axes:
            if labels:
                ps.add_sym_labels(ax, kd, labels)
            else:
                ax.set_xlim(kd[0], kd[-1])
        ps.save_figure(fig, path)
    print(f"  figure: {path}")


# ---------------------------------------------------------------- rashba fit
def _reciprocal_basis(lattice):
    """Reciprocal basis rows b1,b2,b3 (1/Angstrom) for a row-vector lattice."""
    A = np.asarray(lattice, float)
    return 2.0 * np.pi * np.linalg.inv(A).T


def _resolve_rashba_anchor(cfg):
    """Resolve (anchor_frac, direction_frac, anchor_label, direction_label)."""
    rf = cfg.get("rashba_fit") or {}
    kpath = cfg.get("kpath") or []
    labels = [str(p[0]).replace("GAMMA", "\u0393") for p in kpath]
    nodes = [np.asarray(p[1], float) for p in kpath]
    anchor = str(rf.get("anchor", labels[0] if labels else "\u0393"))
    anchor = anchor.replace("GAMMA", "\u0393")
    ia = labels.index(anchor) if anchor in labels else 0
    direction = rf.get("direction")
    if isinstance(direction, str):
        dl = direction.replace("GAMMA", "\u0393")
        if dl not in labels:
            raise SystemExit(
                f"rashba_fit: direction label {direction!r} is not on kpath")
        d_frac = nodes[labels.index(dl)] - nodes[ia]
        dir_label = dl
    elif direction is not None:
        d_frac = np.asarray(direction, float)
        dir_label = str(rf.get("direction_label", ""))
    else:
        ia2 = (ia + 1) % len(nodes)
        d_frac = nodes[ia2] - nodes[ia]
        dir_label = labels[ia2]
    return nodes[ia], d_frac, labels[ia], dir_label


def cmd_rashba(args):
    """Linear Rashba fit of the gap-edge pairs along a line through an anchor.

    Model (same as the user's ``rashba_refined_fitting.py``):
        E(k) = a路k虏 + E0 + 伪路|k|,   k measured from the high-symmetry anchor.
    The fit window grows while both branches keep R虏 >= r2_thr and the two
    fitted |伪| stay consistent; 伪_R = (|伪_lower|+|伪_upper|)/2 in eV路脜.
    Config section ``rashba_fit``: anchor, direction (label or vector),
    nk, kmax, kmin, step, r2_thr, alpha_ratio_thr.
    """
    cfg = load_config(args.config)
    model = build_model(cfg)
    name = cfg.get("name", "material")
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "static")
    os.makedirs(outdir, exist_ok=True)
    rf = dict(cfg.get("rashba_fit") or {})
    ne = cfg.get("num_electrons")
    pairs = gap_edge_pairs(ne, model.nw) if ne else []
    if not pairs:
        pairs = [(0, 1)]
        print("rashba: num_electrons missing -> fitting the lowest adjacent pair")
    nk = int(rf.get("nk", 161))
    kmax = float(rf.get("kmax", 0.15))
    kmin = float(rf.get("kmin", 0.03))
    step = float(rf.get("step", 0.005))
    r2_thr = float(rf.get("r2_thr", 0.999))
    ratio_thr = float(rf.get("alpha_ratio_thr", 0.15))
    split_thr = float(rf.get("split_r2_thr", 0.99))
    fallback_kmax = float(rf.get("fallback_kmax", 0.08))
    lattice = _model_lattice(model)
    t = np.linspace(-kmax, kmax, nk)
    if lattice is None:
        anchor_label, dir_label = "\u0393", "k"
        kpts = np.column_stack([t, np.zeros_like(t), np.zeros_like(t)])
    else:
        a_frac, d_frac, anchor_label, dir_label = _resolve_rashba_anchor(cfg)
        B = _reciprocal_basis(lattice)
        d_cart = np.asarray(d_frac, float) @ B
        nd = float(np.linalg.norm(d_cart))
        if nd < 1e-12:
            raise SystemExit("rashba_fit: direction vector is zero")
        d_hat = d_cart / nd
        kcart = np.asarray(a_frac, float) @ B + t[:, None] * d_hat[None, :]
        kpts = kcart @ np.linalg.inv(B)
    E = model.static_bands(kpts)
    print(f"rashba: anchor {anchor_label}, direction {dir_label}, "
          f"k = [-{kmax:.3f}, +{kmax:.3f}] 1/A ({nk} points), pairs = {pairs}")
    results = []
    for (i, j) in pairs:
        r = fit_gap_edge_pair(t, np.stack([E[:, i], E[:, j]], axis=0),
                              kmin=kmin, kmax_max=kmax, step=step,
                              r2_thr=r2_thr, alpha_ratio_thr=ratio_thr,
                              split_r2_thr=split_thr,
                              fallback_kmax=fallback_kmax)
        r.update(dict(i=int(i), j=int(j)))
        results.append(r)
        if r["status"] == "LINEAR":
            print(f"  pair {i}-{j}: LINEAR  alpha_R = "
                  f"{r['alpha_r'] * 1000:.1f} meV*A ({r['alpha_r']:.4f} eV*A), "
                  f"k_max = {r['kmax']:.3f} 1/A, "
                  f"R2(branch) = {r['r2_lower']:.5f}/{r['r2_upper']:.5f}, "
                  f"R2(splitting) = {r['r2_split']:.5f}")
        elif r["status"] == "FALLBACK":
            print(f"  pair {i}-{j}: FALLBACK (auto window failed; fitting at "
                  f"k_max = {r['kmax']:.3f} like the original script; "
                  f"branch R2 = {r['r2_lower']:.4f}/{r['r2_upper']:.4f}; "
                  f"splitting R2 = {r['r2_split']:.4f} -> not linear) "
                  f"alpha_R = {r['alpha_r'] * 1000:.1f} meV*A")
        else:
            print(f"  pair {i}-{j}: NOT_LINEAR (no window satisfies "
                  f"branch R2 >= {r2_thr}, opposite-sign alpha, "
                  f"|alpha| ratio <= {ratio_thr}, splitting R2 >= {split_thr}; "
                  f"fallback fit not possible)")
    npz = os.path.join(outdir, "rashba_fit.npz")
    save = dict(k=t, anchor=np.array([anchor_label]),
                direction=np.array([dir_label]),
                pairs=np.array(pairs, dtype=int),
                alpha_r=np.array([r["alpha_r"] if r["alpha_r"] is not None
                                  else np.nan for r in results]),
                kmax=np.array([r["kmax"] if r["kmax"] is not None
                               else np.nan for r in results]),
                alpha_from_split=np.array(
                    [r.get("alpha_from_split", np.nan)
                     if r["status"] in ("LINEAR", "FALLBACK") else np.nan
                     for r in results]),
                r2_split=np.array(
                    [r.get("r2_split", np.nan)
                     if r["status"] in ("LINEAR", "FALLBACK") else np.nan
                     for r in results]),
                status=np.array([r["status"] for r in results]))
    for n, r in enumerate(results):
        save[f"E_lower_{n}"] = r["lower"]
        save[f"E_upper_{n}"] = r["upper"]
        if r["status"] in ("LINEAR", "FALLBACK"):
            save[f"popt_lower_{n}"] = r["popt_lower"]
            save[f"popt_upper_{n}"] = r["popt_upper"]
    np.savez(npz, **save)
    print(f"  -> {npz}")
    _plot_rashba_fit(results, anchor_label, dir_label,
                     os.path.join(outdir, "rashba_fit.png"))
    return 0


def _plot_rashba_fit(results, anchor_label, dir_label, path):
    """Rashba fit figure following the user's rashba_refined_fitting.py style.

    * two stacked panels: CBM pair on top, VBM pair at the bottom (when both);
    * data: hollow red circles, N_BRANCH points per branch, mirrored to -k;
    * k.p model: pure blue lines; alpha_R annotation inside the parabola opening;
    * symmetric -k..+k x axis with anchor/direction labels, y in eV (%.3f);
    * direction arrows below the axis; saved at 300 dpi like the original.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FormatStrFormatter

    # global framework typography (same relative sizes as every other figure)
    COLOR_DATA, COLOR_FIT, EDGE = "none", "#0000FF", "#CC0000"
    N_BRANCH, S_POINT = 25, 28.0

    order = list(range(len(results)))
    openings = ["up"] * len(results)
    if len(results) == 2:            # results = [valence, conduction] -> CBM top
        order = [1, 0]
        openings = ["up", "down"]
    with plt.rc_context(ps.STYLE):
        fig, axes = plt.subplots(len(order), 1,
                                 figsize=(10.5 / 2.54, 10.0 / 2.54),
                                 constrained_layout=True, squeeze=False)
        axes = axes[:, 0]
        for ax, idx, opening in zip(axes, order, openings):
            r = results[idx]
            k = np.asarray(r["k"], float)
            kmax = (float(r["kmax"]) if r.get("kmax")
                    else float(np.max(np.abs(k))))
            pos = k >= 0
            kpos = k[pos]
            k_new = np.linspace(0.0, kmax, N_BRANCH)
            for branch in ("lower", "upper"):
                e = np.asarray(r[branch], float)[pos]
                e_new = np.interp(k_new, kpos, e)
                ax.scatter(k_new, e_new, s=S_POINT, facecolors=COLOR_DATA,
                           edgecolors=EDGE, linewidth=1.2, zorder=3)
                m = k_new > 1e-9
                ax.scatter(-k_new[m], e_new[m], s=S_POINT,
                           facecolors=COLOR_DATA, edgecolors=EDGE,
                           linewidth=1.2, zorder=3)
            if r["status"] in ("LINEAR", "FALLBACK"):
                kd = np.linspace(-kmax, kmax, 400)
                for key in ("popt_lower", "popt_upper"):
                    ax.plot(kd, ps_model(r[key], kd), color=COLOR_FIT,
                            lw=1.2, zorder=2)
            ax.axvline(0.0, color="black", ls="--", lw=0.4, alpha=0.4)
            xlim_half = kmax / 0.75
            ax.set_xlim(-xlim_half, xlim_half)
            ax.set_xticks([-xlim_half, 0.0, xlim_half])
            ax.set_xticklabels(["-" + dir_label, anchor_label, dir_label])
            ax.set_ylabel("Energy (eV)", fontsize=ps.FONTS["label"])
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
            # data occupies ~60% of the y range, margin biased to the opening
            e_win = np.concatenate(
                [np.asarray(r[b], float)[np.abs(k) <= kmax]
                 for b in ("lower", "upper")])
            lo, hi = float(e_win.min()), float(e_win.max())
            span = hi - lo if hi > lo else 1.0
            margin = span / 0.60 - span
            if opening == "up":
                m_top, m_bot = 0.7 * margin, 0.3 * margin
                leg_anchor = (0.5, 0.78)
            else:
                m_top, m_bot = 0.3 * margin, 0.7 * margin
                leg_anchor = (0.5, 0.22)
            ax.set_ylim(lo - m_bot, hi + m_top)
            # exactly the original script's legend labels/handles
            h_data = Line2D([], [], marker="o", markerfacecolor="none",
                            markeredgecolor=EDGE, markeredgewidth=1.2,
                            color="none", label="DFT")
            h_fit = Line2D([], [], color=COLOR_FIT, lw=1.2,
                           label="k\u00b7p model")
            handles = [h_data, h_fit] if r["status"] == "LINEAR" else [h_data]
            leg = ax.legend(handles=handles, loc="center",
                            bbox_to_anchor=leg_anchor, frameon=False,
                            fontsize=11)
            if r["status"] in ("LINEAR", "FALLBACK"):
                # alpha_R under the legend, left-aligned with its symbol column
                # (same placement as the original script's place_alpha)
                fig.canvas.draw()
                renderer = fig.canvas.get_renderer()
                bb = leg.get_window_extent(renderer)
                lines = leg.get_lines()
                mx = (min(l.get_window_extent(renderer).x0 for l in lines)
                      if lines else bb.x0)
                inv = ax.transAxes.inverted()
                tx, ty = inv.transform((mx, bb.y0))
                ax.text(tx, ty - 0.018,
                        f"$\\alpha_R$ = {r['alpha_r'] * 1000:.1f} meV\u00b7\u00c5",
                        ha="left", va="top", transform=ax.transAxes,
                        fontsize=ps.FONTS["legend"], clip_on=False)
            # direction arrows below the axis, exactly as the original script:
            # x in data units, y in axes units
            from matplotlib.transforms import blended_transform_factory
            trans = blended_transform_factory(ax.transData, ax.transAxes)
            xh = xlim_half
            for sgn in (-1.0, 1.0):
                ax.annotate("", xy=(sgn * 0.70 * xh, -0.08),
                            xytext=(sgn * 0.30 * xh, -0.08),
                            xycoords=trans, textcoords=trans,
                            annotation_clip=False,
                            arrowprops=dict(arrowstyle="->", color="black",
                                            lw=1.2, mutation_scale=14))
        fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
    print(f"  figure: {path}")


def ps_model(popt, k):
    """Local helper: E(k)=a k虏 + E0 + 伪|k| for plotting (kept in cli)."""
    a, e0, alpha = popt
    k = np.asarray(k, float)
    return a * k * k + e0 + alpha * np.abs(k)


# ------------------------------------------------------------------ pst
def _pst_not_applicable(args, cfg):
    """Pt (PT-symmetric): single-band PST is gauge dependent -> NOT_APPLICABLE."""
    from .observables.kramers import kramers_summary
    model, spin = load_model_spin_cached(cfg)
    spin_at_k = spin.get("spin_at_k") or spin.get("ops")
    kpts = np.asarray(model.win.get("kpoints", []), float)
    if len(kpts) == 0:
        print("BLOCKED: no DFT k-mesh available for the gauge test.")
        return 2
    E = model.static_bands(kpts)
    ks = kramers_summary(E, kpts)
    print("pst: PT-symmetric non-magnetic system -> PST is NOT APPLICABLE.")
    print(f"  DFT-mesh delta_PT: max {ks['delta_max']*1000:.4f} meV, median "
          f"{ks['delta_median']*1000:.4f} meV")
    bands = _pst_bands(cfg, model)
    g = gauge_invariance_test(model, spin_at_k, kpts, bands, deg_tol=1e-6)
    print(f"  gauge test on the true DFT Kramers pairs ({g['n_pairs']} pair "
          f"rotations): collinearity {g['coll_range'][0]:.4f}-"
          f"{g['coll_range'][1]:.4f} (base {g['base']['coll']:.4f}), "
          f"axis dot min {g['axis_dot_min']:.3f}")
    print(f"  gauge-invariant subspace spin: mean |0.5 Tr(P sigma)| = "
          f"{g['invariant_mean_norm']:.4f}, max {g['invariant_max_norm']:.4f}")
    outdir = args.outdir or os.path.join(DEFAULT_OUT, cfg.get("name", "material"), "static")
    os.makedirs(outdir, exist_ok=True)
    npz = os.path.join(outdir, "pst.npz")
    np.savez(npz, kpts=kpts, energies=E,
             eigenvalues=np.array(g["base"]["M"]),
             axes=np.array(g["base"]["axis"]),
             pst_axis=np.array(g["base"]["axis"]),
             collinearity=np.array(g["base"]["coll"]),
             gauge_collinearity=np.array([m["coll"] for m in g["seeds"]]),
             gauge_axis=np.array([m["axis"] for m in g["seeds"]]),
             gauge_invariant_spin=g["invariant_spin"],
             bands=np.array(bands),
             system_type=np.array(["pt_nonmagnetic"]),
             note=np.array(["PST NOT_APPLICABLE: gauge-dependent single-band "
                            "spins in Kramers-degenerate PT system"]))
    print(f"  -> {npz}")
    return 0


def cmd_pst(args):
    cfg = load_config(args.config)
    name = cfg.get("name", "material")
    if (cfg.get("system") or {}).get("type") == "pt_nonmagnetic":
        return _pst_not_applicable(args, cfg)
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "static")
    os.makedirs(outdir, exist_ok=True)
    n_km, kmax_km = resolve_kmesh(cfg, "kmesh_pst", default_n=21,
                                  default_kmax=0.55)
    if getattr(args, "kmesh_n", None):
        n_km = int(args.kmesh_n)
    if getattr(args, "kmax", None) is not None:
        kmax_km = float(args.kmax)
    kpts = plane_kpts(cfg.get("kplane"), n_km, kmax_km)
    jobs = _jobs_for(args, cfg, len(kpts))
    if jobs > 1 and _parallel_ok(cfg):
        lat = _lattice_from_cfg(cfg)
        nw = _nw_from_cfg(cfg)
        print(f"pst mesh: {len(kpts)} k-points (n={n_km}, "
              f"|k_i| <= {kmax_km}), jobs={jobs}")
        load_model_spin_cached(cfg)   # warm cache for workers
        tex = run_chunks("pst", os.path.abspath(args.config), kpts, jobs)
        E, S = tex["energies"], tex["spins"]
        pst_cfg = cfg.get("pst", {})
        bands = ([int(b) for b in pst_cfg["bands"]]
                 if pst_cfg.get("bands") else
                 ([int(x) for x in range(max(0, (cfg.get("num_electrons") or 1) - 2),
                                         min(nw, (cfg.get("num_electrons") or nw) + 2))]
                  if cfg.get("num_electrons") else list(range(min(4, nw)))))
        res = analyse_texture(
            E, S, bands=bands,
            min_polarization=pst_cfg.get("min_polarization", 0.1),
            min_splitting=pst_cfg.get("min_splitting", 1e-6))
        print(f"pst: {len(kpts)} k-points, bands {bands}")
        print(f"  M eigenvalues = {np.round(res['eigenvalues'], 4)}")
        print(f"  collinearity (largest eigenvalue) = {res['collinearity']:.4f}")
        print(f"  PST axis = {np.round(res['pst_axis'], 4)}")
        print(f"  n_states = {res['n_states']}, mean polarization = "
              f"{res['mean_polarization']:.3f}")
        npz = os.path.join(outdir, "pst.npz")
        np.savez(npz, kpts=tex["kpts"], energies=E, spins=S,
                 eigenvalues=res["eigenvalues"], axes=res["axes"],
                 pst_axis=res["pst_axis"], collinearity=res["collinearity"],
                 bands=np.array(bands))
        gauge = None
        try:
            model, spin = load_model_spin_cached(cfg)
            satk = spin.get("spin_at_k")
            op_obj = spin.get("spin_operator")
            if satk is not None and op_obj is not None:
                nsub = min(40, len(kpts))
                idx = np.linspace(0, len(kpts) - 1, nsub).astype(int)
                gauge = gauge_invariance_test(
                    model, satk, kpts[idx], bands,
                    seeds=(1, 2, 3, 4, 5), deg_tol=5e-3)
                print(f"  gauge test (near-degenerate subspaces): "
                      f"collinearity spread {gauge['coll_spread']:.4f}, "
                      f"axis dot min {gauge['axis_dot_min']:.3f}, "
                      f"invariant |s| mean {gauge['invariant_mean_norm']:.3f}")
                if gauge["coll_spread"] > 1e-2 or gauge["axis_dot_min"] < 0.99:
                    print("  WARN: PST metric is gauge sensitive; the axis "
                          "should be quoted only with the subspace caveat.")
        except Exception as exc:                        # noqa: BLE001
            print("  gauge test unavailable:", exc)
        # re-save with gauge fields
        extra = {}
        if gauge is not None:
            extra = dict(
                gauge_coll_spread=np.array(gauge["coll_spread"]),
                gauge_axis_dot_min=np.array(gauge["axis_dot_min"]),
                gauge_invariant_mean_norm=np.array(
                    gauge["invariant_mean_norm"]),
                gauge_invariant_spin=gauge["invariant_spin"])
            np.savez(npz, kpts=tex["kpts"], energies=E, spins=S,
                     eigenvalues=res["eigenvalues"], axes=res["axes"],
                     pst_axis=res["pst_axis"],
                     collinearity=res["collinearity"],
                     bands=np.array(bands), **extra)
        print(f"  -> {npz}")
        ib = bands[0] if bands else 0
        _plot_spin_texture(tex["kpts"], S, os.path.join(outdir, "pst.png"),
                           title=f"{name}: spin texture for PST (band {ib})",
                           iband=ib, lattice=lat, kplane=cfg.get("kplane"))
        return 0
    model, spin = load_model_spin_cached(cfg)
    spin_at_k = spin.get("spin_at_k")
    ops = spin.get("ops")
    if spin_at_k is None and ops is None:
        print("BLOCKED: no gauge-consistent spin operator available.")
        return 2
    lat = _model_lattice(model)
    kplane = cfg.get("kplane")
    op_obj = spin.get("spin_operator")
    if spin_at_k is not None and op_obj is not None and hasattr(op_obj, "S_R"):
        tex = spin_texture_wannier_batch(model, op_obj, kpts, chunk=128)
        E, S = tex["energies"], tex["spins"]
    else:
        E = np.empty((len(kpts), model.nw))
        S = np.empty((len(kpts), model.nw, 3))
        for ik, k in enumerate(kpts):
            Hk = model.Hk(k)
            ops_k = spin_at_k(k) if spin_at_k is not None else ops
            r = spin_texture(Hk, ops_k)
            E[ik] = r["energies"]
            S[ik] = r["spins"]
    bands = _pst_bands(cfg, model)
    res = analyse_texture(
        E, S, bands=bands,
        min_polarization=cfg.get("pst", {}).get("min_polarization", 0.1),
        min_splitting=cfg.get("pst", {}).get("min_splitting", 1e-6))
    print(f"pst: {len(kpts)} k-points, bands {bands}")
    print(f"  M eigenvalues = {np.round(res['eigenvalues'], 4)}")
    print(f"  collinearity (largest eigenvalue) = {res['collinearity']:.4f}")
    print(f"  PST axis = {np.round(res['pst_axis'], 4)}")
    print(f"  n_states = {res['n_states']}, "
          f"mean polarization = {res['mean_polarization']:.3f}")
    gauge = None
    try:
        op_obj = spin.get("spin_operator")
        satk = spin_at_k if spin_at_k is not None else None
        if satk is not None and op_obj is not None:
            nsub = min(40, len(kpts))
            idx = np.linspace(0, len(kpts) - 1, nsub).astype(int)
            gauge = gauge_invariance_test(
                model, satk, kpts[idx], bands,
                seeds=(1, 2, 3, 4, 5), deg_tol=5e-3)
            print(f"  gauge test (near-degenerate subspaces, {nsub} k): "
                  f"collinearity spread {gauge['coll_spread']:.4f}, "
                  f"axis dot min {gauge['axis_dot_min']:.3f}, "
                  f"invariant |s| mean {gauge['invariant_mean_norm']:.3f}")
    except Exception as exc:                            # noqa: BLE001
        print("  gauge test unavailable:", exc)
    npz = os.path.join(outdir, "pst.npz")
    extra = {}
    if gauge is not None:
        extra = dict(gauge_coll_spread=np.array(gauge["coll_spread"]),
                     gauge_axis_dot_min=np.array(gauge["axis_dot_min"]),
                     gauge_invariant_mean_norm=np.array(
                         gauge["invariant_mean_norm"]),
                     gauge_invariant_spin=gauge["invariant_spin"])
    np.savez(npz, kpts=kpts, energies=E, spins=S,
             eigenvalues=res["eigenvalues"], axes=res["axes"],
             pst_axis=res["pst_axis"], collinearity=res["collinearity"],
             bands=np.array(bands), **extra)
    print(f"  -> {npz}")
    ib = bands[0] if bands else 0
    _plot_spin_texture(kpts, S, os.path.join(outdir, "pst.png"),
                       title=f"{name}: spin texture for PST (band {ib})",
                       iband=ib, lattice=lat, kplane=kplane)
    return 0


# ------------------------------------------------------------------ stage D: scan
def cmd_scan(args):
    cfg = load_config(args.config)
    model = build_model(cfg)
    name = cfg.get("name", "material")
    spin = _load_spin_or_exit(cfg, model)
    if spin is None:
        return 2
    spin_arg = spin.get("spin_at_k") if spin.get("spin_at_k") is not None \
        else spin.get("ops")
    if spin_arg is None:
        print("BLOCKED: no spin operator available for the scan observable.")
        return 2
    drive_cfg = dict(cfg.get("drive") or dict(polarization="circular",
                                              sign=-1))
    drive_cfg.pop("hw", None)
    drive_cfg.pop("gamma", None)
    sc = dict(cfg.get("scan") or {})
    hw = args.hw if args.hw else sc.get("hw", [0.2, 0.4, 0.6])
    gamma = args.gamma if args.gamma else sc.get("gamma", [0.0, 0.05, 0.10])
    kpt = np.atleast_2d(np.asarray(sc.get("kpt", [[0.05, 0.03, 0.0]]), float))
    n_step = args.n_step or int(sc.get("n_step", 150))
    m_sp = args.m_sp or int(sc.get("m_sp", 16))
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "floquet")
    os.makedirs(outdir, exist_ok=True)
    print(f"scan: {model.nw} bands, hw={hw} eV, gamma={gamma}, "
          f"{len(kpt)} k-point(s), n_step={n_step}, m_sp={m_sp}")
    window = tuple(sc.get("window", (-2.0, 2.0)))
    rows = run_scan(model, drive_cfg, hw, gamma, kpt, spin_arg,
                    n_step=n_step, m_sp=m_sp, window=window)
    arr = rows_to_arrays(rows)
    npz = os.path.join(outdir, "scan.npz")
    np.savez(npz, kpts=kpt, **arr)
    print("    hw    gamma   eps_spread  mean|S|  max|S|  wsum_err")
    for r in rows:
        if r["ik"] == 0:
            print(f"  {r['hw']:.2f}  {r['gamma']:.3f}   {r['eps_spread']:.4f}"
                  f"     {r['mean_pol']:.3f}   {r['max_pol']:.3f}"
                  f"   {r['wsum_err']:.1e}")
    print(f"  -> {npz}")
    _plot_scan(rows, os.path.join(outdir, "scan.png"))
    return 0


def _plot_scan(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ps.enforce_style()
    hws = np.array(sorted({r["hw"] for r in rows}), float)
    gammas = np.array(sorted({r["gamma"] for r in rows}), float)
    is_grid = len(hws) >= 2 and len(gammas) >= 2
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2),
                             constrained_layout=True)
    for ax, q, lab in ((axes[0], "mean_pol", r"mean $|\bar S|$"),
                       (axes[1], "wsum_err", r"max $|\sum_m w - 1|$")):
        if is_grid:
            _, _, Z = observable_map(rows, q, ik=0)
            pcm = ax.pcolormesh(hws, gammas, Z, cmap="viridis",
                                shading="auto")
            ps.add_colorbar(fig, pcm, ax, label=lab)
            ax.set_ylabel(r"$\gamma = A_0 a_{ref}$",
                          fontsize=ps.FONTS["label"])
            ax.set_xlabel(r"$\hbar\omega$ (eV)", fontsize=ps.FONTS["label"])
        else:
            # 1-D or single-point scan: pcolormesh has no cells to draw and
            # would silently produce an empty figure, so plot markers instead.
            _plot_scan_1d(ax, rows, q, lab, hws, gammas)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  figure: {path}")


def _plot_scan_1d(ax, rows, quantity, label, hws, gammas):
    """Marker plot for scans that vary one parameter (or none)."""
    rows = [r for r in rows if r.get("ik", 0) == 0]
    if len(hws) >= 2:
        x = np.array([r["hw"] for r in rows], float)
        xlabel = r"$\hbar\omega$ (eV)"
    elif len(gammas) >= 2:
        x = np.array([r["gamma"] for r in rows], float)
        xlabel = r"$\gamma = A_0 a_{ref}$"
    else:
        x = np.zeros(len(rows))
        xlabel = "scan point"
    y = np.array([r[quantity] for r in rows], float)
    ok = np.isfinite(y)
    ax.plot(x[ok], y[ok], "o-", color="#0072B2", lw=1.2, ms=5)
    if ok.sum() == 1:
        ax.annotate(f"{y[ok][0]:.4g}", (x[ok][0], y[ok][0]),
                    textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax.set_xlabel(xlabel, fontsize=ps.FONTS["label"])
    ax.set_ylabel(label, fontsize=ps.FONTS["label"])


# ------------------------------------------------------------------ stage D: converge
def cmd_converge(args):
    cfg = load_config(args.config)
    model = build_model(cfg)
    name = cfg.get("name", "material")
    spin = _load_spin_or_exit(cfg, model)
    if spin is None:
        return 2
    spin_arg = spin.get("spin_at_k") if spin.get("spin_at_k") is not None \
        else spin.get("ops")
    if spin_arg is None:
        print("BLOCKED: no spin operator available.")
        return 2
    drive = drive_from_cfg(cfg, model)
    cv = dict(cfg.get("converge") or {})
    kpt = np.atleast_2d(np.asarray(cv.get("kpt", [[0.05, 0.03, 0.0]]), float))
    k = kpt[0]
    n_steps = args.n_step or cv.get("n_step", [60, 100, 200, 400])
    m_sps = args.m_sp or cv.get("m_sp", [4, 8, 16, 32])
    outdir = args.outdir or os.path.join(DEFAULT_OUT, name, "floquet")
    os.makedirs(outdir, exist_ok=True)
    ops_k = spin_arg(k) if callable(spin_arg) else spin_arg
    print(f"converge: {model.nw} bands at k={np.round(k, 3)}, "
          f"hw={drive.hw:g}, gamma={drive.gamma:g}")
    window = tuple(cv.get("window", (-2.0, 2.0)))
    res = run_convergence(model, drive, k, ops_k, n_steps, m_sps,
                          window=window,
                          tol_ebar=float(cv.get("tol_ebar", 1e-5)),
                          tol_spin=float(cv.get("tol_spin", 1e-4)),
                          tol_wsum=float(cv.get("tol_wsum", 1e-6)))
    ref = res["reference"]
    print(f"  reference: n_step={ref['n_step']}, m_sp={ref['m_sp']}")
    print("  n_step convergence:")
    for r in res["n_step_rows"]:
        print(f"    n_step={r['n_step']:4d}  dEbar={r['dev_ebar']:.2e} eV  "
              f"dS={r['dev_spin']:.2e}  wsum_err={r['wsum_err']:.1e}  "
              f"matched={r['n_matched']}")
    print("  m_sp convergence:")
    for r in res["m_sp_rows"]:
        print(f"    m_sp={r['m_sp']:3d}    dEbar={r['dev_ebar']:.2e} eV  "
              f"dS={r['dev_spin']:.2e}  dWsum={r['dev_wsum']:.2e}  "
              f"w_def(max/med/n>1e-3)={r.get('w_deficit_max', float('nan')):.2e}/"
              f"{r.get('w_deficit_median', float('nan')):.2e}/"
              f"{r.get('n_deficit_gt_1e3', -1)}  matched={r['n_matched']}")
    npz = os.path.join(outdir, "converge.npz")
    np.savez(npz,
             n_step=np.array([r["n_step"] for r in res["n_step_rows"]]),
             dev_ebar=[r["dev_ebar"] for r in res["n_step_rows"]],
             dev_spin=[r["dev_spin"] for r in res["n_step_rows"]],
             m_sp=np.array([r["m_sp"] for r in res["m_sp_rows"]]),
             dev_ebar_m=[r["dev_ebar"] for r in res["m_sp_rows"]],
             dev_spin_m=[r["dev_spin"] for r in res["m_sp_rows"]],
             dev_wsum_n=[r["dev_wsum"] for r in res["n_step_rows"]],
             dev_wsum_m=[r["dev_wsum"] for r in res["m_sp_rows"]],
             max_dev_wsum=np.array([res.get("max_dev_wsum", np.nan)]),
             w_deficit_max_m=[r.get("w_deficit_max", np.nan)
                              for r in res["m_sp_rows"]],
             w_deficit_median_m=[r.get("w_deficit_median", np.nan)
                                 for r in res["m_sp_rows"]],
             n_deficit_gt_1e3_m=[r.get("n_deficit_gt_1e3", -1)
                                  for r in res["m_sp_rows"]],
             accepted_m_sp=np.array([res.get("accepted_m_sp") or -1]),
             status=np.array(["PASS" if res["pass_"] else "FAIL"]))
    print(f"  -> {npz}")
    _plot_converge(res, os.path.join(outdir, "converge.png"))
    print(f"  invariants: sumEbar={res['invariants']['sum_ebar']:.6f} eV, "
          f"sumS={np.round(res['invariants']['sum_sbar'], 6)}, "
          f"sumW={res['invariants']['sum_w']:.6f}")
    print(f"  convergence: {'PASS' if res['pass_'] else 'FAIL'} "
          f"(production={res.get('production')}, "
          f"max dEbar={res['max_dev_ebar']:.2e} eV, "
          f"max dS={res['max_dev_spin']:.2e}, "
          f"max dWsum={res.get('max_dev_wsum', float('nan')):.2e}, "
          f"accepted m_sp>={res.get('accepted_m_sp')}, "
          f"max invariant dev={res['max_dev_invariant']:.2e})")
    return 0 if res["pass_"] else 1


def _plot_converge(res, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ps.enforce_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0),
                             constrained_layout=True)
    ax = axes[0]
    ax.loglog([r["n_step"] for r in res["n_step_rows"]],
              [max(r["dev_ebar"], 1e-18) for r in res["n_step_rows"]],
              "o-", label=r"$\Delta\bar E$ (eV)")
    ax.loglog([r["n_step"] for r in res["n_step_rows"]],
              [max(r["dev_spin"], 1e-18) for r in res["n_step_rows"]],
              "s-", label=r"$\Delta\bar S$")
    ax.set_xlabel("n_step", fontsize=ps.FONTS["label"])
    ax.legend(fontsize=ps.FONTS["legend"], frameon=False)
    ax = axes[1]
    ax.loglog([r["m_sp"] for r in res["m_sp_rows"]],
              [max(r["dev_ebar"], 1e-18) for r in res["m_sp_rows"]],
              "o-", label=r"$\Delta\bar E$ (eV)")
    ax.loglog([r["m_sp"] for r in res["m_sp_rows"]],
              [max(r["dev_spin"], 1e-18) for r in res["m_sp_rows"]],
              "s-", label=r"$\Delta\bar S$")
    ax.set_xlabel("m_sp", fontsize=ps.FONTS["label"])
    ax.legend(fontsize=ps.FONTS["legend"], frameon=False)
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  figure: {path}")


# ------------------------------------------------------------------ benchmark
def cmd_benchmark(args):
    from .benchmarks import BENCHMARKS, run_all
    outdir = args.outdir or os.path.join(DEFAULT_OUT, "benchmarks")
    os.makedirs(outdir, exist_ok=True)
    if args.which == "all":
        todo = list(BENCHMARKS)
    else:
        if args.which not in BENCHMARKS:
            print(f"unknown benchmark {args.which!r}; choose from "
                  f"{sorted(BENCHMARKS)} or 'all'")
            return 1
        todo = [args.which]
    failures = 0
    for name in todo:
        res = BENCHMARKS[name](outdir=os.path.join(outdir, name))
        rep = res["report"]
        rep.print()
        md, js = rep.save(os.path.join(outdir, name), "report")
        failures += rep.count(FAIL)
        if res.get("figure"):
            print(f"  figure: {res['figure']}")
        print(f"  report: {md}")
    print(f"\nbenchmarks finished: {'PASS' if failures == 0 else f'{failures} FAIL'}")
    return 1 if failures else 0


# ------------------------------------------------------------------ _kchunk worker
def cmd_kchunk(args):
    """Hidden worker: compute one k-point chunk and write it to --out."""
    from .workflows import parallel as _par
    job = (args.kind, np.load(args.kpts))
    res = _par._worker((job, args.cfg))
    np.savez(args.out, **res)
    return 0


# ------------------------------------------------------------------ status
STATUS = [
    ("core/wannier.py", "Wannier90 win/hr/centres parsing, H(k)", "stage A, benchmarked"),
    ("core/spin_operator.py", "VASP binary spn decoding + quality metrics", "stage A, structural checks (see report)"),
    ("core/spin_wannier.py", "band -> Wannier spin rotation (chk) + interpolation", "stage A, validated on PbNV (H_W-H_hr 9e-6 eV)"),
    ("observables/spin_texture.py", "static spin texture", "stage A, benchmarked on analytic models"),
    ("observables/pst.py", "PST collinearity metric", "stage A, benchmarked (Rashba/RD)"),
    ("drives/peierls.py", "general vector potential + Peierls H(k,t)", "stage A, internal consistency only"),
    ("solvers/evolution.py", "U(T) quasienergies, Ebar, photon weights", "stage B, benchmarked (stage-B report)"),
    ("solvers/sambe.py", "Sambe Floquet Hamiltonian cross-check", "stage B, benchmarked vs evolution"),
    ("solvers/effective.py", "high-frequency van Vleck expansion", "stage B, benchmarked (convergence)"),
    ("observables/spectral.py", "Floquet spectral function A(k,E)", "stage B, implemented"),
    ("observables/floquet_spin.py", "periodic-averaged spin + spin-resolved spectra", "stage B, benchmarked (Rashba)"),
    ("observables/soc_splitting.py", "SOC splitting / effective-term fitting", "stage C, benchmarked (Rashba)"),
    ("observables/pst.py", "PST collinearity / axis / rotation analysis", "stage C, benchmarked (Rashba/RD)"),
    ("observables/rashba_fit.py", "linear Rashba fit of gap-edge pairs (alpha_R)", "ported from user script; model-level"),
    ("workflows/scan.py", "drive-parameter scans (hw/gamma)", "stage D, implemented + benchmarked (Rashba)"),
    ("workflows/converge.py", "n_step / m_sp convergence", "stage D, implemented + benchmarked (Rashba/PbNV)"),
    ("benchmarks/pt.py", "Benchmark 7: fcc Pt (18 Wannier, chk+eig input)", "stage D, PASS/WARN report"),
    ("benchmarks_data/Pt/", "public tutorial data for Benchmark 7 (provenance in README)", "data, 6.7 MB"),
]


def cmd_status(args):
    print(f"{'module':34s} {'feature':44s} status")
    print("-" * 110)
    for mod, feat, st in STATUS:
        print(f"{mod:34s} {feat:44s} {st}")
    return 0


# ------------------------------------------------------------------ main
def build_parser():
    p = argparse.ArgumentParser(
        prog="floquet.py",
        description="Wannier/Floquet/SOC engineering framework (stage A)")
    p.add_argument("--version", action="version",
                   version=f"floquet_soc {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    v = sub.add_parser("validate", help="parse inputs and run numerical checks")
    v.add_argument("--config", required=True)
    v.add_argument("--outdir", default=None)
    v.set_defaults(func=cmd_validate)

    b = sub.add_parser("bands", help="static bands along a k path")
    b.add_argument("--config", required=True)
    b.add_argument("--nk", type=int, default=None)
    b.add_argument("--outdir", default=None)
    b.set_defaults(func=cmd_bands)

    s = sub.add_parser("spin_static", help="static spin texture")
    s.add_argument("--config", required=True)
    s.add_argument("--kmesh-n", type=int, default=None,
                   help="k-mesh points per direction (default from config)")
    s.add_argument("--kmax", type=float, default=None,
                   help="fractional |k_i| upper bound (0.55 ~ full 2D BZ)")
    s.add_argument("--jobs", type=int, default=None,
                   help="k-point worker processes (default: config 'jobs')")
    s.add_argument("--outdir", default=None)
    s.set_defaults(func=cmd_spin_static)

    stx = sub.add_parser("spin_texture",
                         help="generic static spin-texture map (full/window)")
    stx.add_argument("--config", required=True)
    stx.add_argument("--mode", choices=["full", "window"], default="full")
    stx.add_argument("--band", default=None,
                     help="band index or vbm/cbm (default: config/VBM)")
    stx.add_argument("--center", nargs=3, type=float, default=None,
                     help="window center in fractional coords (K=0.333 0.333 0)")
    stx.add_argument("--half", type=float, default=0.05,
                     help="window half-width in 1/Angstrom")
    stx.add_argument("--n", type=int, default=41,
                     help="points per direction")
    stx.add_argument("--kmax", type=float, default=1.65,
                     help="full mode: Cartesian |k| bound in 1/Angstrom")
    stx.add_argument("--mask-split", type=float, default=0.0,
                     help="suppress in-plane arrows where the neighbour gap "
                          "is below this value (eV)")
    stx.add_argument("--inplane-ratio", type=float, default=None,
                     help="|S_xy|/|S_z| below this cutoff -> out-of-plane "
                          "(Sx=Sy=0, dot); default: config spin_inplane_ratio "
                          "or 0.05")
    stx.add_argument("--chunk", type=int, default=64)
    stx.add_argument("--figsize", nargs=2, type=float, default=None)
    stx.add_argument("--outdir", default=None)
    stx.set_defaults(func=cmd_spin_texture)

    dsc = sub.add_parser("dft-spin-check",
                         help="compare a DFT spin-texture table with the model")
    dsc.add_argument("--config", required=True)
    dsc.add_argument("--dft", required=True, help="DFT spin texture .dat")
    dsc.add_argument("--band", default=None,
                     help="model band index / vbm / cbm (default: band_dft-1)")
    dsc.add_argument("--mask-split", type=float, default=0.05,
                     help="exclude points with neighbour gap below this (eV)")
    dsc.add_argument("--prefix", default="dft_spin_compare")
    dsc.add_argument("--outdir", default=None)
    dsc.set_defaults(func=cmd_dft_spin_check)

    pv = sub.add_parser("provenance",
                        help="record version/config/input hashes/settings")
    pv.add_argument("--config", required=True)
    pv.add_argument("--outdir", default=None)
    pv.set_defaults(func=cmd_provenance)

    ed = sub.add_parser("export_dressed",
                        help="write light-dressed Wannier H_p(R) as _hr.dat "
                             "files for downstream calculations")
    ed.add_argument("--config", required=True)
    ed.add_argument("--pmax", type=int, default=None,
                    help="photon order (default: config numerics.pmax or 3)")
    ed.add_argument("--outdir", default=None)
    ed.set_defaults(func=cmd_export_dressed)

    fp = sub.add_parser("floquet_points",
                        help="band-resolved Floquet spin scan at k points")
    fp.add_argument("--config", required=True)
    fp.add_argument("--points", default=None,
                    help="comma labels (K,K′,Γ) or semicolon-separated "
                         "fractional triples")
    fp.add_argument("--hw", type=float, default=None)
    fp.add_argument("--gamma", type=float, nargs="+", default=None)
    fp.add_argument("--pol", nargs="+", default=None,
                    help="circ+, circ-, linear")
    fp.add_argument("--n-step", type=int, default=150)
    fp.add_argument("--m-sp", type=int, default=16)
    fp.add_argument("--rashba", action="store_true",
                    help="also scan the Gamma->K effective Rashba slope")
    fp.add_argument("--rashba-tmax", type=float, default=0.03)
    fp.add_argument("--rashba-n", type=int, default=7)
    fp.add_argument("--grid", action="store_true",
                    help="also compute small valley grids")
    fp.add_argument("--grid-center", default="K",
                    help="K, Kp, Gamma or a triple")
    fp.add_argument("--grid-half", type=float, default=0.01)
    fp.add_argument("--grid-n", type=int, default=5)
    fp.add_argument("--grid-gamma", type=float, default=0.05)
    fp.add_argument("--grid-pol", default="circ+")
    fp.add_argument("--subspace", action="store_true",
                    help="also report degenerate-subspace averaged spin")
    fp.add_argument("--subspace-tol", type=float, default=5e-3,
                    help="energy tolerance for subspace grouping (eV)")
    fp.add_argument("--outdir", default=None)
    fp.set_defaults(func=cmd_floquet_points)

    fb = sub.add_parser("floquet_bands",
                        help="static vs Floquet quasienergy/Ebar along k path")
    fb.add_argument("--config", required=True)
    fb.add_argument("--nk-seg", type=int, default=None)
    fb.add_argument("--hw", type=float, default=None)
    fb.add_argument("--gamma", type=float, default=None)
    fb.add_argument("--pol", default=None,
                    help="circ+, circ-, linear, elliptical (default: config)")
    fb.add_argument("--n-step", type=int, default=100)
    fb.add_argument("--m-sp", type=int, default=16)
    fb.add_argument("--outdir", default=None)
    fb.set_defaults(func=cmd_floquet_bands)

    ss = sub.add_parser("soc_scan",
                        help="standard scan: alpha_R, Delta_spin, PST vs gamma")
    ss.add_argument("--config", required=True)
    ss.add_argument("--points", default="K,Kp,Gamma")
    ss.add_argument("--hw", type=float, default=None)
    ss.add_argument("--gamma", type=float, nargs="+", default=None)
    ss.add_argument("--pol", nargs="+", default=None)
    ss.add_argument("--pst-n", type=int, default=7)
    ss.add_argument("--pst-kmax", type=float, default=0.35)
    ss.add_argument("--n-step", type=int, default=100)
    ss.add_argument("--m-sp", type=int, default=16)
    ss.add_argument("--outdir", default=None)
    ss.set_defaults(func=cmd_soc_scan)

    bm = sub.add_parser("benchmark", help="analytic model benchmarks")
    bm.add_argument("--which", default="all")
    bm.add_argument("--outdir", default=None)
    bm.set_defaults(func=cmd_benchmark)

    kc = sub.add_parser("_kchunk", help=argparse.SUPPRESS)
    kc.add_argument("--cfg", required=True)
    kc.add_argument("--kind", required=True)
    kc.add_argument("--kpts", required=True)
    kc.add_argument("--out", required=True)
    kc.set_defaults(func=cmd_kchunk)

    st = sub.add_parser("status", help="implementation status table")
    st.set_defaults(func=cmd_status)

    sp = sub.add_parser("spectrum", help="Floquet spectral function A(k,E)")
    sp.add_argument("--config", required=True)
    sp.add_argument("--nk", type=int, default=None)
    sp.add_argument("--kmax", type=float, default=0.3)
    sp.add_argument("--n-step", type=int, default=None)
    sp.add_argument("--emin", type=float, default=None)
    sp.add_argument("--emax", type=float, default=None)
    sp.add_argument("--outdir", default=None)
    sp.set_defaults(func=cmd_spectrum)

    sf = sub.add_parser("spin_floquet", help="Floquet spin texture")
    sf.add_argument("--config", required=True)
    sf.add_argument("--kmesh-n", type=int, default=None,
                    help="k-mesh points per direction (default from config)")
    sf.add_argument("--kmax", type=float, default=None)
    sf.add_argument("--jobs", type=int, default=None,
                    help="k-point worker processes (default: config 'jobs')")
    sf.add_argument("--outdir", default=None)
    sf.set_defaults(func=cmd_spin_floquet)

    soc = sub.add_parser("soc", help="SOC-related band splitting")
    soc.add_argument("--config", required=True)
    soc.add_argument("--nk", type=int, default=None)
    soc.add_argument("--kmax", type=float, default=0.3)
    soc.add_argument("--outdir", default=None)
    soc.set_defaults(func=cmd_soc)

    pst = sub.add_parser("pst", help="persistent spin texture analysis")
    pst.add_argument("--config", required=True)
    pst.add_argument("--kmesh-n", type=int, default=None)
    pst.add_argument("--kmax", type=float, default=None)
    pst.add_argument("--jobs", type=int, default=None,
                     help="k-point worker processes (default: config 'jobs')")
    pst.add_argument("--outdir", default=None)
    pst.set_defaults(func=cmd_pst)

    rf = sub.add_parser("rashba",
                        help="linear Rashba fit of the gap-edge band pairs")
    rf.add_argument("--config", required=True)
    rf.add_argument("--outdir", default=None)
    rf.set_defaults(func=cmd_rashba)

    scan = sub.add_parser("scan", help="drive-parameter scan (stage D)")
    scan.add_argument("--config", required=True)
    scan.add_argument("--hw", type=float, nargs="+", default=None)
    scan.add_argument("--gamma", type=float, nargs="+", default=None)
    scan.add_argument("--n-step", type=int, default=None)
    scan.add_argument("--m-sp", type=int, default=None)
    scan.add_argument("--outdir", default=None)
    scan.set_defaults(func=cmd_scan)

    conv = sub.add_parser("converge", help="solver convergence (stage D)")
    conv.add_argument("--config", required=True)
    conv.add_argument("--n-step", type=int, nargs="+", default=None)
    conv.add_argument("--m-sp", type=int, nargs="+", default=None)
    conv.add_argument("--outdir", default=None)
    conv.set_defaults(func=cmd_converge)

    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())


