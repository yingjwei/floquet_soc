# -*- coding: utf-8 -*-
"""Floquet-SOC framework version information."""
__version__ = "0.1.7-rc1"
__release__ = (
    "Benchmarked Wannier-Floquet-SOC research framework - validated "
    "numerical core, material-specific quantitative validation required. "
    "v0.1.7-rc1 makes static spin-texture maps and DFT spin-texture "
    "cross-checks generic framework commands; v0.1.6-rc1 is the Wannier90 "
    "use_ws_distance patch: wannier90_wsvec.dat "
    "pair-dependent supercell shifts are applied to H(k), the Wannier-basis "
    "spin operator and the Peierls-driven Hamiltonian (unless "
    "write_ndegen_applied is set); the previous plain-Fourier interpolation "
    "was wrong off-grid for MoSSe (Rashba alpha 64.8 -> 44.3 meV*A). "
    "v0.1.5 fixed the odd-mp_grid WS R set; v0.1.6 builds on it. "
    "Pt/PbNV regressions unchanged (no wsvec files there).")
