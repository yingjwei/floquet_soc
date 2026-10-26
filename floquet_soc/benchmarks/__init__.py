# -*- coding: utf-8 -*-
"""Analytic benchmark registry (stage A + stage B)."""
from __future__ import annotations

from . import (rashba, rashba_dresselhaus, floquet_rashba,
               pbnv_spin, pt, spin_degenerate, magnetic_soc)

BENCHMARKS = {
    "rashba": rashba.run,
    "rashba_dresselhaus": rashba_dresselhaus.run,
    "floquet_rashba": floquet_rashba.run,
    "pbnv_spin": pbnv_spin.run,
    "pt": pt.run,
    "spin_degenerate": spin_degenerate.run,
    "magnetic_soc": magnetic_soc.run,
}


def run_all(outdir=None):
    results = []
    for name, fn in BENCHMARKS.items():
        sub = None if outdir is None else f"{outdir}/{name}"
        results.append((name, fn(outdir=sub)))
    return results
