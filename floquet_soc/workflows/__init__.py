# -*- coding: utf-8 -*-
"""Workflows (stage D; placeholders for parameter scans and convergence)."""
from __future__ import annotations

STAGE = "D"


def _not_yet(what):
    raise NotImplementedError(
        f"{what} belongs to stage {STAGE} (scan/convergence workflows) and is "
        f"not implemented yet.")


def parameter_scan(*a, **k):                 # pragma: no cover
    _not_yet("parameter scan")


def convergence(*a, **k):                    # pragma: no cover
    _not_yet("convergence study")
