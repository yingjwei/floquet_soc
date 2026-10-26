#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Entry point: ``python floquet.py <command> ...``.

Stage A commands: validate, bands, spin_static, benchmark, status.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from floquet_soc.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
