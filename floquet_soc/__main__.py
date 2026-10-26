# -*- coding: utf-8 -*-
"""``python -m floquet_soc`` entry point (mirrors ``floquet.py``)."""
import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
