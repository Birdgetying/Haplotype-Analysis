#!/usr/bin/env python3
"""Check or run Wheat Nature 2024 WatSeq star-gene validation targets."""

import io
import sys

if sys.platform == 'win32':
    try:
        if not isinstance(sys.stdout, io.TextIOWrapper) or sys.stdout.encoding != 'utf-8':
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        if not isinstance(sys.stderr, io.TextIOWrapper) or sys.stderr.encoding != 'utf-8':
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    except (ValueError, AttributeError):
        pass

from star_gene_validation import main


if __name__ == '__main__':
    raise SystemExit(main(['--paper', 'wheat2024', *sys.argv[1:]]))
