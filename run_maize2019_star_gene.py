#!/usr/bin/env python3
"""Check or run Maize Nature Genetics 2019 qHKW1/ZmBAM1d validation targets."""

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
    raise SystemExit(main(['--paper', 'maize2019', *sys.argv[1:]]))
