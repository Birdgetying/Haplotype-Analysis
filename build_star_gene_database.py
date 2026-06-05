#!/usr/bin/env python3
"""Build a star-gene precomputed database from a small marker matrix."""

import argparse
import io
import sys
from pathlib import Path
from typing import Dict, Optional, Sequence

if sys.platform == 'win32':
    try:
        if not isinstance(sys.stdout, io.TextIOWrapper) or sys.stdout.encoding != 'utf-8':
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        if not isinstance(sys.stderr, io.TextIOWrapper) or sys.stderr.encoding != 'utf-8':
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    except (ValueError, AttributeError):
        pass

from star_gene_data import build_database_from_marker_matrix


def parse_marker_positions(values: Optional[Sequence[str]]) -> Dict[str, int]:
    positions: Dict[str, int] = {}
    for value in values or []:
        if "=" not in value:
            raise ValueError(f"marker position must use marker=position format: {value}")
        marker, pos = value.split("=", 1)
        marker = marker.strip()
        if not marker:
            raise ValueError(f"marker name is empty in marker-position: {value}")
        positions[marker] = int(pos)
    return positions


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert a small sample-by-marker matrix into the star-gene database format.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--marker-matrix", required=True, help="CSV/TSV marker matrix with one sample per row")
    parser.add_argument("--phenotype-table", required=True, help="CSV/TSV phenotype table")
    parser.add_argument("--output-root", required=True, help="Output database root")
    parser.add_argument("--target-id", required=True, help="Target gene/locus id, e.g. qHKW1")
    parser.add_argument("--chrom", required=True, help="Chromosome/contig name")
    parser.add_argument("--start", required=True, type=int, help="Target start coordinate")
    parser.add_argument("--end", required=True, type=int, help="Target end coordinate")
    parser.add_argument("--phenotype-column", action="append", required=True,
                        help="Phenotype column to retain; repeatable")
    parser.add_argument("--marker-column", action="append", default=[],
                        help="Marker column to use; repeatable. Defaults to all non-sample columns.")
    parser.add_argument("--marker-position", action="append", default=[],
                        help="Marker position as marker=position; repeatable")
    parser.add_argument("--sample-column", default="SampleID", help="Sample ID column in both tables")
    parser.add_argument("--expected-direction", default="unknown",
                        choices=["unknown", "increases_trait", "decreases_trait"],
                        help="Expected phenotype direction for the beneficial marker/haplotype")
    parser.add_argument("--min-haplotype-count", type=int, default=1,
                        help="Minimum sample count required to retain a haplotype")
    parser.add_argument("--phenotype-missing-value", action="append", default=[],
                        help="Phenotype value to exclude before building the database; repeatable")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    marker_positions = parse_marker_positions(args.marker_position)
    db_dir = build_database_from_marker_matrix(
        marker_matrix=Path(args.marker_matrix),
        phenotype_table=Path(args.phenotype_table),
        output_root=Path(args.output_root),
        target_id=args.target_id,
        chrom=args.chrom,
        start=args.start,
        end=args.end,
        phenotype_columns=args.phenotype_column,
        marker_columns=args.marker_column or None,
        marker_positions=marker_positions or None,
        expected_direction=args.expected_direction,
        sample_column=args.sample_column,
        min_haplotype_count=args.min_haplotype_count,
        phenotype_missing_values=args.phenotype_missing_value,
    )
    print(f"[INFO] Built star-gene database: {db_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
