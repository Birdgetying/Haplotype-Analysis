#!/usr/bin/env python3
"""Prepare the MaizeGo paper-genotype qHKW1/ZmBAM1d positive control."""

import argparse
import csv
import io
import json
import sys
import zipfile
from pathlib import Path
from typing import Optional, Sequence

if sys.platform == 'win32':
    try:
        if not isinstance(sys.stdout, io.TextIOWrapper) or sys.stdout.encoding != 'utf-8':
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        if not isinstance(sys.stderr, io.TextIOWrapper) or sys.stderr.encoding != 'utf-8':
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    except (ValueError, AttributeError):
        pass

from star_gene_data import (
    build_database_from_marker_matrix,
    extract_maizego_marker_matrix,
    scan_maizego_sv_catalogue_candidates,
    scan_maizego_sv_candidates,
)


DEFAULT_SV_ROOT = Path("external_data/maize_natgenet_2019/maizego")
DEFAULT_PHENOTYPE = DEFAULT_SV_ROOT / "blup_traits_final.csv"
DEFAULT_OUTPUT_ROOT = Path("star_gene_database/maize_natgenet_2019")
DEFAULT_CANDIDATE_TSV = DEFAULT_SV_ROOT / "qHKW1_paper_8p9kb_candidates.tsv"
DEFAULT_CATALOGUE_TSV = DEFAULT_SV_ROOT / "qHKW1_paper_8p9kb_catalogue_candidates.tsv"
DEFAULT_MARKER_MATRIX = Path("external_data/maize_natgenet_2019/qHKW1_exact_indel_variant_matrix.tsv")
DEFAULT_NORMALIZED_PHENOTYPE = Path("external_data/maize_natgenet_2019/qHKW1_100grainweight_phenotype.tsv")

QHKW1_CHROM = "1"
QHKW1_FIG_WINDOW_START = 30_440_000
QHKW1_FIG_WINDOW_END = 30_540_000
QHKW1_LENGTH_MIN = 8_500
QHKW1_LENGTH_MAX = 9_500


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Find the paper qHKW1/ZmBAM1d 8.9-kb indel in the full MaizeGo "
            "SV.386014 genotype package and build the star-gene database only "
            "after an exact candidate marker is present."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--sv-root", default=str(DEFAULT_SV_ROOT),
                        help="Directory containing SV.386014.zip or its extracted files")
    parser.add_argument("--sv-zip", default=None,
                        help="Explicit SV.386014.zip path; defaults to <sv-root>/SV.386014.zip")
    parser.add_argument("--phenotype-table", default=str(DEFAULT_PHENOTYPE),
                        help="MaizeGo BLUP phenotype table")
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT),
                        help="Output star-gene database root")
    parser.add_argument("--candidate-output", default=str(DEFAULT_CANDIDATE_TSV),
                        help="TSV file to write candidate marker scan results")
    parser.add_argument("--catalogue-output", default=str(DEFAULT_CATALOGUE_TSV),
                        help="TSV file to write no-sample SV catalogue candidates")
    parser.add_argument("--marker-output", default=str(DEFAULT_MARKER_MATRIX),
                        help="Sample-by-marker matrix written when exactly one candidate is found")
    parser.add_argument("--normalized-phenotype-output", default=str(DEFAULT_NORMALIZED_PHENOTYPE),
                        help="Standardized SampleID phenotype table written for database construction")
    parser.add_argument("--target-id", default="qHKW1", help="Target id for the database")
    parser.add_argument("--trait-column", default="100grainweight", help="Phenotype column for HKW")
    parser.add_argument("--phenotype-missing-value", action="append", default=["-999"],
                        help="Phenotype value to exclude before building the database; repeatable")
    parser.add_argument("--chrom", default=QHKW1_CHROM, help="Chromosome for qHKW1 scan")
    parser.add_argument("--window-start", type=int, default=QHKW1_FIG_WINDOW_START,
                        help="Paper figure qHKW1 window start")
    parser.add_argument("--window-end", type=int, default=QHKW1_FIG_WINDOW_END,
                        help="Paper figure qHKW1 window end")
    parser.add_argument("--length-min", type=int, default=QHKW1_LENGTH_MIN,
                        help="Minimum SV length for the paper 8.9-kb indel")
    parser.add_argument("--length-max", type=int, default=QHKW1_LENGTH_MAX,
                        help="Maximum SV length for the paper 8.9-kb indel")
    parser.add_argument("--force-extract", action="store_true",
                        help="Extract SV.386014.zip even if extracted files already exist")
    return parser


def find_candidate_matrices(sv_root: Path) -> Sequence[Path]:
    patterns = ("*.hmp.txt", "*.org.txt", "*.txt")
    paths = []
    for pattern in patterns:
        paths.extend(sv_root.rglob(pattern))
    unique = []
    seen = set()
    for path in paths:
        normalized = str(path.resolve())
        if normalized in seen:
            continue
        seen.add(normalized)
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                header = f.readline().rstrip("\n\r").split("\t")
        except (IndexError, OSError):
            continue
        if {"rs#", "chrom", "pos"}.issubset(set(header)):
            unique.append(path)
    return unique


def find_sv_catalogues(sv_root: Path) -> Sequence[Path]:
    patterns = ("svs.final.*.txt",)
    paths = []
    for pattern in patterns:
        paths.extend(sv_root.rglob(pattern))
    unique = []
    seen = set()
    for path in paths:
        normalized = str(path.resolve())
        if normalized in seen:
            continue
        seen.add(normalized)
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as f:
                first_line = f.readline().rstrip("\n\r").split("\t")
        except OSError:
            continue
        if len(first_line) >= 5 and first_line[3].lower() in {"insertion", "deletion"}:
            unique.append(path)
    return unique


def _full_package_extracted(sv_root: Path) -> bool:
    return bool(find_sv_catalogues(sv_root))


def ensure_extracted(sv_root: Path, sv_zip: Path, force_extract: bool = False) -> int:
    existing = find_candidate_matrices(sv_root)
    existing_catalogues = find_sv_catalogues(sv_root)
    if existing_catalogues and not force_extract:
        return len(existing)
    if existing and not force_extract and (not sv_zip.exists() or _full_package_extracted(sv_root)):
        if not sv_zip.exists():
            print(f"[WARN] Full paper genotype package not found: {sv_zip}")
            print("[WARN] Scanning currently available MaizeGo matrices only; qHKW1 remains data-blocked if no exact marker is found.")
        return len(existing)
    if not sv_zip.exists():
        print("[BLOCKED] paper genotype package is missing")
        print(f"Expected file: {sv_zip}")
        print("Manual source: https://pan.baidu.com/s/10ieQpWGTEC805K4sI4RHOg")
        print("After downloading, place SV.386014.zip at the expected path and rerun this script.")
        return -1
    sv_root.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(sv_zip) as zf:
        zf.extractall(sv_root / "SV.386014")
    return len(find_candidate_matrices(sv_root))


def write_candidates(path: Path, candidates: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "source_path", "line_number", "marker", "alleles", "chrom", "pos",
        "marker_start", "marker_end", "variant_type", "sv_length",
        "valid_sample_count", "state_count", "counts",
    ]
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in candidates:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_catalogue_candidates(path: Path, candidates: Sequence[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "source_path", "line_number", "record_id", "raw_chrom", "chrom",
        "start", "end", "variant_type", "sv_length", "strand", "anchor_id",
        "anchor_start", "anchor_end", "score", "sequence_length",
        "sample_genotype_columns", "has_sample_genotypes",
    ]
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in candidates:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def write_normalized_phenotype_table(
    phenotype_table: Path,
    output_path: Path,
    trait_column: str,
    source_sample_column: str = "<Trait>",
) -> Path:
    with open(phenotype_table, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or source_sample_column not in reader.fieldnames:
            raise ValueError(f"phenotype table missing sample column: {source_sample_column}")
        if trait_column not in reader.fieldnames:
            raise ValueError(f"phenotype table missing trait column: {trait_column}")
        rows = [
            {"SampleID": row[source_sample_column], trait_column: row[trait_column]}
            for row in reader
        ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["SampleID", trait_column],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    sv_root = Path(args.sv_root)
    sv_zip = Path(args.sv_zip) if args.sv_zip else sv_root / "SV.386014.zip"
    phenotype_table = Path(args.phenotype_table)
    candidate_output = Path(args.candidate_output)
    catalogue_output = Path(args.catalogue_output)
    marker_output = Path(args.marker_output)
    normalized_phenotype_output = Path(args.normalized_phenotype_output)
    output_root = Path(args.output_root)

    extracted_count = ensure_extracted(sv_root, sv_zip, force_extract=args.force_extract)
    if extracted_count < 0:
        return 2

    matrix_paths = find_candidate_matrices(sv_root)
    print(f"[INFO] Candidate matrix files: {len(matrix_paths)}")
    catalogue_paths = find_sv_catalogues(sv_root)
    print(f"[INFO] SV catalogue files without sample genotypes: {len(catalogue_paths)}")
    candidates = scan_maizego_sv_candidates(
        matrix_paths=matrix_paths,
        chrom=args.chrom,
        window_start=args.window_start,
        window_end=args.window_end,
        length_min=args.length_min,
        length_max=args.length_max,
        variant_types=["insertion", "deletion"],
    )
    write_candidates(candidate_output, candidates)
    print(f"[INFO] qHKW1 paper-window 8.5-9.5 kb candidates: {len(candidates)}")
    print(f"[INFO] Candidate TSV: {candidate_output}")

    catalogue_candidates = scan_maizego_sv_catalogue_candidates(
        catalogue_paths=catalogue_paths,
        chrom=args.chrom,
        window_start=args.window_start,
        window_end=args.window_end,
        length_min=args.length_min,
        length_max=args.length_max,
        variant_types=["insertion", "deletion"],
    )
    write_catalogue_candidates(catalogue_output, catalogue_candidates)
    print(f"[INFO] SV catalogue candidate records: {len(catalogue_candidates)}")
    print(f"[INFO] Catalogue TSV: {catalogue_output}")

    if not candidates:
        print("[BLOCKED] no exact qHKW1 8.9-kb candidate marker found in available SV matrices")
        if catalogue_candidates:
            print("[INFO] Matching/nearby SV catalogue rows were found, but these files are not sample-level genotype matrices.")
        print("This means the currently available local matrices are not sufficient for the paper-genotype positive control.")
        return 3
    if len(candidates) > 1:
        print("[BLOCKED] multiple qHKW1 8.9-kb candidate markers found; inspect the candidate TSV before choosing one")
        return 4
    if not phenotype_table.exists():
        print("[BLOCKED] phenotype table is missing")
        print(f"Expected file: {phenotype_table}")
        return 5

    candidate = candidates[0]
    marker_id = str(candidate["marker"])
    source_path = Path(str(candidate["source_path"]))
    extract_maizego_marker_matrix(
        matrix_path=source_path,
        marker_id=marker_id,
        output_path=marker_output,
    )
    write_normalized_phenotype_table(
        phenotype_table=phenotype_table,
        output_path=normalized_phenotype_output,
        trait_column=args.trait_column,
    )
    db_dir = build_database_from_marker_matrix(
        marker_matrix=marker_output,
        phenotype_table=normalized_phenotype_output,
        output_root=output_root,
        target_id=args.target_id,
        chrom=args.chrom,
        start=int(candidate["marker_start"]),
        end=int(candidate["marker_end"]),
        phenotype_columns=[args.trait_column],
        marker_columns=[marker_id],
        marker_positions={marker_id: int(candidate["pos"] or candidate["marker_start"])},
        expected_direction="increases_trait",
        sample_column="SampleID",
        min_haplotype_count=1,
        phenotype_missing_values=args.phenotype_missing_value,
    )
    gene_info_path = db_dir / "gene_info.json"
    with open(gene_info_path, "r", encoding="utf-8") as f:
        gene_info = json.load(f)
    gene_info["source"] = "maizego_paper_marker"
    gene_info["source_marker"] = marker_id
    gene_info["source_matrix"] = str(source_path)
    with open(gene_info_path, "w", encoding="utf-8") as f:
        json.dump(gene_info, f, ensure_ascii=False, indent=2)
    print(f"[INFO] Marker matrix: {marker_output}")
    print(f"[INFO] Normalized phenotype table: {normalized_phenotype_output}")
    print(f"[INFO] Built star-gene database: {db_dir}")
    print("[NEXT] python run_star_gene_validation.py --run-analysis --paper maize2019 --target qHKW1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
