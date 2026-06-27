#!/usr/bin/env python3
"""Prepare a VRN-B1 structural-variant positive control from Frontiers 2022."""

import argparse
import io
import json
import sys
from pathlib import Path
from typing import Optional, Sequence

if sys.platform == "win32":
    try:
        if not isinstance(sys.stdout, io.TextIOWrapper) or sys.stdout.encoding != "utf-8":
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
        if not isinstance(sys.stderr, io.TextIOWrapper) or sys.stderr.encoding != "utf-8":
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
    except (ValueError, AttributeError):
        pass

import pandas as pd

from star_gene_data import build_database_from_marker_matrix


DEFAULT_SOURCE_WORKBOOK = Path("external_data/literature/vrn_b1_structural/PMC9676936_Table_1.xlsx")
DEFAULT_OUTPUT_ROOT = Path("star_gene_database/wheat_nature_2024")
DEFAULT_INTERMEDIATE_ROOT = Path("external_data/wheat_nature_2024/vrn_b1_frontiers2022")
TARGET_ID = "VRN-B1-Frontiers2022"
SINGLE_MARKER_TARGET_ID = "VRN-B1a-Frontiers2022"

MARKER_COLUMNS = [
    "VRN-B1_insertion_838_duplication",
    "VRN-B1_deletion_6851",
    "VRN-B1_deletion_37",
]
MARKER_TO_SOURCE_COLUMN = {
    "VRN-B1_insertion_838_duplication": "Duplication(838bp)",
    "VRN-B1_deletion_6851": "Deletion(6851bp)",
    "VRN-B1_deletion_37": "Deletion(37bp)",
}
MARKER_POSITIONS = {
    "VRN-B1_insertion_838_duplication": 1,
    "VRN-B1_deletion_6851": 2,
    "VRN-B1_deletion_37": 3,
}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert Makhoul et al. 2022 Frontiers VRN1 supplementary tables "
            "into a VRN-B1 structural-variant star-gene database."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--source-workbook", default=str(DEFAULT_SOURCE_WORKBOOK))
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument("--intermediate-root", default=str(DEFAULT_INTERMEDIATE_ROOT))
    parser.add_argument("--target-id", default=TARGET_ID)
    parser.add_argument("--single-marker-target-id", default=SINGLE_MARKER_TARGET_ID)
    parser.add_argument("--skip-single-marker-target", action="store_true")
    parser.add_argument("--min-haplotype-count", type=int, default=1)
    return parser


def _read_table_with_header(workbook: Path, sheet_name: str, header_row: int) -> pd.DataFrame:
    raw = pd.read_excel(workbook, sheet_name=sheet_name, header=None)
    if raw.shape[0] <= header_row:
        raise ValueError(f"{sheet_name} does not contain header row {header_row}")
    columns = [str(value).strip() if not pd.isna(value) else "" for value in raw.iloc[header_row]]
    df = raw.iloc[header_row + 1:].copy()
    df.columns = columns
    df = df.dropna(how="all")
    return df


def _normal_sample(value: object) -> str:
    return str(value).strip()


def _normal_key(value: object) -> str:
    return _normal_sample(value).casefold()


def _clean_marker(value: object) -> str:
    if pd.isna(value):
        return "-"
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return "-"
    return text


def build_frontiers2022_marker_tables(
    workbook: Path,
    marker_columns: Optional[Sequence[str]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return marker and phenotype tables from Frontiers 2022 S1/S12/S13."""
    workbook = Path(workbook)
    if not workbook.exists():
        raise FileNotFoundError(f"Frontiers 2022 workbook is missing: {workbook}")

    s1 = _read_table_with_header(workbook, "S1", header_row=1)
    s12 = _read_table_with_header(workbook, "S12", header_row=2)
    # S13 is loaded to make malformed workbooks fail early; the marker target is
    # based on the structural calls in S12, not on paper haplotype labels.
    s13 = _read_table_with_header(workbook, "S13", header_row=1)

    required_s1 = {"Cultivar name", "Heading date"}
    required_s12 = {"Genotype name", *MARKER_TO_SOURCE_COLUMN.values()}
    required_s13 = {"Cultivar name", "Haplotype of VRN-B1"}
    missing = sorted(required_s1 - set(s1.columns))
    if missing:
        raise ValueError("S1 missing columns: " + ", ".join(missing))
    missing = sorted(required_s12 - set(s12.columns))
    if missing:
        raise ValueError("S12 missing columns: " + ", ".join(missing))
    missing = sorted(required_s13 - set(s13.columns))
    if missing:
        raise ValueError("S13 missing columns: " + ", ".join(missing))

    s12 = s12.copy()
    s12["_key"] = s12["Genotype name"].map(_normal_key)
    s1 = s1.copy()
    s1["_key"] = s1["Cultivar name"].map(_normal_key)
    merged = s1.merge(s12[["_key", *MARKER_TO_SOURCE_COLUMN.values()]], on="_key", how="left")
    merged["Heading_date"] = pd.to_numeric(merged["Heading date"], errors="coerce")
    merged = merged.dropna(subset=["Cultivar name", "Heading_date"]).copy()

    selected_markers = list(marker_columns or MARKER_COLUMNS)
    unknown_markers = [marker for marker in selected_markers if marker not in MARKER_TO_SOURCE_COLUMN]
    if unknown_markers:
        raise ValueError("unknown Frontiers 2022 VRN-B1 markers: " + ", ".join(unknown_markers))

    marker_df = pd.DataFrame({"SampleID": merged["Cultivar name"].map(_normal_sample)})
    for marker in selected_markers:
        source_col = MARKER_TO_SOURCE_COLUMN[marker]
        marker_df[marker] = merged[source_col].map(_clean_marker)

    phenotype_df = pd.DataFrame({
        "SampleID": merged["Cultivar name"].map(_normal_sample),
        "Heading_date": merged["Heading_date"],
    })
    return marker_df, phenotype_df


def write_intermediate_tables(
    marker_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    intermediate_root: Path,
    target_id: str,
    marker_columns: Sequence[str],
) -> tuple[Path, Path]:
    target_dir = intermediate_root / target_id
    target_dir.mkdir(parents=True, exist_ok=True)
    marker_path = target_dir / "marker_matrix.tsv"
    phenotype_path = target_dir / "phenotype.tsv"
    summary_path = target_dir / "marker_summary.tsv"

    marker_df.to_csv(marker_path, sep="\t", index=False)
    phenotype_df.to_csv(phenotype_path, sep="\t", index=False)

    summary_rows = []
    for marker in marker_columns:
        for state, count in marker_df[marker].value_counts(dropna=False).sort_index().items():
            summary_rows.append({"marker": marker, "state": state, "count": int(count)})
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)
    return marker_path, phenotype_path


def update_database_metadata(
    db_dir: Path,
    source_workbook: Path,
    marker_df: pd.DataFrame,
    target_id: str,
    marker_columns: Sequence[str],
) -> None:
    gene_info_path = db_dir / "gene_info.json"
    with gene_info_path.open("r", encoding="utf-8") as f:
        gene_info = json.load(f)
    gene_info.update({
        "gene_id": target_id,
        "gene_symbol": "VRN-B1",
        "chrom": "VRN-B1_structural_marker_panel",
        "start": 1,
        "end": len(marker_columns),
        "gene_start": 1,
        "gene_end": len(marker_columns),
        "source": "frontiers2022_vrn1_structural_variants",
        "source_workbook": str(source_workbook),
        "source_note": (
            "Makhoul et al. 2022 Frontiers Supplementary Table S12/S1; "
            "VRN-B1 intron structural variants with heading-date phenotypes."
        ),
        "marker_panel": list(marker_columns),
        "n_source_samples": int(len(marker_df)),
        "causal_validation_marker": "VRN-B1_deletion_6851",
        "causal_validation_allele": "Vrn-B1a",
    })
    with gene_info_path.open("w", encoding="utf-8") as f:
        json.dump(gene_info, f, ensure_ascii=False, indent=2)

    variant_path = db_dir / "variant_info.csv"
    variant_df = pd.read_csv(variant_path)
    for idx, row in variant_df.iterrows():
        marker_id = str(row["marker_id"])
        variant_df.at[idx, "position"] = MARKER_POSITIONS.get(marker_id, idx + 1)
        variant_df.at[idx, "annotation"] = "intron_structural_variant"
        variant_df.at[idx, "source_column"] = MARKER_TO_SOURCE_COLUMN.get(marker_id, "")
        variant_df.at[idx, "source"] = "frontiers2022_table_s12"
        if marker_id == "VRN-B1_deletion_6851":
            variant_df.at[idx, "literature_allele"] = "Vrn-B1a"
            variant_df.at[idx, "literature_note"] = "6,851 bp first-intron deletion in VRN-B1"
    variant_df.to_csv(variant_path, index=False)


def build_database(
    source_workbook: Path,
    output_root: Path,
    intermediate_root: Path,
    target_id: str,
    min_haplotype_count: int,
    marker_columns: Optional[Sequence[str]] = None,
) -> Path:
    selected_markers = list(marker_columns or MARKER_COLUMNS)
    marker_df, phenotype_df = build_frontiers2022_marker_tables(source_workbook, selected_markers)
    marker_path, phenotype_path = write_intermediate_tables(
        marker_df,
        phenotype_df,
        intermediate_root,
        target_id,
        selected_markers,
    )
    db_dir = build_database_from_marker_matrix(
        marker_matrix=marker_path,
        phenotype_table=phenotype_path,
        output_root=output_root,
        target_id=target_id,
        chrom="VRN-B1_structural_marker_panel",
        start=1,
        end=len(selected_markers),
        phenotype_columns=["Heading_date"],
        marker_columns=selected_markers,
        marker_positions={marker: MARKER_POSITIONS[marker] for marker in selected_markers},
        expected_direction="decreases_trait",
        min_haplotype_count=min_haplotype_count,
    )
    update_database_metadata(db_dir, source_workbook, marker_df, target_id, selected_markers)
    return db_dir


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    db_dir = build_database(
        source_workbook=Path(args.source_workbook),
        output_root=Path(args.output_root),
        intermediate_root=Path(args.intermediate_root),
        target_id=args.target_id,
        min_haplotype_count=args.min_haplotype_count,
    )
    print(f"[INFO] Built Frontiers 2022 VRN-B1 structural database: {db_dir}")
    if not args.skip_single_marker_target:
        single_db_dir = build_database(
            source_workbook=Path(args.source_workbook),
            output_root=Path(args.output_root),
            intermediate_root=Path(args.intermediate_root),
            target_id=args.single_marker_target_id,
            min_haplotype_count=args.min_haplotype_count,
            marker_columns=["VRN-B1_deletion_6851"],
        )
        print(f"[INFO] Built Frontiers 2022 Vrn-B1a single-marker database: {single_db_dir}")
    print("[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 "
          f"--target {args.target_id} --target {args.single_marker_target_id} --score-mode robust_discovery")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
