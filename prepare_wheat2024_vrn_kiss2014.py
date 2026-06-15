#!/usr/bin/env python3
"""Prepare a VRN diagnostic-marker positive-control database from Kiss 2014."""

import argparse
import io
import json
import sys
from pathlib import Path
from typing import Optional, Sequence

LOCAL_PYDEPS = Path(__file__).resolve().parent / "external_data" / "_pydeps"
if LOCAL_PYDEPS.exists():
    sys.path.insert(0, str(LOCAL_PYDEPS))

if sys.platform == 'win32':
    try:
        if not isinstance(sys.stdout, io.TextIOWrapper) or sys.stdout.encoding != 'utf-8':
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        if not isinstance(sys.stderr, io.TextIOWrapper) or sys.stderr.encoding != 'utf-8':
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    except (ValueError, AttributeError):
        pass

import pandas as pd

from star_gene_data import build_database_from_marker_matrix


DEFAULT_SOURCE_WORKBOOK = Path("external_data/literature/vrn_kiss2014/Kiss2014_embedded_workbook.xls")
DEFAULT_OUTPUT_ROOT = Path("star_gene_database/wheat_nature_2024")
DEFAULT_INTERMEDIATE_ROOT = Path("external_data/wheat_nature_2024/vrn_kiss2014")

TARGET_ID = "VRN-Kiss2014"
DEV49_SHEET = "stepwise_reg2011_2012a"
DEV59_SHEET = "stepwise_reg2011_2012b"
VRN_MARKERS = ["VRN-A1", "VRN-B1", "VRN-D1"]
PPD_MARKERS = ["PPD-D1", "PPDB1"]
ALL_DIAGNOSTIC_MARKERS = [*VRN_MARKERS, *PPD_MARKERS]
PHENOTYPE_COLUMNS = ["DEV49_2011", "DEV49_2012", "DEV49_mean", "DEV59_2011", "DEV59_2012", "DEV59_mean"]

# These are marker-order pseudo-positions, not physical coordinates. The source
# is diagnostic PCR/CNV marker states rather than a VCF interval.
MARKER_POSITIONS = {
    "VRN-A1": 1,
    "VRN-B1": 2,
    "VRN-D1": 3,
}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert Kiss et al. 2014 supplementary VRN/PPD diagnostic marker "
            "workbook into the precomputed star-gene database format."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--source-workbook", default=str(DEFAULT_SOURCE_WORKBOOK),
                        help="Extracted Kiss 2014 embedded workbook (.xls or .xlsx)")
    parser.add_argument("--dev49-sheet", default=DEV49_SHEET,
                        help="Worksheet containing DEV49 2011/2012 and diagnostic markers")
    parser.add_argument("--dev59-sheet", default=DEV59_SHEET,
                        help="Worksheet containing DEV59 2011/2012 and diagnostic markers")
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT),
                        help="Output star-gene database root")
    parser.add_argument("--intermediate-root", default=str(DEFAULT_INTERMEDIATE_ROOT),
                        help="Directory for marker/phenotype TSV intermediates")
    parser.add_argument("--target-id", default=TARGET_ID,
                        help="Target id for the output database")
    parser.add_argument("--min-haplotype-count", type=int, default=1,
                        help="Minimum sample count for retaining a haplotype")
    parser.add_argument("--include-ppd-in-haplotype", action="store_true",
                        help="Use VRN plus PPD markers for haplotype definition")
    parser.add_argument("--single-marker-targets", action="store_true",
                        help="Also build one-marker targets for VRN-A1, VRN-B1, and VRN-D1")
    return parser


def _read_workbook_sheet(workbook: Path, sheet_name: str) -> pd.DataFrame:
    if not workbook.exists():
        raise FileNotFoundError(f"source workbook is missing: {workbook}")
    suffix = workbook.suffix.lower()
    engine = "xlrd" if suffix == ".xls" else None
    return pd.read_excel(workbook, sheet_name=sheet_name, engine=engine)


def _clean_marker_value(value: object) -> Optional[str]:
    if pd.isna(value):
        return None
    if isinstance(value, str):
        value = value.strip()
        if not value:
            return None
        if ":" in value:
            value = value.split(":", 1)[0].strip()
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value).strip()
    if pd.isna(number):
        return None
    if number.is_integer():
        return str(int(number))
    return str(number)


def _require_columns(df: pd.DataFrame, columns: Sequence[str], sheet_name: str) -> None:
    missing = [col for col in columns if col not in df.columns]
    if missing:
        raise ValueError(f"{sheet_name} missing columns: {', '.join(missing)}")


def load_kiss2014_table(
    workbook: Path,
    dev49_sheet: str = DEV49_SHEET,
    dev59_sheet: str = DEV59_SHEET,
) -> pd.DataFrame:
    dev49 = _read_workbook_sheet(workbook, dev49_sheet)
    dev59 = _read_workbook_sheet(workbook, dev59_sheet)

    _require_columns(
        dev49,
        ["Pedigre", *ALL_DIAGNOSTIC_MARKERS, "DEV49_2011", "DEV49_2012"],
        dev49_sheet,
    )
    _require_columns(
        dev59,
        ["Pedigre", *ALL_DIAGNOSTIC_MARKERS, "DEV59_2011", "DEV59_2012"],
        dev59_sheet,
    )

    dev49 = dev49[["Pedigre", *ALL_DIAGNOSTIC_MARKERS, "DEV49_2011", "DEV49_2012"]].copy()
    dev59 = dev59[["Pedigre", "DEV59_2011", "DEV59_2012"]].copy()
    dev49["Pedigre"] = dev49["Pedigre"].astype(str).str.strip()
    dev59["Pedigre"] = dev59["Pedigre"].astype(str).str.strip()

    for col in ALL_DIAGNOSTIC_MARKERS:
        dev49[col] = dev49[col].map(_clean_marker_value)
    for col in ["DEV49_2011", "DEV49_2012"]:
        dev49[col] = pd.to_numeric(dev49[col], errors="coerce")
    for col in ["DEV59_2011", "DEV59_2012"]:
        dev59[col] = pd.to_numeric(dev59[col], errors="coerce")

    merged = dev49.merge(dev59, on="Pedigre", how="inner")
    merged = merged.rename(columns={"Pedigre": "SampleID"})
    merged = merged.dropna(subset=["SampleID", *VRN_MARKERS, "DEV49_2011", "DEV49_2012", "DEV59_2011", "DEV59_2012"])
    for col in ["DEV49_2011", "DEV49_2012", "DEV59_2011", "DEV59_2012"]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce")
    merged = merged.dropna(subset=["DEV49_2011", "DEV49_2012", "DEV59_2011", "DEV59_2012"])
    merged["DEV49_mean"] = merged[["DEV49_2011", "DEV49_2012"]].mean(axis=1)
    merged["DEV59_mean"] = merged[["DEV59_2011", "DEV59_2012"]].mean(axis=1)

    if merged.empty:
        raise ValueError("no complete Kiss 2014 rows after filtering")
    if merged["SampleID"].duplicated().any():
        duplicated = merged.loc[merged["SampleID"].duplicated(), "SampleID"].tolist()
        raise ValueError("duplicated Kiss 2014 sample IDs: " + ", ".join(duplicated[:10]))
    return merged.sort_values("SampleID")


def write_precomputed_inputs(
    table: pd.DataFrame,
    intermediate_root: Path,
    target_id: str,
    marker_columns: Sequence[str],
) -> tuple[Path, Path]:
    target_dir = intermediate_root / target_id
    target_dir.mkdir(parents=True, exist_ok=True)
    marker_output = target_dir / "marker_matrix.tsv"
    phenotype_output = target_dir / "phenotype.tsv"
    summary_output = target_dir / "diagnostic_marker_summary.tsv"

    table[["SampleID", *marker_columns]].to_csv(marker_output, sep="\t", index=False)
    table[["SampleID", *PHENOTYPE_COLUMNS]].to_csv(phenotype_output, sep="\t", index=False)

    summary_rows = []
    for marker in ALL_DIAGNOSTIC_MARKERS:
        counts = table[marker].value_counts(dropna=False).sort_index()
        for state, count in counts.items():
            summary_rows.append({"marker": marker, "state": state, "count": int(count)})
    pd.DataFrame(summary_rows).to_csv(summary_output, sep="\t", index=False)
    return marker_output, phenotype_output


def update_metadata(
    db_dir: Path,
    source_workbook: Path,
    marker_columns: Sequence[str],
    table: pd.DataFrame,
    target_id: str,
    gene_symbol: str = "VRN",
) -> None:
    gene_info_path = db_dir / "gene_info.json"
    with gene_info_path.open("r", encoding="utf-8") as f:
        gene_info = json.load(f)
    gene_info.update({
        "gene_id": target_id,
        "gene_symbol": gene_symbol,
        "chrom": "diagnostic_marker_panel",
        "start": 1,
        "end": len(marker_columns),
        "gene_start": 1,
        "gene_end": len(marker_columns),
        "source": "kiss2014_vrn_diagnostic_markers",
        "source_workbook": str(source_workbook),
        "source_note": (
            "Kiss et al. 2014 supplementary embedded workbook; diagnostic marker "
            "states for VRN-A1/VRN-B1/VRN-D1 with 2011/2012 DEV49 and DEV59 heading dates."
        ),
        "marker_panel": list(marker_columns),
        "available_diagnostic_markers": ALL_DIAGNOSTIC_MARKERS,
        "n_source_samples": int(len(table)),
    })
    with gene_info_path.open("w", encoding="utf-8") as f:
        json.dump(gene_info, f, ensure_ascii=False, indent=2)

    variant_path = db_dir / "variant_info.csv"
    variant_df = pd.read_csv(variant_path)
    for col in ["ref", "alt", "annotation", "marker_id"]:
        if col in variant_df.columns:
            variant_df[col] = variant_df[col].astype("object")
    for idx, row in variant_df.iterrows():
        marker_id = str(row["marker_id"])
        variant_df.at[idx, "position"] = int(MARKER_POSITIONS.get(marker_id, idx + 1))
        variant_df.at[idx, "ref"] = "0"
        variant_df.at[idx, "alt"] = "1"
        variant_df.at[idx, "len_diff"] = 0
        variant_df.at[idx, "is_sv"] = False
        variant_df.at[idx, "annotation"] = "diagnostic_marker"
        variant_df.at[idx, "marker_id"] = marker_id
        variant_df.at[idx, "marker_meaning"] = "0=winter/non-dominant state; 1=spring/dominant diagnostic state"
    variant_df.to_csv(variant_path, index=False)


def build_vrn_database(
    source_workbook: Path,
    output_root: Path,
    intermediate_root: Path,
    target_id: str,
    min_haplotype_count: int,
    include_ppd_in_haplotype: bool,
    dev49_sheet: str = DEV49_SHEET,
    dev59_sheet: str = DEV59_SHEET,
) -> Path:
    table = load_kiss2014_table(source_workbook, dev49_sheet=dev49_sheet, dev59_sheet=dev59_sheet)
    marker_columns = ALL_DIAGNOSTIC_MARKERS if include_ppd_in_haplotype else VRN_MARKERS
    db_dir = build_vrn_database_from_table(
        table=table,
        source_workbook=source_workbook,
        output_root=output_root,
        intermediate_root=intermediate_root,
        target_id=target_id,
        marker_columns=marker_columns,
        min_haplotype_count=min_haplotype_count,
        gene_symbol="VRN",
    )
    return db_dir


def build_vrn_database_from_table(
    table: pd.DataFrame,
    source_workbook: Path,
    output_root: Path,
    intermediate_root: Path,
    target_id: str,
    marker_columns: Sequence[str],
    min_haplotype_count: int,
    gene_symbol: str,
) -> Path:
    marker_output, phenotype_output = write_precomputed_inputs(
        table=table,
        intermediate_root=intermediate_root,
        target_id=target_id,
        marker_columns=marker_columns,
    )
    marker_positions = {marker: i + 1 for i, marker in enumerate(marker_columns)}
    db_dir = build_database_from_marker_matrix(
        marker_matrix=marker_output,
        phenotype_table=phenotype_output,
        output_root=output_root,
        target_id=target_id,
        chrom="diagnostic_marker_panel",
        start=1,
        end=len(marker_columns),
        phenotype_columns=PHENOTYPE_COLUMNS,
        marker_columns=marker_columns,
        marker_positions=marker_positions,
        expected_direction="unknown",
        sample_column="SampleID",
        min_haplotype_count=min_haplotype_count,
    )
    update_metadata(db_dir, source_workbook, marker_columns, table, target_id=target_id, gene_symbol=gene_symbol)
    return db_dir


def build_single_marker_databases(
    source_workbook: Path,
    output_root: Path,
    intermediate_root: Path,
    min_haplotype_count: int,
    dev49_sheet: str = DEV49_SHEET,
    dev59_sheet: str = DEV59_SHEET,
) -> list[Path]:
    table = load_kiss2014_table(source_workbook, dev49_sheet=dev49_sheet, dev59_sheet=dev59_sheet)
    db_dirs = []
    for marker in VRN_MARKERS:
        db_dirs.append(
            build_vrn_database_from_table(
                table=table,
                source_workbook=source_workbook,
                output_root=output_root,
                intermediate_root=intermediate_root,
                target_id=f"{marker}-Kiss2014",
                marker_columns=[marker],
                min_haplotype_count=min_haplotype_count,
                gene_symbol=marker,
            )
        )
    return db_dirs


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    source_workbook = Path(args.source_workbook)
    output_root = Path(args.output_root)
    intermediate_root = Path(args.intermediate_root)

    db_dir = build_vrn_database(
        source_workbook=source_workbook,
        output_root=output_root,
        intermediate_root=intermediate_root,
        target_id=args.target_id,
        min_haplotype_count=args.min_haplotype_count,
        include_ppd_in_haplotype=args.include_ppd_in_haplotype,
        dev49_sheet=args.dev49_sheet,
        dev59_sheet=args.dev59_sheet,
    )
    phenotype_rows = len(pd.read_csv(db_dir / "phenotype_data.csv"))
    variant_rows = len(pd.read_csv(db_dir / "variant_info.csv"))
    print(f"[INFO] Source workbook: {source_workbook}")
    print(f"[INFO] Built {args.target_id}: {variant_rows} markers, {phenotype_rows} samples -> {db_dir}")
    print(f"[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 --target {args.target_id}")
    if args.single_marker_targets:
        single_dirs = build_single_marker_databases(
            source_workbook=source_workbook,
            output_root=output_root,
            intermediate_root=intermediate_root,
            min_haplotype_count=args.min_haplotype_count,
            dev49_sheet=args.dev49_sheet,
            dev59_sheet=args.dev59_sheet,
        )
        for single_dir in single_dirs:
            print(f"[INFO] Built single-marker target: {single_dir.name} -> {single_dir}")
        print("[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 "
              "--target VRN-A1-Kiss2014 --target VRN-B1-Kiss2014 --target VRN-D1-Kiss2014")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
