#!/usr/bin/env python3
"""Prepare an external Rht positive-control database from Zanke et al. 2014."""

import argparse
import io
import json
import sys
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

import pandas as pd

from star_gene_data import build_database_from_marker_matrix


DEFAULT_SOURCE_WORKBOOK = Path("external_data/literature/rht_zanke2014/Table_S2_candidate_genes_phenotypes.xlsx")
DEFAULT_OUTPUT_ROOT = Path("star_gene_database/wheat_nature_2024")
DEFAULT_INTERMEDIATE_ROOT = Path("external_data/wheat_nature_2024/rht_zanke2014")
TARGET_ID = "Rht-Zanke2014"
SHEET_NAME = "Table_S2"
MARKER_COLUMNS = ["Rht-B1", "Rht-D1"]
MARKER_REF_ALT = {
    "Rht-B1": ("B1a", "B1b"),
    "Rht-D1": ("Rht-D1a", "D1b"),
}
MARKER_MEANING = {
    "Rht-B1": "B1b=canonical Rht-B1b semi-dwarf mutant marker; B1a=wild type",
    "Rht-D1": "D1b=canonical Rht-D1b semi-dwarf mutant marker; Rht-D1a=wild type",
}
PHENOTYPE_COLUMNS = [
    "PlantHeight_BLUE",
    "PlantHeight_GAT_2012",
    "PlantHeight_09_AND",
    "PlantHeight_09_SEL",
    "PlantHeight_09_WOH",
    "PlantHeight_10_AND",
    "PlantHeight_10_JAN",
    "PlantHeight_10_SAU",
    "PlantHeight_10_SEL",
    "PlantHeight_10_WOH",
]

RAW_COLUMN_MAP = {
    0: "VarietyName",
    1: "GW_no",
    2: "Habit",
    3: "Rht-B1b",
    4: "Rht-B1a",
    5: "Rht-D1b",
    6: "Rht-D1a",
    16: "PlantHeight_09_AND",
    17: "PlantHeight_09_SEL",
    18: "PlantHeight_09_WOH",
    19: "PlantHeight_10_AND",
    20: "PlantHeight_10_JAN",
    21: "PlantHeight_10_SAU",
    22: "PlantHeight_10_SEL",
    23: "PlantHeight_10_WOH",
    24: "PlantHeight_BLUE",
    25: "PlantHeight_GAT_2012",
}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert Zanke et al. 2014 PLOS ONE Table S2 Rht-B1/Rht-D1 marker "
            "states and plant-height phenotypes into the precomputed database format."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--source-workbook", default=str(DEFAULT_SOURCE_WORKBOOK),
                        help="PLOS ONE Table S2 workbook")
    parser.add_argument("--sheet-name", default=SHEET_NAME,
                        help="Worksheet containing Table S2")
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT),
                        help="Output star-gene database root")
    parser.add_argument("--intermediate-root", default=str(DEFAULT_INTERMEDIATE_ROOT),
                        help="Directory for marker/phenotype TSV intermediates")
    parser.add_argument("--target-id", default=TARGET_ID,
                        help="Target id for the output database")
    parser.add_argument("--min-haplotype-count", type=int, default=1,
                        help="Minimum sample count for retaining a haplotype")
    parser.add_argument("--single-marker-targets", action="store_true",
                        help="Also build one-marker targets for Rht-B1 and Rht-D1")
    return parser


def _marker_state(mutant_value: object, wild_type_value: object, mutant_label: str, wild_type_label: str) -> Optional[str]:
    mutant = pd.to_numeric(pd.Series([mutant_value]), errors="coerce").iloc[0]
    wild_type = pd.to_numeric(pd.Series([wild_type_value]), errors="coerce").iloc[0]
    has_mutant = pd.notna(mutant) and float(mutant) != 0
    has_wild_type = pd.notna(wild_type) and float(wild_type) != 0
    if has_mutant and not has_wild_type:
        return mutant_label
    if has_wild_type and not has_mutant:
        return wild_type_label
    if has_mutant and has_wild_type:
        return f"{mutant_label}/{wild_type_label}"
    return None


def load_zanke2014_table(source_workbook: Path, sheet_name: str = SHEET_NAME) -> pd.DataFrame:
    if not source_workbook.exists():
        raise FileNotFoundError(f"source workbook is missing: {source_workbook}")

    raw = pd.read_excel(source_workbook, sheet_name=sheet_name, header=None, dtype=object)
    data = raw.iloc[4:].copy()
    data = data.rename(columns=RAW_COLUMN_MAP)
    required = list(RAW_COLUMN_MAP.values())
    missing = [col for col in required if col not in data.columns]
    if missing:
        raise ValueError(f"{sheet_name} missing expected columns: {', '.join(missing)}")

    data["VarietyName"] = data["VarietyName"].astype(str).str.strip()
    data = data[~data["VarietyName"].isin(["", "nan", "NaN", "Sum", "n. d. = not detected"])]
    data["SampleID"] = data["GW_no"].astype(str).str.strip()
    data.loc[data["SampleID"].isin(["", "nan", "NaN"]), "SampleID"] = data["VarietyName"]

    data["Rht-B1"] = [
        _marker_state(mutant, wild_type, "B1b", "B1a")
        for mutant, wild_type in zip(data["Rht-B1b"], data["Rht-B1a"])
    ]
    data["Rht-D1"] = [
        _marker_state(mutant, wild_type, "D1b", "Rht-D1a")
        for mutant, wild_type in zip(data["Rht-D1b"], data["Rht-D1a"])
    ]

    for col in PHENOTYPE_COLUMNS:
        data[col] = pd.to_numeric(data[col], errors="coerce")

    keep_cols = ["SampleID", "VarietyName", "Habit", *MARKER_COLUMNS, *PHENOTYPE_COLUMNS]
    out = data[keep_cols].dropna(subset=["SampleID", *MARKER_COLUMNS, "PlantHeight_BLUE"]).copy()
    out = out[out["SampleID"].astype(str).str.len() > 0]
    out = out.groupby("SampleID", as_index=False).agg({
        "VarietyName": "first",
        "Habit": "first",
        "Rht-B1": "first",
        "Rht-D1": "first",
        **{col: "mean" for col in PHENOTYPE_COLUMNS},
    })
    if out.empty:
        raise ValueError("no complete Zanke 2014 Rht rows after filtering")
    return out.sort_values("SampleID")


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
    summary_output = target_dir / "marker_summary.tsv"

    table[["SampleID", *marker_columns]].to_csv(marker_output, sep="\t", index=False)
    table[["SampleID", *PHENOTYPE_COLUMNS]].to_csv(phenotype_output, sep="\t", index=False)

    rows = []
    for marker in marker_columns:
        for state, count in table[marker].value_counts(dropna=False).sort_index().items():
            rows.append({"marker": marker, "state": state, "count": int(count)})
    pd.DataFrame(rows).to_csv(summary_output, sep="\t", index=False)
    return marker_output, phenotype_output


def update_metadata(
    db_dir: Path,
    source_workbook: Path,
    table: pd.DataFrame,
    target_id: str,
    marker_columns: Sequence[str],
    gene_symbol: str = "Rht1",
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
        "source": "zanke2014_rht_candidate_gene_table_s2",
        "source_workbook": str(source_workbook),
        "source_note": (
            "Zanke et al. 2014 PLOS ONE Table S2: varieties, candidate-gene "
            "genotypes and multi-environment plant-height phenotypes."
        ),
        "marker_panel": list(marker_columns),
        "available_rht_markers": MARKER_COLUMNS,
        "n_source_samples": int(len(table)),
    })
    with gene_info_path.open("w", encoding="utf-8") as f:
        json.dump(gene_info, f, ensure_ascii=False, indent=2)

    variant_path = db_dir / "variant_info.csv"
    variant_df = pd.read_csv(variant_path)
    for col in ["ref", "alt", "annotation", "marker_id"]:
        if col in variant_df.columns:
            variant_df[col] = variant_df[col].astype("object")
    for idx, marker in enumerate(marker_columns):
        ref, alt = MARKER_REF_ALT[marker]
        variant_df.at[idx, "position"] = idx + 1
        variant_df.at[idx, "ref"] = ref
        variant_df.at[idx, "alt"] = alt
        variant_df.at[idx, "marker_id"] = marker
        variant_df.at[idx, "marker_meaning"] = MARKER_MEANING[marker]
    variant_df["len_diff"] = 0
    variant_df["is_sv"] = False
    variant_df["annotation"] = "functional_marker"
    variant_df.to_csv(variant_path, index=False)


def build_rht_zanke_database(
    source_workbook: Path,
    output_root: Path,
    intermediate_root: Path,
    target_id: str = TARGET_ID,
    min_haplotype_count: int = 1,
    sheet_name: str = SHEET_NAME,
) -> Path:
    table = load_zanke2014_table(source_workbook, sheet_name=sheet_name)
    return build_rht_zanke_database_from_table(
        table=table,
        source_workbook=source_workbook,
        output_root=output_root,
        intermediate_root=intermediate_root,
        target_id=target_id,
        marker_columns=MARKER_COLUMNS,
        min_haplotype_count=min_haplotype_count,
        gene_symbol="Rht1",
    )


def build_rht_zanke_database_from_table(
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
        marker_positions={marker: i + 1 for i, marker in enumerate(marker_columns)},
        expected_direction="decreases_trait",
        sample_column="SampleID",
        min_haplotype_count=min_haplotype_count,
    )
    update_metadata(
        db_dir=db_dir,
        source_workbook=source_workbook,
        table=table,
        target_id=target_id,
        marker_columns=marker_columns,
        gene_symbol=gene_symbol,
    )
    return db_dir


def build_single_marker_databases(
    source_workbook: Path,
    output_root: Path,
    intermediate_root: Path,
    min_haplotype_count: int,
    sheet_name: str = SHEET_NAME,
) -> list[Path]:
    table = load_zanke2014_table(source_workbook, sheet_name=sheet_name)
    db_dirs = []
    for marker in MARKER_COLUMNS:
        db_dirs.append(
            build_rht_zanke_database_from_table(
                table=table,
                source_workbook=source_workbook,
                output_root=output_root,
                intermediate_root=intermediate_root,
                target_id=f"{marker}-Zanke2014",
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

    db_dir = build_rht_zanke_database(
        source_workbook=source_workbook,
        output_root=output_root,
        intermediate_root=intermediate_root,
        target_id=args.target_id,
        min_haplotype_count=args.min_haplotype_count,
        sheet_name=args.sheet_name,
    )
    gene_info = json.loads((db_dir / "gene_info.json").read_text(encoding="utf-8"))
    print(f"[OK] Built {args.target_id}: {gene_info.get('n_source_samples')} samples")
    print(f"[OK] Database: {db_dir}")
    print(f"[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 --target {args.target_id} --score-mode robust_discovery")
    if args.single_marker_targets:
        single_dirs = build_single_marker_databases(
            source_workbook=source_workbook,
            output_root=output_root,
            intermediate_root=intermediate_root,
            min_haplotype_count=args.min_haplotype_count,
            sheet_name=args.sheet_name,
        )
        for single_dir in single_dirs:
            print(f"[OK] Built single-marker target: {single_dir.name} -> {single_dir}")
        print("[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 "
              "--target Rht-D1-Zanke2014 --score-mode robust_discovery")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
