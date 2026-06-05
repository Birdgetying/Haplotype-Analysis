#!/usr/bin/env python3
"""Prepare the Wheat Nature 2024 Q7B-PH Figure 3g positive control."""

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


DEFAULT_SOURCE_XLSX = Path(
    "external_data/wheat_nature_2024/wwwg2b/q7b_ph/"
    "Watseq_Figure_3g_NIL_Q7B-PH_field_data.xlsx"
)
DEFAULT_OUTPUT_ROOT = Path("star_gene_database/wheat_nature_2024")
DEFAULT_MARKER_OUTPUT = Path("external_data/wheat_nature_2024/q7b_ph_figure3g_marker_matrix.tsv")
DEFAULT_PHENOTYPE_OUTPUT = Path("external_data/wheat_nature_2024/q7b_ph_figure3g_phenotype.tsv")

TARGET_ID = "Q7B-HT"
MARKER_ID = "Q7B-PH_allele"
PHENOTYPE_COLUMN = "PH_M_cm"

# Figure 3d places the Q7B-PH peak on chromosome 7B near genetic position 67 cM.
# This is a paper-source haplotype positive control, not a physical-coordinate
# chr7B VCF extraction, so use an explicit pseudo-coordinate for the marker DB.
Q7B_PH_CHROM = "7B"
Q7B_PH_START = 67
Q7B_PH_END = 68


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert WWWG2B WatSeq Figure 3g Q7B-PH NIL field data into the "
            "precomputed star-gene database format."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--source-xlsx", default=str(DEFAULT_SOURCE_XLSX),
                        help="WWWG2B Figure 3g Q7B-PH field-data workbook")
    parser.add_argument("--sheet-name", default="Figure_3_NIL_Q7B-PH_W141_field_",
                        help="Worksheet containing Accession, allele, and PH_M_cm")
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT),
                        help="Output star-gene database root")
    parser.add_argument("--marker-output", default=str(DEFAULT_MARKER_OUTPUT),
                        help="Sample-by-marker TSV to write")
    parser.add_argument("--phenotype-output", default=str(DEFAULT_PHENOTYPE_OUTPUT),
                        help="Sample phenotype TSV to write")
    parser.add_argument("--target-id", default=TARGET_ID,
                        help="Target id for the database")
    parser.add_argument("--sample-column", default="Accession",
                        help="Sample/accession column in the workbook")
    parser.add_argument("--allele-column", default="allele",
                        help="Paper-defined Q7B-PH haplotype/allele column")
    parser.add_argument("--phenotype-column", default=PHENOTYPE_COLUMN,
                        help="Plant-height phenotype column")
    parser.add_argument("--unit", choices=["plot", "accession"], default="plot",
                        help="Use each Figure 3g plot row, or aggregate replicates to accession means")
    return parser


def load_figure3g_table(
    source_xlsx: Path,
    sheet_name: str,
    sample_column: str,
    allele_column: str,
    phenotype_column: str,
) -> pd.DataFrame:
    if not source_xlsx.exists():
        raise FileNotFoundError(f"source workbook is missing: {source_xlsx}")

    df = pd.read_excel(source_xlsx, sheet_name=sheet_name)
    required = [sample_column, allele_column, phenotype_column]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"source workbook missing columns: {', '.join(missing)}")

    subset = df[required].copy()
    subset[sample_column] = subset[sample_column].astype(str).str.strip()
    subset[allele_column] = subset[allele_column].astype(str).str.strip()
    subset[phenotype_column] = pd.to_numeric(subset[phenotype_column], errors="coerce")
    subset = subset.replace({"": pd.NA, "nan": pd.NA, "NaN": pd.NA})
    subset = subset.dropna(subset=[sample_column, allele_column, phenotype_column])
    if subset.empty:
        raise ValueError("no complete Q7B-PH rows after filtering source workbook")
    return subset


def write_precomputed_inputs(
    source_df: pd.DataFrame,
    marker_output: Path,
    phenotype_output: Path,
    sample_column: str,
    allele_column: str,
    phenotype_column: str,
    unit: str = "plot",
) -> tuple[Path, Path, pd.DataFrame]:
    if unit == "accession":
        normalized = (
            source_df
            .groupby([sample_column, allele_column], as_index=False)[phenotype_column]
            .mean()
        )
        sample_allele_counts = normalized.groupby(sample_column)[allele_column].nunique()
        ambiguous = sample_allele_counts[sample_allele_counts > 1]
        if not ambiguous.empty:
            raise ValueError(
                "samples assigned to multiple Q7B-PH alleles: "
                + ", ".join(str(x) for x in ambiguous.index.tolist())
            )
        normalized = normalized.rename(columns={
            sample_column: "SampleID",
            allele_column: MARKER_ID,
        })
    else:
        normalized = source_df[[sample_column, allele_column, phenotype_column]].copy()
        normalized.insert(
            0,
            "SampleID",
            [
                f"{sample}__plot{i + 1:04d}"
                for i, sample in enumerate(normalized[sample_column].astype(str).str.strip())
            ],
        )
        normalized = normalized.rename(columns={allele_column: MARKER_ID})

    normalized = normalized[["SampleID", MARKER_ID, phenotype_column]].sort_values("SampleID")

    marker_output.parent.mkdir(parents=True, exist_ok=True)
    phenotype_output.parent.mkdir(parents=True, exist_ok=True)
    normalized[["SampleID", MARKER_ID]].to_csv(marker_output, sep="\t", index=False)
    normalized[["SampleID", phenotype_column]].to_csv(phenotype_output, sep="\t", index=False)
    return marker_output, phenotype_output, normalized


def build_q7b_database(
    marker_output: Path,
    phenotype_output: Path,
    output_root: Path,
    target_id: str,
    phenotype_column: str,
    source_xlsx: Path,
) -> Path:
    db_dir = build_database_from_marker_matrix(
        marker_matrix=marker_output,
        phenotype_table=phenotype_output,
        output_root=output_root,
        target_id=target_id,
        chrom=Q7B_PH_CHROM,
        start=Q7B_PH_START,
        end=Q7B_PH_END,
        phenotype_columns=[phenotype_column],
        marker_columns=[MARKER_ID],
        marker_positions={MARKER_ID: Q7B_PH_START},
        expected_direction="increases_trait",
        sample_column="SampleID",
        min_haplotype_count=1,
    )

    gene_info_path = db_dir / "gene_info.json"
    with open(gene_info_path, "r", encoding="utf-8") as f:
        gene_info = json.load(f)
    gene_info.update({
        "source": "wwwg2b_q7b_ph_figure3g",
        "source_file": str(source_xlsx),
        "source_marker": MARKER_ID,
        "coordinate_note": (
            "Figure 3d Q7B-PH genetic-map peak near 67 cM; this precomputed "
            "database uses the paper-defined Figure 3g allele labels, not a "
            "physical chr7B VCF interval."
        ),
    })
    with open(gene_info_path, "w", encoding="utf-8") as f:
        json.dump(gene_info, f, ensure_ascii=False, indent=2)
    return db_dir


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    source_xlsx = Path(args.source_xlsx)
    marker_output = Path(args.marker_output)
    phenotype_output = Path(args.phenotype_output)
    output_root = Path(args.output_root)

    source_df = load_figure3g_table(
        source_xlsx=source_xlsx,
        sheet_name=args.sheet_name,
        sample_column=args.sample_column,
        allele_column=args.allele_column,
        phenotype_column=args.phenotype_column,
    )
    marker_output, phenotype_output, normalized = write_precomputed_inputs(
        source_df=source_df,
        marker_output=marker_output,
        phenotype_output=phenotype_output,
        sample_column=args.sample_column,
        allele_column=args.allele_column,
        phenotype_column=args.phenotype_column,
        unit=args.unit,
    )
    db_dir = build_q7b_database(
        marker_output=marker_output,
        phenotype_output=phenotype_output,
        output_root=output_root,
        target_id=args.target_id,
        phenotype_column=args.phenotype_column,
        source_xlsx=source_xlsx,
    )

    print(f"[INFO] Source rows: {len(source_df)}")
    print(f"[INFO] Analysis rows ({args.unit} unit): {len(normalized)}")
    print(f"[INFO] Marker matrix: {marker_output}")
    print(f"[INFO] Phenotype table: {phenotype_output}")
    print(f"[INFO] Built star-gene database: {db_dir}")
    print("[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Q7B-HT")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
