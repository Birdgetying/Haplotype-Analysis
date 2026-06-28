#!/usr/bin/env python3
"""Prepare a base-level VRN-B1 full-sequence validation database."""

import argparse
import io
import json
import re
import sys
from pathlib import Path
from typing import Callable, Dict, Optional, Sequence

if sys.platform == "win32":
    try:
        if not isinstance(sys.stdout, io.TextIOWrapper) or sys.stdout.encoding != "utf-8":
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
        if not isinstance(sys.stderr, io.TextIOWrapper) or sys.stderr.encoding != "utf-8":
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
    except (ValueError, AttributeError):
        pass

import pandas as pd

from star_gene_data import (
    build_database_from_marker_matrix,
    build_marker_matrix_from_aligned_fasta,
)


DEFAULT_GENE_FASTA = Path(
    "external_data/literature/vrn1_full_sequence/esm1_alignments/ESM1/VRNB1_gene.fasta"
)
DEFAULT_ESM2_WORKBOOK = Path(
    "external_data/literature/vrn1_full_sequence/s001_extracted/ESM2.xlsx"
)
DEFAULT_OUTPUT_ROOT = Path("star_gene_database/wheat_nature_2024")
DEFAULT_INTERMEDIATE_ROOT = Path("external_data/wheat_nature_2024/vrn_b1_full_sequence_ijms2021")
TARGET_ID = "VRN-B1-fullSequence-IJMS2021"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert the IJMS 2021 VRN-B1 complete gene-body FASTA alignment "
            "into a base/indel marker database for discovery-style validation."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--gene-fasta", default=str(DEFAULT_GENE_FASTA))
    parser.add_argument("--esm2-workbook", default=str(DEFAULT_ESM2_WORKBOOK))
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument("--intermediate-root", default=str(DEFAULT_INTERMEDIATE_ROOT))
    parser.add_argument("--target-id", default=TARGET_ID)
    parser.add_argument("--max-missing-rate", type=float, default=0.0)
    parser.add_argument("--min-minor-count", type=int, default=1)
    parser.add_argument("--min-haplotype-count", type=int, default=1)
    return parser


def _clean_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _normal_code(value: object) -> str:
    text = _clean_text(value)
    if not text:
        return ""
    try:
        numeric = float(text)
    except ValueError:
        return text.upper().replace(" ", "")
    if numeric.is_integer():
        return str(int(numeric))
    return text.upper().replace(" ", "")


def build_sample_maps_from_esm2(workbook: Path) -> tuple[Dict[str, str], pd.DataFrame]:
    """Return FASTA-code-to-cultivar mapping and a spring/winter phenotype table."""
    workbook = Path(workbook)
    if not workbook.exists():
        raise FileNotFoundError(f"IJMS 2021 ESM2 workbook is missing: {workbook}")

    s1 = pd.read_excel(workbook, sheet_name="Table S1", header=4)
    required = {"Code", "Name", "Habit"}
    missing = sorted(required - set(s1.columns))
    if missing:
        raise ValueError("Table S1 missing columns: " + ", ".join(missing))

    code_to_sample: Dict[str, str] = {}
    phenotype_rows = []
    for _, row in s1.dropna(subset=["Code", "Name", "Habit"]).iterrows():
        code = _normal_code(row["Code"])
        sample = _clean_text(row["Name"])
        habit = _clean_text(row["Habit"]).upper()
        if not code or not sample:
            continue
        code_to_sample[code] = sample
        if habit.startswith("S"):
            spring_score = 1.0
        elif habit.startswith("W"):
            spring_score = 0.0
        else:
            continue
        phenotype_rows.append({
            "SampleID": sample,
            "GrowthHabitSpringScore": spring_score,
        })

    if not code_to_sample:
        raise ValueError("No cultivar codes found in Table S1")
    phenotype_df = pd.DataFrame(phenotype_rows).drop_duplicates(subset=["SampleID"])
    if phenotype_df.empty:
        raise ValueError("No spring/winter phenotype rows found in Table S1")
    return code_to_sample, phenotype_df


def make_vrn_b1_fasta_sample_id_fn(code_to_sample: Dict[str, str]) -> Callable[[str], str]:
    def sample_id(header: str) -> str:
        text = header.strip()
        match = re.match(r"^(S\d+)", text, flags=re.IGNORECASE)
        if match:
            code = match.group(1).upper()
            return code_to_sample.get(code, code)
        match = re.match(r"^(\d+)(?:_|$)", text)
        if match:
            code = match.group(1)
            return code_to_sample.get(code, code)
        token = text.split()[0]
        parts = token.split("_")
        cut = len(parts)
        for marker in ("A", "B", "D", "NODE"):
            if marker in parts:
                cut = min(cut, parts.index(marker))
        name_parts = parts[:cut]
        if len(name_parts) > 1 and name_parts[-1].isdigit():
            name_parts = name_parts[:-1]
        return " ".join(name_parts) if name_parts else token

    return sample_id


def write_intermediate_tables(
    marker_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    marker_metadata: Sequence[Dict[str, object]],
    intermediate_root: Path,
    target_id: str,
    source_summary: Dict[str, object],
) -> tuple[Path, Path, Path]:
    target_dir = intermediate_root / target_id
    target_dir.mkdir(parents=True, exist_ok=True)
    marker_path = target_dir / "marker_matrix.tsv"
    phenotype_path = target_dir / "phenotype.tsv"
    metadata_path = target_dir / "marker_metadata.tsv"
    summary_path = target_dir / "source_summary.json"

    marker_df.to_csv(marker_path, sep="\t", index=False)
    phenotype_df.to_csv(phenotype_path, sep="\t", index=False)
    pd.DataFrame(marker_metadata).to_csv(metadata_path, sep="\t", index=False)
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(source_summary, f, ensure_ascii=False, indent=2)
    return marker_path, phenotype_path, metadata_path


def update_database_metadata(
    db_dir: Path,
    target_id: str,
    gene_fasta: Path,
    esm2_workbook: Path,
    marker_metadata: Sequence[Dict[str, object]],
    marker_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
) -> None:
    metadata_by_id = {str(row["marker_id"]): dict(row) for row in marker_metadata}

    gene_info_path = db_dir / "gene_info.json"
    with gene_info_path.open("r", encoding="utf-8") as f:
        gene_info = json.load(f)
    gene_info.update({
        "gene_id": target_id,
        "gene_symbol": "VRN-B1",
        "chrom": "VRN-B1_IJMS2021_alignment",
        "start": 1,
        "end": int(max(row["alignment_end"] for row in marker_metadata)) if marker_metadata else 1,
        "gene_start": 1,
        "gene_end": int(max(row["alignment_end"] for row in marker_metadata)) if marker_metadata else 1,
        "source": "ijms2021_vrn_b1_full_gene_alignment",
        "source_fasta": str(gene_fasta),
        "source_workbook": str(esm2_workbook),
        "source_note": (
            "In-Depth Sequence Analysis of Bread Wheat VRN1 Genes, IJMS 2021 "
            "ESM1 VRNB1_gene.fasta plus ESM2 Table S1 growth habit. "
            "Discovery markers are all polymorphic A/C/G/T columns and indel blocks "
            "retained from the full gene-body alignment."
        ),
        "n_alignment_records": int(len(marker_df)),
        "n_phenotype_samples": int(len(phenotype_df)),
        "n_base_markers": int(len(marker_metadata)),
        "validation_policy": (
            "Literature variants are post hoc audit only; no literature allele or "
            "haplotype label is used to build discovery markers."
        ),
    })
    with gene_info_path.open("w", encoding="utf-8") as f:
        json.dump(gene_info, f, ensure_ascii=False, indent=2)

    variant_path = db_dir / "variant_info.csv"
    variant_df = pd.read_csv(variant_path)
    for idx, row in variant_df.iterrows():
        marker_id = str(row["marker_id"])
        meta = metadata_by_id.get(marker_id, {})
        kind = meta.get("kind", "")
        variant_df.at[idx, "annotation"] = "base_indel" if kind == "indel" else "base_snp"
        variant_df.at[idx, "alignment_start"] = meta.get("alignment_start", "")
        variant_df.at[idx, "alignment_end"] = meta.get("alignment_end", "")
        variant_df.at[idx, "alignment_length"] = meta.get("length", "")
        variant_df.at[idx, "source"] = "ijms2021_vrnb1_full_gene_alignment"
    variant_df.to_csv(variant_path, index=False)


def build_database(
    gene_fasta: Path,
    esm2_workbook: Path,
    output_root: Path,
    intermediate_root: Path,
    target_id: str,
    max_missing_rate: float,
    min_minor_count: int,
    min_haplotype_count: int,
) -> Path:
    code_to_sample, phenotype_df = build_sample_maps_from_esm2(esm2_workbook)
    sample_id_fn = make_vrn_b1_fasta_sample_id_fn(code_to_sample)
    marker_df, marker_metadata = build_marker_matrix_from_aligned_fasta(
        gene_fasta,
        marker_prefix="VRNB1gene",
        sample_id_fn=sample_id_fn,
        max_missing_rate=max_missing_rate,
        min_minor_count=min_minor_count,
    )
    if not marker_metadata:
        raise ValueError("No polymorphic base/indel markers were retained from the VRN-B1 alignment")

    source_summary = {
        "gene_fasta": str(gene_fasta),
        "esm2_workbook": str(esm2_workbook),
        "alignment_records": int(len(marker_df)),
        "phenotype_samples": int(len(phenotype_df)),
        "retained_markers": int(len(marker_metadata)),
        "max_missing_rate": max_missing_rate,
        "min_minor_count": min_minor_count,
        "unmatched_alignment_records": sorted(
            set(marker_df["SampleID"]) - set(phenotype_df["SampleID"])
        ),
    }
    marker_path, phenotype_path, _ = write_intermediate_tables(
        marker_df,
        phenotype_df,
        marker_metadata,
        intermediate_root,
        target_id,
        source_summary,
    )
    marker_columns = [str(row["marker_id"]) for row in marker_metadata]
    marker_positions = {
        str(row["marker_id"]): int(row["alignment_start"])
        for row in marker_metadata
    }
    db_dir = build_database_from_marker_matrix(
        marker_matrix=marker_path,
        phenotype_table=phenotype_path,
        output_root=output_root,
        target_id=target_id,
        chrom="VRN-B1_IJMS2021_alignment",
        start=1,
        end=max(marker_positions.values()),
        phenotype_columns=["GrowthHabitSpringScore"],
        marker_columns=marker_columns,
        marker_positions=marker_positions,
        expected_direction="increases_trait",
        min_haplotype_count=min_haplotype_count,
    )
    update_database_metadata(db_dir, target_id, gene_fasta, esm2_workbook, marker_metadata, marker_df, phenotype_df)
    return db_dir


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    db_dir = build_database(
        gene_fasta=Path(args.gene_fasta),
        esm2_workbook=Path(args.esm2_workbook),
        output_root=Path(args.output_root),
        intermediate_root=Path(args.intermediate_root),
        target_id=args.target_id,
        max_missing_rate=args.max_missing_rate,
        min_minor_count=args.min_minor_count,
        min_haplotype_count=args.min_haplotype_count,
    )
    print(f"[INFO] Built IJMS 2021 VRN-B1 full-sequence database: {db_dir}")
    print("[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 "
          f"--target {args.target_id} --score-mode robust_discovery")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
