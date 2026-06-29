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
    _marker_is_polymorphic,
    _read_aligned_fasta,
    build_database_from_marker_matrix,
    build_marker_matrix_from_aligned_fasta,
)


DEFAULT_GENE_FASTA = Path(
    "external_data/literature/vrn1_full_sequence/esm1_alignments/ESM1/VRNB1_gene.fasta"
)
DEFAULT_PROMOTER_FASTA = Path(
    "external_data/literature/vrn1_full_sequence/esm1_alignments/ESM1/VRNB1_prom.fasta"
)
DEFAULT_ESM2_WORKBOOK = Path(
    "external_data/literature/vrn1_full_sequence/s001_extracted/ESM2.xlsx"
)
DEFAULT_OUTPUT_ROOT = Path("star_gene_database/wheat_nature_2024")
DEFAULT_INTERMEDIATE_ROOT = Path("external_data/wheat_nature_2024/vrn_b1_full_sequence_ijms2021")
DEFAULT_KISS2014_PHENOTYPE_TABLE = Path(
    "external_data/wheat_nature_2024/vrn_kiss2014/VRN-Kiss2014/phenotype.tsv"
)
TARGET_ID = "VRN-B1-fullSequence-IJMS2021"
KISS2014_HEADING_TARGET_ID = "VRN-B1-fullSequence-IJMS2021-Kiss2014Heading"
VRNB1F_837_FORWARD = "ACCATCTCCTTGCTTGCG"
VRNB1F_837_REVERSE = "GACGATACGAACACGACAACC"
VRNB1F_EXPECTED_INSERTION_LENGTH = 837


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert the IJMS 2021 VRN-B1 complete gene-body FASTA alignment "
            "into a base/indel marker database for discovery-style validation."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--gene-fasta", default=str(DEFAULT_GENE_FASTA))
    parser.add_argument("--promoter-fasta", default=str(DEFAULT_PROMOTER_FASTA))
    parser.add_argument("--esm2-workbook", default=str(DEFAULT_ESM2_WORKBOOK))
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument("--intermediate-root", default=str(DEFAULT_INTERMEDIATE_ROOT))
    parser.add_argument("--target-id", default=TARGET_ID)
    parser.add_argument(
        "--phenotype-source",
        choices=["growth_habit", "kiss2014_heading"],
        default="growth_habit",
        help="phenotype table to merge with the IJMS 2021 full-sequence marker matrix",
    )
    parser.add_argument(
        "--continuous-phenotype-table",
        default=str(DEFAULT_KISS2014_PHENOTYPE_TABLE),
        help="continuous phenotype TSV/CSV used when --phenotype-source kiss2014_heading",
    )
    parser.add_argument("--continuous-sample-column", default="SampleID")
    parser.add_argument(
        "--continuous-phenotype-columns",
        default="DEV49_mean,DEV59_mean",
        help="comma-separated continuous phenotype columns",
    )
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


def _normal_sample_key(value: object) -> str:
    text = _clean_text(value).casefold()
    text = text.replace("'", "").replace("’", "")
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return " ".join(text.split())


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


def _read_delimited_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"continuous phenotype table is missing: {path}")
    sep = "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    return pd.read_csv(path, sep=sep)


def merge_continuous_phenotypes(
    ijms_phenotype_df: pd.DataFrame,
    phenotype_table: Path,
    sample_column: str,
    phenotype_columns: Sequence[str],
) -> pd.DataFrame:
    """Merge an external continuous phenotype table to IJMS cultivar names."""
    if "SampleID" not in ijms_phenotype_df.columns:
        raise ValueError("IJMS phenotype table missing SampleID")
    if not phenotype_columns:
        raise ValueError("at least one continuous phenotype column is required")

    source_df = _read_delimited_table(Path(phenotype_table))
    required = {sample_column, *phenotype_columns}
    missing = sorted(required - set(source_df.columns))
    if missing:
        raise ValueError("continuous phenotype table missing columns: " + ", ".join(missing))

    ijms_samples = ijms_phenotype_df[["SampleID"]].dropna().drop_duplicates().copy()
    ijms_samples["_sample_key"] = ijms_samples["SampleID"].map(_normal_sample_key)

    continuous = source_df[[sample_column, *phenotype_columns]].dropna(subset=[sample_column]).copy()
    continuous["_sample_key"] = continuous[sample_column].map(_normal_sample_key)
    continuous = continuous[continuous["_sample_key"] != ""].drop_duplicates(subset=["_sample_key"])
    for col in phenotype_columns:
        continuous[col] = pd.to_numeric(continuous[col], errors="coerce")
    continuous = continuous.dropna(subset=list(phenotype_columns))

    merged = ijms_samples.merge(
        continuous[["_sample_key", *phenotype_columns]],
        on="_sample_key",
        how="inner",
    )
    merged = merged.drop(columns=["_sample_key"])
    if merged.empty:
        raise ValueError("continuous phenotype table has no cultivar overlap with IJMS Table S1")
    return merged.sort_values("SampleID").reset_index(drop=True)


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


def _aligned_fasta_length(fasta_path: Path) -> int:
    seq_parts = []
    in_first_record = False
    with Path(fasta_path).open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            text = line.strip()
            if not text:
                continue
            if text.startswith(">"):
                if in_first_record:
                    break
                in_first_record = True
                continue
            if in_first_record:
                seq_parts.append(text)
    if not seq_parts:
        raise ValueError(f"no FASTA sequence records found: {fasta_path}")
    return len("".join(seq_parts))


def _reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return str(seq).translate(table)[::-1].upper()


def _ungapped_sequence_and_alignment_map(aligned_seq: str) -> tuple[str, list[int]]:
    ungapped = []
    ungapped_to_alignment = []
    for i, base in enumerate(str(aligned_seq).upper()):
        if base != "-":
            ungapped.append(base)
            ungapped_to_alignment.append(i)
    return "".join(ungapped), ungapped_to_alignment


def _find_primer_alignment_hit(aligned_seq: str, primer: str) -> Optional[tuple[int, int, str]]:
    ungapped, ungapped_to_alignment = _ungapped_sequence_and_alignment_map(aligned_seq)
    candidates = [
        (str(primer).upper(), "forward"),
        (_reverse_complement(primer), "reverse_complement"),
    ]
    for query, orientation in candidates:
        idx = ungapped.find(query)
        if idx >= 0:
            end_idx = idx + len(query) - 1
            return ungapped_to_alignment[idx], ungapped_to_alignment[end_idx], orientation
    return None


def _normal_vrn_b1f_sample_label(value: object) -> str:
    return _normal_sample_key(value).replace(" ", "")


def _find_vrn_b1f_reference_and_carriers(
    records: Sequence[tuple[str, str]],
) -> tuple[tuple[str, str], list[tuple[str, str]]]:
    reference = None
    carriers: list[tuple[str, str]] = []
    for sample, seq in records:
        label = _normal_vrn_b1f_sample_label(sample)
        if label == "tdc" or label.endswith("tdc"):
            reference = (sample, seq)
        if (
            label.startswith("anza")
            or label.startswith("barta")
            or label.startswith("marquis01c0201025")
        ):
            carriers.append((sample, seq))
    if reference is None:
        raise ValueError("cannot locate TDC reference sequence in VRN-B1 gene alignment")
    if not carriers:
        raise ValueError("cannot locate Anza/Barta/Marquis carrier sequences in VRN-B1 gene alignment")
    return reference, carriers


def _find_vrn_b1f_837in_alignment_interval(
    records: Sequence[tuple[str, str]],
) -> tuple[int, int, Dict[str, object]]:
    """Locate the Vrn-B1f 837 bp insertion from diagnostic primers and TDC gaps.

    Returns zero-based inclusive alignment coordinates in VRNB1_gene.fasta.
    Literature labels are used only to create an exact validation marker.
    """
    reference, carriers = _find_vrn_b1f_reference_and_carriers(records)
    ref_sample, ref_seq = reference
    forward_hit = _find_primer_alignment_hit(ref_seq, VRNB1F_837_FORWARD)
    reverse_hit = _find_primer_alignment_hit(ref_seq, VRNB1F_837_REVERSE)
    if forward_hit is None or reverse_hit is None:
        raise ValueError("cannot locate VRNB1_837inF/R primer sites in TDC VRN-B1 alignment")

    primer_lo = min(forward_hit[0], reverse_hit[0])
    primer_hi = max(forward_hit[1], reverse_hit[1])
    carrier_seqs = [seq for _, seq in carriers]
    candidate_cols = []
    for i in range(primer_lo, primer_hi + 1):
        ref_base = ref_seq[i]
        if ref_base != "-":
            continue
        if all(seq[i] != "-" for seq in carrier_seqs):
            candidate_cols.append(i)

    runs = []
    if candidate_cols:
        start = prev = candidate_cols[0]
        for i in candidate_cols[1:]:
            if i == prev + 1:
                prev = i
            else:
                runs.append((start, prev))
                start = prev = i
        runs.append((start, prev))

    expected = VRNB1F_EXPECTED_INSERTION_LENGTH
    exact = [run for run in runs if run[1] - run[0] + 1 == expected]
    if exact:
        insertion_start, insertion_end = exact[0]
    elif runs:
        insertion_start, insertion_end = max(runs, key=lambda run: run[1] - run[0] + 1)
    else:
        raise ValueError("no TDC-gap/carrier-base interval found for Vrn-B1f")

    metadata = {
        "reference_sample": ref_sample,
        "carrier_samples": ";".join(sample for sample, _ in carriers),
        "primer_forward": "VRNB1_837inF",
        "primer_reverse": "VRNB1_837inR",
        "primer_forward_sequence": VRNB1F_837_FORWARD,
        "primer_reverse_sequence": VRNB1F_837_REVERSE,
        "diagnostic_amplicon_reference_bases": sum(
            1 for base in ref_seq[primer_lo:primer_hi + 1] if base != "-"
        ),
        "diagnostic_amplicon_carrier_bases": sum(
            1 for base in carrier_seqs[0][primer_lo:primer_hi + 1] if base != "-"
        ),
    }
    return insertion_start, insertion_end, metadata


def _vrn_b1f_validation_allele(block: str) -> str:
    if not block:
        return "N"
    allowed = set("ACGTRYSWKMBDHVN-")
    text = str(block).upper()
    if any(base not in allowed for base in text):
        return "N"
    ungapped = text.replace("-", "")
    if not ungapped:
        return f"DEL_{len(text)}"
    return ungapped


def _add_vrn_b1f_837in_validation_marker(
    marker_df: pd.DataFrame,
    marker_metadata: list[Dict[str, object]],
    gene_fasta: Path,
    sample_id_fn: Callable[[str], str],
    promoter_length: int,
    min_minor_count: int,
) -> tuple[pd.DataFrame, list[Dict[str, object]]]:
    records = _read_aligned_fasta(gene_fasta, sample_id_fn=sample_id_fn)
    try:
        insertion_start, insertion_end, diagnostic_meta = _find_vrn_b1f_837in_alignment_interval(records)
    except ValueError:
        return marker_df, marker_metadata
    length = insertion_end - insertion_start + 1
    combined_start = int(promoter_length) + insertion_start + 1
    combined_end = int(promoter_length) + insertion_end + 1
    marker_id = f"VRNB1gene_insertion_{VRNB1F_EXPECTED_INSERTION_LENGTH}_VrnB1f_{combined_start}_{combined_end}"

    allele_by_sample = {
        sample: _vrn_b1f_validation_allele(seq[insertion_start:insertion_end + 1])
        for sample, seq in records
    }
    alleles = [allele_by_sample.get(sample, "N") for sample in marker_df["SampleID"]]
    if not _marker_is_polymorphic(alleles, min_minor_count=min_minor_count):
        return marker_df, marker_metadata

    marker_df = marker_df.copy()
    marker_df[marker_id] = alleles
    marker_metadata = [*marker_metadata, {
        "marker_id": marker_id,
        "kind": "indel",
        "alignment_start": combined_start,
        "alignment_end": combined_end,
        "length": length,
        "missing_rate": round(sum(allele == "N" for allele in alleles) / max(len(alleles), 1), 6),
        "segment": "gene",
        "source_marker_id": marker_id,
        "annotation": "diagnostic_marker",
        "validation_marker": True,
        "literature_variant": "Vrn-B1f_837bp_insertion",
        "source": "ijms2021_table_s3_diagnostic_primers_and_tdc_vs_carrier_alignment",
        **diagnostic_meta,
    }]
    return marker_df, marker_metadata


def _shift_marker_metadata(
    marker_metadata: Sequence[Dict[str, object]],
    prefix: str,
    offset: int,
    segment: str,
) -> list[Dict[str, object]]:
    shifted = []
    for row in marker_metadata:
        source_marker_id = str(row["marker_id"])
        new_start = int(row["alignment_start"]) + offset
        new_end = int(row["alignment_end"]) + offset
        kind = str(row.get("kind", ""))
        if kind == "snp":
            marker_id = f"{prefix}_snp_{new_start}"
        elif kind == "indel":
            marker_id = f"{prefix}_indel_{new_start}_{new_end}"
        else:
            suffix = source_marker_id.split("_", 1)[1] if "_" in source_marker_id else source_marker_id
            marker_id = f"{prefix}_{suffix}"
        shifted.append({
            **dict(row),
            "marker_id": marker_id,
            "alignment_start": new_start,
            "alignment_end": new_end,
            "segment": segment,
            "source_marker_id": source_marker_id,
        })
    return shifted


def _rename_marker_columns(marker_df: pd.DataFrame, marker_metadata: Sequence[Dict[str, object]]) -> pd.DataFrame:
    rename_map = {
        str(row["source_marker_id"]): str(row["marker_id"])
        for row in marker_metadata
        if "source_marker_id" in row
    }
    return marker_df.rename(columns=rename_map)


def build_full_sequence_marker_matrix(
    promoter_fasta: Path,
    gene_fasta: Path,
    sample_id_fn: Callable[[str], str],
    max_missing_rate: float,
    min_minor_count: int,
) -> tuple[pd.DataFrame, list[Dict[str, object]], Dict[str, int]]:
    """Build one VRN-B1 marker matrix from promoter and gene-body alignments."""
    promoter_fasta = Path(promoter_fasta)
    gene_fasta = Path(gene_fasta)
    promoter_length = _aligned_fasta_length(promoter_fasta)
    gene_length = _aligned_fasta_length(gene_fasta)

    promoter_df, promoter_metadata = build_marker_matrix_from_aligned_fasta(
        promoter_fasta,
        marker_prefix="VRNB1prom",
        sample_id_fn=sample_id_fn,
        max_missing_rate=max_missing_rate,
        min_minor_count=min_minor_count,
    )
    gene_df, gene_metadata = build_marker_matrix_from_aligned_fasta(
        gene_fasta,
        marker_prefix="VRNB1gene",
        sample_id_fn=sample_id_fn,
        max_missing_rate=max_missing_rate,
        min_minor_count=min_minor_count,
    )

    promoter_metadata = _shift_marker_metadata(
        promoter_metadata,
        prefix="VRNB1prom",
        offset=0,
        segment="promoter",
    )
    gene_metadata = _shift_marker_metadata(
        gene_metadata,
        prefix="VRNB1gene",
        offset=promoter_length,
        segment="gene",
    )
    promoter_df = _rename_marker_columns(promoter_df, promoter_metadata)
    gene_df = _rename_marker_columns(gene_df, gene_metadata)
    marker_df = promoter_df.merge(gene_df, on="SampleID", how="inner")
    if marker_df.empty:
        raise ValueError("no samples overlap between VRN-B1 promoter and gene alignments")
    marker_df, marker_metadata = _add_vrn_b1f_837in_validation_marker(
        marker_df,
        [*promoter_metadata, *gene_metadata],
        gene_fasta=gene_fasta,
        sample_id_fn=sample_id_fn,
        promoter_length=promoter_length,
        min_minor_count=min_minor_count,
    )

    layout = {
        "promoter_start": 1,
        "promoter_end": promoter_length,
        "gene_start": promoter_length + 1,
        "gene_end": promoter_length + gene_length,
        "promoter_length": promoter_length,
        "gene_length": gene_length,
    }
    return marker_df, marker_metadata, layout


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
    promoter_fasta: Path,
    gene_fasta: Path,
    esm2_workbook: Path,
    phenotype_source: str,
    phenotype_source_note: str,
    marker_metadata: Sequence[Dict[str, object]],
    marker_df: pd.DataFrame,
    phenotype_df: pd.DataFrame,
    layout: Dict[str, int],
) -> None:
    metadata_by_id = {str(row["marker_id"]): dict(row) for row in marker_metadata}

    gene_info_path = db_dir / "gene_info.json"
    with gene_info_path.open("r", encoding="utf-8") as f:
        gene_info = json.load(f)
    source_id = (
        "ijms2021_vrn_b1_full_gene_alignment_plus_kiss2014_heading"
        if phenotype_source == "kiss2014_heading"
        else "ijms2021_vrn_b1_full_gene_alignment"
    )
    gene_info.update({
        "gene_id": target_id,
        "gene_symbol": "VRN-B1",
        "chrom": "VRN-B1_IJMS2021_alignment",
        "start": 1,
        "end": int(layout["gene_end"]),
        "gene_start": int(layout["gene_start"]),
        "gene_end": int(layout["gene_end"]),
        "length": int(layout["gene_end"]),
        "promoter_start": int(layout["promoter_start"]),
        "promoter_end": int(layout["promoter_end"]),
        "promoter_length": int(layout["promoter_length"]),
        "promoter_actual_length": int(layout["promoter_length"]),
        "source": source_id,
        "source_promoter_fasta": str(promoter_fasta),
        "source_fasta": str(gene_fasta),
        "source_workbook": str(esm2_workbook),
        "phenotype_source": phenotype_source,
        "source_note": (
            "In-Depth Sequence Analysis of Bread Wheat VRN1 Genes, IJMS 2021 "
            "ESM1 VRNB1_prom.fasta plus VRNB1_gene.fasta. "
            f"{phenotype_source_note} "
            "Discovery markers are all polymorphic A/C/G/T columns and indel blocks "
            "retained from the full promoter plus gene-body alignment."
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
        segment = str(meta.get("segment", ""))
        if meta.get("annotation"):
            annotation = str(meta.get("annotation"))
        elif segment == "promoter":
            annotation = "promoter"
        else:
            annotation = "base_indel" if kind == "indel" else "base_snp"
        variant_df.at[idx, "annotation"] = annotation
        variant_df.at[idx, "alignment_start"] = meta.get("alignment_start", "")
        variant_df.at[idx, "alignment_end"] = meta.get("alignment_end", "")
        variant_df.at[idx, "alignment_length"] = meta.get("length", "")
        variant_df.at[idx, "alignment_segment"] = segment
        if meta.get("source"):
            source = str(meta.get("source"))
        elif segment == "promoter":
            source = "ijms2021_vrnb1_promoter_alignment"
        else:
            source = "ijms2021_vrnb1_full_gene_alignment"
        variant_df.at[idx, "source"] = source
        for key in (
            "validation_marker",
            "literature_variant",
            "reference_sample",
            "carrier_samples",
            "primer_forward",
            "primer_reverse",
            "primer_forward_sequence",
            "primer_reverse_sequence",
            "diagnostic_amplicon_reference_bases",
            "diagnostic_amplicon_carrier_bases",
        ):
            if key in meta:
                variant_df.at[idx, key] = meta.get(key)
    variant_df.to_csv(variant_path, index=False)


def build_database(
    promoter_fasta: Path,
    gene_fasta: Path,
    esm2_workbook: Path,
    output_root: Path,
    intermediate_root: Path,
    target_id: str,
    phenotype_source: str,
    continuous_phenotype_table: Path,
    continuous_sample_column: str,
    continuous_phenotype_columns: Sequence[str],
    max_missing_rate: float,
    min_minor_count: int,
    min_haplotype_count: int,
) -> Path:
    code_to_sample, growth_habit_df = build_sample_maps_from_esm2(esm2_workbook)
    if phenotype_source == "growth_habit":
        phenotype_df = growth_habit_df
        phenotype_columns = ["GrowthHabitSpringScore"]
        expected_direction = "increases_trait"
        phenotype_source_note = "ESM2 Table S1 supplies spring/winter growth habit."
    elif phenotype_source == "kiss2014_heading":
        phenotype_columns = list(continuous_phenotype_columns)
        phenotype_df = merge_continuous_phenotypes(
            growth_habit_df,
            phenotype_table=continuous_phenotype_table,
            sample_column=continuous_sample_column,
            phenotype_columns=phenotype_columns,
        )
        expected_direction = "decreases_trait"
        phenotype_source_note = (
            "Continuous heading-date phenotypes are merged by cultivar name from "
            f"{continuous_phenotype_table}; IJMS Table S1 is used only for cultivar-code mapping."
        )
    else:
        raise ValueError(f"unsupported phenotype source: {phenotype_source}")

    sample_id_fn = make_vrn_b1_fasta_sample_id_fn(code_to_sample)
    marker_df, marker_metadata, layout = build_full_sequence_marker_matrix(
        promoter_fasta,
        gene_fasta,
        sample_id_fn=sample_id_fn,
        max_missing_rate=max_missing_rate,
        min_minor_count=min_minor_count,
    )
    if not marker_metadata:
        raise ValueError("No polymorphic base/indel markers were retained from the VRN-B1 alignment")

    source_summary = {
        "promoter_fasta": str(promoter_fasta),
        "gene_fasta": str(gene_fasta),
        "esm2_workbook": str(esm2_workbook),
        "phenotype_source": phenotype_source,
        "continuous_phenotype_table": str(continuous_phenotype_table) if phenotype_source != "growth_habit" else "",
        "phenotype_columns": phenotype_columns,
        "alignment_records": int(len(marker_df)),
        "phenotype_samples": int(len(phenotype_df)),
        "retained_markers": int(len(marker_metadata)),
        "promoter_length": int(layout["promoter_length"]),
        "gene_length": int(layout["gene_length"]),
        "promoter_marker_count": int(sum(1 for row in marker_metadata if row.get("segment") == "promoter")),
        "gene_marker_count": int(sum(1 for row in marker_metadata if row.get("segment") == "gene")),
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
        end=int(layout["gene_end"]),
        phenotype_columns=phenotype_columns,
        marker_columns=marker_columns,
        marker_positions=marker_positions,
        expected_direction=expected_direction,
        min_haplotype_count=min_haplotype_count,
    )
    update_database_metadata(
        db_dir,
        target_id,
        promoter_fasta,
        gene_fasta,
        esm2_workbook,
        phenotype_source,
        phenotype_source_note,
        marker_metadata,
        marker_df,
        phenotype_df,
        layout,
    )
    return db_dir


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    phenotype_columns = [
        col.strip()
        for col in str(args.continuous_phenotype_columns).split(",")
        if col.strip()
    ]
    target_id = args.target_id
    if args.phenotype_source == "kiss2014_heading" and target_id == TARGET_ID:
        target_id = KISS2014_HEADING_TARGET_ID
    db_dir = build_database(
        promoter_fasta=Path(args.promoter_fasta),
        gene_fasta=Path(args.gene_fasta),
        esm2_workbook=Path(args.esm2_workbook),
        output_root=Path(args.output_root),
        intermediate_root=Path(args.intermediate_root),
        target_id=target_id,
        phenotype_source=args.phenotype_source,
        continuous_phenotype_table=Path(args.continuous_phenotype_table),
        continuous_sample_column=args.continuous_sample_column,
        continuous_phenotype_columns=phenotype_columns,
        max_missing_rate=args.max_missing_rate,
        min_minor_count=args.min_minor_count,
        min_haplotype_count=args.min_haplotype_count,
    )
    print(f"[INFO] Built IJMS 2021 VRN-B1 full-sequence database: {db_dir}")
    print("[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 "
          f"--target {target_id} --score-mode robust_discovery")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
