#!/usr/bin/env python3
"""Prepare TaGW2-A1 and TaGW2-D1 WatSeq VCF positive-control databases."""

import argparse
import gzip
import io
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

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


GW2_DATA_DIR = Path(r"D:\Desktop\data\GW2")
DEFAULT_PHENOTYPE_XLSX = Path(
    "external_data/wheat_nature_2024/watseq/"
    "Watkins_Collection_WGIN_WISP_DFW_watseq_phenotype_data_JIC.xlsx"
)
DEFAULT_OUTPUT_ROOT = Path("star_gene_database/wheat_nature_2024")
DEFAULT_INTERMEDIATE_ROOT = Path("external_data/wheat_nature_2024/tagw2_ad")

DEFAULT_VCF_A = GW2_DATA_DIR / (
    "chr6A.SNP.Missing-unphasing.ID.ann.finalSID.allele2_retain."
    "hard_retain.InbreedingCoeff_retain.missing_retain.maf_retain.vcf.gz"
)
DEFAULT_VCF_D = GW2_DATA_DIR / (
    "chr6D.SNP.Missing-unphasing.ID.ann.finalSID.allele2_retain."
    "hard_retain.InbreedingCoeff_retain.missing_retain.maf_retain.vcf.gz"
)

PHENOTYPE_SHEETS = {
    "WISP_Watkins_JIC_CFLN10": ("GW_M_g1000grn-CFLN10", "TGW_CFLN10"),
    "DFW_Watkins_JI_CFLN14": ("GW_M_g1000grn-CFLN14", "TGW_CFLN14"),
}
PHENOTYPE_COLUMNS = ["TGW_mean", "TGW_CFLN10", "TGW_CFLN14"]

TARGETS = {
    "TaGW2-A1": {
        "gene_id": "TraesCS6A02G189300",
        "chrom": "6A",
        "gene_start": 237_734_651,
        "gene_end": 237_760_058,
        "strand": "+",
        "vcf_arg": "vcf_a",
    },
    "TaGW2-D1": {
        "gene_id": "TraesCS6D02G176900",
        "chrom": "6D",
        "gene_start": 175_712_228,
        "gene_end": 175_721_507,
        "strand": "+",
        "vcf_arg": "vcf_d",
    },
}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build precomputed star-gene databases for TaGW2-A1 and TaGW2-D1 "
            "from local WatSeq chr6A/chr6D VCFs and Watkins grain-weight phenotypes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--vcf-a", default=str(DEFAULT_VCF_A), help="chr6A WatSeq SNP VCF")
    parser.add_argument("--vcf-d", default=str(DEFAULT_VCF_D), help="chr6D WatSeq SNP VCF")
    parser.add_argument("--phenotype-xlsx", default=str(DEFAULT_PHENOTYPE_XLSX),
                        help="Watkins JIC phenotype workbook")
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT),
                        help="Output star-gene database root")
    parser.add_argument("--intermediate-root", default=str(DEFAULT_INTERMEDIATE_ROOT),
                        help="Directory for marker and phenotype TSV intermediates")
    parser.add_argument("--target", action="append", choices=sorted(TARGETS),
                        help="Target to prepare; repeatable. Defaults to both A1 and D1.")
    parser.add_argument("--promoter-length", type=int, default=2000,
                        help="Upstream promoter length to include for plus-strand genes")
    parser.add_argument("--min-haplotype-count", type=int, default=3,
                        help="Minimum sample count for retained haplotypes")
    parser.add_argument("--max-missing-rate", type=float, default=0.2,
                        help="Maximum per-variant missing rate before complete-sample filtering")
    return parser


def normalize_chrom(value: object) -> str:
    text = str(value).strip()
    return text[3:] if text.lower().startswith("chr") else text


def load_watkins_tgw_phenotypes(phenotype_xlsx: Path) -> pd.DataFrame:
    if not phenotype_xlsx.exists():
        raise FileNotFoundError(f"phenotype workbook is missing: {phenotype_xlsx}")

    merged: Optional[pd.DataFrame] = None
    for sheet_name, (source_col, output_col) in PHENOTYPE_SHEETS.items():
        df = pd.read_excel(phenotype_xlsx, sheet_name=sheet_name)
        required = ["StoreCode", source_col]
        missing = [col for col in required if col not in df.columns]
        if missing:
            raise ValueError(f"{sheet_name} missing columns: {', '.join(missing)}")
        subset = df[required].copy()
        subset = subset.rename(columns={"StoreCode": "SampleID", source_col: output_col})
        subset["SampleID"] = subset["SampleID"].astype(str).str.strip()
        subset[output_col] = pd.to_numeric(subset[output_col], errors="coerce")
        subset = subset.dropna(subset=["SampleID"])
        subset = subset[subset["SampleID"] != ""]
        subset = subset.groupby("SampleID", as_index=False)[output_col].mean()
        merged = subset if merged is None else pd.merge(merged, subset, on="SampleID", how="outer")

    if merged is None or merged.empty:
        raise ValueError("no phenotype rows were loaded")
    merged["TGW_mean"] = merged[["TGW_CFLN10", "TGW_CFLN14"]].mean(axis=1, skipna=True)
    merged = merged.dropna(subset=["TGW_mean"]).sort_values("SampleID").reset_index(drop=True)
    if merged.empty:
        raise ValueError("no complete grain-weight phenotypes after filtering")
    return merged[["SampleID", "TGW_mean", "TGW_CFLN10", "TGW_CFLN14"]]


def genotype_to_state(sample_value: str, fmt: str, ref: str, alt: str) -> Optional[str]:
    fmt_parts = fmt.split(":")
    try:
        gt_index = fmt_parts.index("GT")
    except ValueError:
        gt_index = 0
    sample_parts = str(sample_value).split(":")
    if gt_index >= len(sample_parts):
        return None
    gt = sample_parts[gt_index]
    if not gt or "." in gt:
        return None

    alleles = [ref] + alt.split(",")
    indices = re.split(r"[\/|]", gt)
    states: List[str] = []
    for index in indices:
        if index == ".":
            return None
        try:
            allele_index = int(index)
        except ValueError:
            return None
        if allele_index < 0 or allele_index >= len(alleles):
            return None
        states.append(alleles[allele_index])
    if not states:
        return None
    if len(set(states)) == 1:
        return states[0]
    return "/".join(states)


def marker_id_for_record(chrom: str, pos: int, record_id: str, ref: str, alt: str) -> str:
    if record_id and record_id != ".":
        return record_id
    safe_alt = alt.replace(",", "_")
    return f"chr{normalize_chrom(chrom)}_{pos}_{ref}_{safe_alt}"


def extract_vcf_region_markers(
    vcf_path: Path,
    chrom: str,
    start: int,
    end: int,
    phenotype_samples: Iterable[str],
    max_missing_rate: float,
) -> tuple[pd.DataFrame, Dict[str, int], Dict[str, Dict[str, object]]]:
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF is missing: {vcf_path}")

    target_chrom = normalize_chrom(chrom)
    phenotype_sample_set = {str(sample) for sample in phenotype_samples}
    sample_names: List[str] = []
    sample_indices: List[int] = []
    variants: List[Dict[str, object]] = []

    with gzip.open(vcf_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n\r").split("\t")
                sample_names = header[9:]
                sample_indices = [
                    i for i, sample in enumerate(sample_names, start=9)
                    if sample in phenotype_sample_set
                ]
                if not sample_indices:
                    raise ValueError(f"no VCF samples overlap phenotype samples in {vcf_path}")
                break
        else:
            raise ValueError(f"VCF header not found: {vcf_path}")

        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 10:
                continue
            rec_chrom = normalize_chrom(parts[0])
            if rec_chrom != target_chrom:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if pos < start:
                continue
            if pos > end:
                break

            ref = parts[3]
            alt = parts[4]
            fmt = parts[8]
            marker_id = marker_id_for_record(rec_chrom, pos, parts[2], ref, alt)
            calls: Dict[str, str] = {}
            missing = 0
            for sample_idx in sample_indices:
                sample = sample_names[sample_idx - 9]
                state = genotype_to_state(parts[sample_idx], fmt, ref, alt)
                if state is None:
                    missing += 1
                    continue
                calls[sample] = state

            denominator = len(sample_indices)
            missing_rate = missing / denominator if denominator else 1.0
            if missing_rate > max_missing_rate:
                continue
            if len(set(calls.values())) < 2:
                continue
            variants.append({
                "marker_id": marker_id,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "calls": calls,
                "missing_rate": missing_rate,
            })

    if not variants:
        raise ValueError(f"no polymorphic variants found in {vcf_path} for {chrom}:{start}-{end}")

    complete_samples = sorted(
        set.intersection(*(set(v["calls"].keys()) for v in variants))
    )
    if len(complete_samples) < 2:
        raise ValueError(f"fewer than two complete samples after VCF filtering for {vcf_path}")

    rows = []
    for sample in complete_samples:
        row = {"SampleID": sample}
        for variant in variants:
            row[str(variant["marker_id"])] = variant["calls"][sample]
        rows.append(row)
    marker_df = pd.DataFrame(rows)
    marker_positions = {str(v["marker_id"]): int(v["position"]) for v in variants}
    marker_metadata = {
        str(v["marker_id"]): {
            "position": int(v["position"]),
            "ref": str(v["ref"]),
            "alt": str(v["alt"]),
            "missing_rate": float(v["missing_rate"]),
        }
        for v in variants
    }
    return marker_df, marker_positions, marker_metadata


def region_for_target(target: Dict[str, object], promoter_length: int) -> tuple[int, int]:
    gene_start = int(target["gene_start"])
    gene_end = int(target["gene_end"])
    strand = str(target.get("strand", "+"))
    if strand == "-":
        return gene_start, gene_end + promoter_length
    return max(1, gene_start - promoter_length), gene_end


def update_gene_and_variant_metadata(
    db_dir: Path,
    target_id: str,
    target: Dict[str, object],
    source_vcf: Path,
    phenotype_xlsx: Path,
    marker_metadata: Dict[str, Dict[str, object]],
    promoter_length: int,
) -> None:
    gene_info_path = db_dir / "gene_info.json"
    with open(gene_info_path, "r", encoding="utf-8") as f:
        gene_info = json.load(f)

    region_start, region_end = region_for_target(target, promoter_length)
    gene_start = int(target["gene_start"])
    gene_end = int(target["gene_end"])
    strand = str(target.get("strand", "+"))
    if strand == "-":
        promoter_start = gene_end + 1
        promoter_end = gene_end + promoter_length
    else:
        promoter_start = region_start
        promoter_end = gene_start - 1

    gene_info.update({
        "gene_id": target_id,
        "gene_symbol": target_id,
        "traes_id": target["gene_id"],
        "chrom": target["chrom"],
        "start": region_start,
        "end": region_end,
        "gene_start": gene_start,
        "gene_end": gene_end,
        "strand": strand,
        "promoter_length": promoter_length,
        "promoter_start": promoter_start,
        "promoter_end": promoter_end,
        "source": "watseq_tagw2_ad_vcf",
        "source_vcf": str(source_vcf),
        "source_phenotype": str(phenotype_xlsx),
        "source_note": (
            "WatSeq Chinese Spring RefSeq v1.0 chr6A/chr6D SNP VCF region; "
            "TaGW2-B1 is not included in this A/D-only first pass."
        ),
    })
    with open(gene_info_path, "w", encoding="utf-8") as f:
        json.dump(gene_info, f, ensure_ascii=False, indent=2)

    variant_path = db_dir / "variant_info.csv"
    variant_df = pd.read_csv(variant_path)
    if "marker_id" in variant_df.columns:
        for idx, row in variant_df.iterrows():
            marker_id = str(row["marker_id"])
            meta = marker_metadata.get(marker_id, {})
            pos = int(meta.get("position", row["position"]))
            annotation = "promoter" if promoter_start <= pos <= promoter_end else "intron"
            variant_df.at[idx, "position"] = pos
            variant_df.at[idx, "annotation"] = annotation
            if meta:
                variant_df.at[idx, "ref"] = meta["ref"]
                variant_df.at[idx, "alt"] = meta["alt"]
                variant_df.at[idx, "missing_rate"] = meta["missing_rate"]
    variant_df.to_csv(variant_path, index=False)


def build_target_database(
    target_id: str,
    target: Dict[str, object],
    vcf_path: Path,
    phenotype_df: pd.DataFrame,
    phenotype_xlsx: Path,
    output_root: Path,
    intermediate_root: Path,
    promoter_length: int,
    min_haplotype_count: int,
    max_missing_rate: float,
) -> Path:
    region_start, region_end = region_for_target(target, promoter_length)
    marker_df, marker_positions, marker_metadata = extract_vcf_region_markers(
        vcf_path=vcf_path,
        chrom=str(target["chrom"]),
        start=region_start,
        end=region_end,
        phenotype_samples=phenotype_df["SampleID"].tolist(),
        max_missing_rate=max_missing_rate,
    )
    phenotype_subset = phenotype_df[phenotype_df["SampleID"].isin(marker_df["SampleID"])].copy()
    phenotype_subset = phenotype_subset.sort_values("SampleID")
    marker_df = marker_df[marker_df["SampleID"].isin(phenotype_subset["SampleID"])].sort_values("SampleID")

    if len(marker_df) < 2:
        raise ValueError(f"{target_id}: fewer than two samples after phenotype/VCF merge")

    target_intermediate = intermediate_root / target_id
    marker_output = target_intermediate / "marker_matrix.tsv"
    phenotype_output = target_intermediate / "phenotype.tsv"
    target_intermediate.mkdir(parents=True, exist_ok=True)
    marker_df.to_csv(marker_output, sep="\t", index=False)
    phenotype_subset[["SampleID", *PHENOTYPE_COLUMNS]].to_csv(phenotype_output, sep="\t", index=False)

    db_dir = build_database_from_marker_matrix(
        marker_matrix=marker_output,
        phenotype_table=phenotype_output,
        output_root=output_root,
        target_id=target_id,
        chrom=str(target["chrom"]),
        start=region_start,
        end=region_end,
        phenotype_columns=PHENOTYPE_COLUMNS,
        marker_columns=[c for c in marker_df.columns if c != "SampleID"],
        marker_positions=marker_positions,
        expected_direction="unknown",
        sample_column="SampleID",
        min_haplotype_count=min_haplotype_count,
    )
    update_gene_and_variant_metadata(
        db_dir=db_dir,
        target_id=target_id,
        target=target,
        source_vcf=vcf_path,
        phenotype_xlsx=phenotype_xlsx,
        marker_metadata=marker_metadata,
        promoter_length=promoter_length,
    )
    return db_dir


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    phenotype_xlsx = Path(args.phenotype_xlsx)
    output_root = Path(args.output_root)
    intermediate_root = Path(args.intermediate_root)
    phenotype_df = load_watkins_tgw_phenotypes(phenotype_xlsx)

    selected_targets = args.target or sorted(TARGETS)
    vcf_paths = {
        "vcf_a": Path(args.vcf_a),
        "vcf_d": Path(args.vcf_d),
    }

    print(f"[INFO] Loaded grain-weight phenotypes: {len(phenotype_df)} samples")
    for target_id in selected_targets:
        target = TARGETS[target_id]
        vcf_path = vcf_paths[str(target["vcf_arg"])]
        print(f"[INFO] Preparing {target_id} from {vcf_path}")
        db_dir = build_target_database(
            target_id=target_id,
            target=target,
            vcf_path=vcf_path,
            phenotype_df=phenotype_df,
            phenotype_xlsx=phenotype_xlsx,
            output_root=output_root,
            intermediate_root=intermediate_root,
            promoter_length=args.promoter_length,
            min_haplotype_count=args.min_haplotype_count,
            max_missing_rate=args.max_missing_rate,
        )
        variant_count = len(pd.read_csv(db_dir / "variant_info.csv"))
        sample_count = len(pd.read_csv(db_dir / "phenotype_data.csv"))
        print(f"[INFO] Built {target_id}: {variant_count} variants, {sample_count} samples -> {db_dir}")

    print("[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-A1")
    print("[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-D1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
