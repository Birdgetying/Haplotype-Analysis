#!/usr/bin/env python3
"""Prepare VRN-A1/VRN-B1/VRN-D1 SNP-VCF star-gene databases."""

import argparse
import gzip
import io
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

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


DEFAULT_PHENOTYPE_XLSX = Path(
    "external_data/wheat_nature_2024/watseq/"
    "Watkins_Collection_WGIN_WISP_DFW_watseq_phenotype_data_JIC.xlsx"
)
DEFAULT_OUTPUT_ROOT = Path("star_gene_database/wheat_nature_2024")
DEFAULT_INTERMEDIATE_ROOT = Path("external_data/wheat_nature_2024/vrn_remote")

PHENOTYPE_COLUMNS = ["GrowthHabitSpringScore_CFLN06", "PlantHeight_CFLN06"]

TARGETS = {
    "VRN-A1-remoteSNP": {
        "gene_or_locus": "VRN-A1",
        "gene_id": "TraesCS5A01G391700",
        "chrom": "5A",
        "region_chrom": "chr5A",
        "gene_start": 587_411_454,
        "gene_end": 587_423_416,
        "region_start": 587_409_454,
        "region_end": 587_425_416,
        "strand": "+",
        "vcf": "VRN-A1.wheatomics_snp.vcf.gz",
    },
    "VRN-B1-remoteSNP": {
        "gene_or_locus": "VRN-B1",
        "gene_id": "TraesCS5B01G396600",
        "chrom": "5B",
        "region_chrom": "chr5B",
        "gene_start": 573_802_883,
        "gene_end": 573_816_070,
        "region_start": 573_800_883,
        "region_end": 573_818_070,
        "strand": "+",
        "vcf": "VRN-B1.wheatomics_snp.vcf.gz",
    },
    "VRN-D1-remoteSNP": {
        "gene_or_locus": "VRN-D1",
        "gene_id": "TraesCS5D01G401500",
        "chrom": "5D",
        "region_chrom": "chr5D",
        "gene_start": 467_176_609,
        "gene_end": 467_184_508,
        "region_start": 467_174_609,
        "region_end": 467_186_508,
        "strand": "+",
        "vcf": "VRN-D1.wheatomics_snp.vcf.gz",
    },
}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build separate VRN-A1/VRN-B1/VRN-D1 star-gene databases from "
            "WheatOmics SNP-only micro-VCFs and Watkins growth-habit phenotypes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--phenotype-xlsx", default=str(DEFAULT_PHENOTYPE_XLSX),
                        help="Watkins JIC phenotype workbook")
    parser.add_argument("--vcf-dir", default=str(DEFAULT_INTERMEDIATE_ROOT),
                        help="Directory containing VRN *.wheatomics_snp.vcf.gz files")
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT),
                        help="Output star-gene database root")
    parser.add_argument("--intermediate-root", default=str(DEFAULT_INTERMEDIATE_ROOT),
                        help="Directory for marker and phenotype TSV intermediates")
    parser.add_argument("--target", action="append", choices=sorted(TARGETS),
                        help="Target to prepare; repeatable. Defaults to all VRN targets.")
    parser.add_argument("--min-haplotype-count", type=int, default=3,
                        help="Minimum sample count for retained haplotypes")
    parser.add_argument("--max-missing-rate", type=float, default=0.2,
                        help="Maximum per-variant missing rate before complete-sample filtering")
    return parser


def normalize_chrom(value: object) -> str:
    text = str(value).strip()
    return text[3:] if text.lower().startswith("chr") else text


def spring_score(value: object) -> Optional[float]:
    text = str(value).strip().lower()
    if not text or text in {"nan", "na", ".", "*"}:
        return None
    if not re.fullmatch(r"[sw]+", text):
        return None
    return text.count("s") / len(text)


def load_watkins_growth_habit_phenotypes(phenotype_xlsx: Path) -> pd.DataFrame:
    if not phenotype_xlsx.exists():
        raise FileNotFoundError(f"phenotype workbook is missing: {phenotype_xlsx}")

    df = pd.read_excel(phenotype_xlsx, sheet_name="WGIN_Watkins_JIC_CFLN06")
    required = ["StoreCode", "GrwHabit_E_sw-CFLN06", "PH_M_cm-CFLN06"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"WGIN_Watkins_JIC_CFLN06 missing columns: {', '.join(missing)}")

    out = df[required].copy()
    out = out.rename(columns={
        "StoreCode": "SampleID",
        "PH_M_cm-CFLN06": "PlantHeight_CFLN06",
    })
    out["SampleID"] = out["SampleID"].astype(str).str.strip()
    out["GrowthHabitSpringScore_CFLN06"] = out["GrwHabit_E_sw-CFLN06"].map(spring_score)
    out["PlantHeight_CFLN06"] = pd.to_numeric(out["PlantHeight_CFLN06"], errors="coerce")
    out = out.dropna(subset=["SampleID", "GrowthHabitSpringScore_CFLN06"])
    out = out[out["SampleID"] != ""]
    out = out.groupby("SampleID", as_index=False).agg({
        "GrowthHabitSpringScore_CFLN06": "mean",
        "PlantHeight_CFLN06": "mean",
    })
    if out.empty:
        raise ValueError("no complete Watkins growth-habit phenotypes after filtering")
    return out[["SampleID", *PHENOTYPE_COLUMNS]]


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
    states: List[str] = []
    for index in re.split(r"[\/|]", gt):
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
    target: Dict[str, object],
    phenotype_samples: Iterable[str],
    max_missing_rate: float,
) -> tuple[pd.DataFrame, Dict[str, int], Dict[str, Dict[str, object]]]:
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF is missing: {vcf_path}")

    target_chrom = normalize_chrom(target["chrom"])
    region_start = int(target["region_start"])
    region_end = int(target["region_end"])
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
            if pos < region_start:
                continue
            if pos > region_end:
                break

            ref = parts[3]
            alt = parts[4]
            fmt = parts[8]
            marker_id = marker_id_for_record(rec_chrom, pos, parts[2], ref, alt)
            calls: Dict[str, str] = {}
            missing = 0
            allele_counts: Counter[str] = Counter()
            for sample_idx in sample_indices:
                sample = sample_names[sample_idx - 9]
                state = genotype_to_state(parts[sample_idx], fmt, ref, alt)
                if state is None:
                    missing += 1
                    continue
                calls[sample] = state
                allele_counts[state] += 1

            denominator = len(sample_indices)
            missing_rate = missing / denominator if denominator else 1.0
            if missing_rate > max_missing_rate:
                continue
            if len(allele_counts) < 2:
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
        raise ValueError(
            f"no polymorphic variants found in {vcf_path} for "
            f"{target['region_chrom']}:{region_start}-{region_end}"
        )

    complete_samples = sorted(set.intersection(*(set(v["calls"].keys()) for v in variants)))
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


def annotation_for_position(target: Dict[str, object], pos: int) -> str:
    if int(target["gene_start"]) <= pos <= int(target["gene_end"]):
        return "gene_region"
    return "promoter_or_flank"


def update_gene_and_variant_metadata(
    db_dir: Path,
    target_id: str,
    target: Dict[str, object],
    source_vcf: Path,
    phenotype_xlsx: Path,
    marker_metadata: Dict[str, Dict[str, object]],
) -> None:
    gene_info_path = db_dir / "gene_info.json"
    with gene_info_path.open("r", encoding="utf-8") as f:
        gene_info = json.load(f)

    gene_info.update({
        "gene_id": target_id,
        "gene_symbol": target["gene_or_locus"],
        "gene_or_locus": target["gene_or_locus"],
        "traes_id": target["gene_id"],
        "chrom": target["chrom"],
        "start": int(target["region_start"]),
        "end": int(target["region_end"]),
        "gene_start": int(target["gene_start"]),
        "gene_end": int(target["gene_end"]),
        "strand": target["strand"],
        "expected_direction": "increases_trait",
        "source": "wheatomics_remote_vrn_snp_vcf",
        "source_vcf": str(source_vcf),
        "source_phenotype": str(phenotype_xlsx),
        "source_note": (
            "Single-gene VRN SNP-only micro-VCF from WheatOmics merged SNP VCF; "
            "known VRN promoter/intron-1 deletion and CNV/SV alleles require a separate INDEL/SV source."
        ),
    })
    with gene_info_path.open("w", encoding="utf-8") as f:
        json.dump(gene_info, f, ensure_ascii=False, indent=2)

    variant_path = db_dir / "variant_info.csv"
    variant_df = pd.read_csv(variant_path)
    if "marker_id" in variant_df.columns:
        for idx, row in variant_df.iterrows():
            marker_id = str(row["marker_id"])
            meta = marker_metadata.get(marker_id, {})
            pos = int(meta.get("position", row["position"]))
            variant_df.at[idx, "position"] = pos
            variant_df.at[idx, "annotation"] = annotation_for_position(target, pos)
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
    min_haplotype_count: int,
    max_missing_rate: float,
) -> Path:
    marker_df, marker_positions, marker_metadata = extract_vcf_region_markers(
        vcf_path=vcf_path,
        target=target,
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
        start=int(target["region_start"]),
        end=int(target["region_end"]),
        phenotype_columns=PHENOTYPE_COLUMNS,
        marker_columns=[c for c in marker_df.columns if c != "SampleID"],
        marker_positions=marker_positions,
        expected_direction="increases_trait",
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
    )
    return db_dir


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    phenotype_xlsx = Path(args.phenotype_xlsx)
    vcf_dir = Path(args.vcf_dir)
    output_root = Path(args.output_root)
    intermediate_root = Path(args.intermediate_root)
    phenotype_df = load_watkins_growth_habit_phenotypes(phenotype_xlsx)
    selected_targets = args.target or sorted(TARGETS)

    print(f"[INFO] Loaded Watkins growth-habit phenotypes: {len(phenotype_df)} samples")
    status_rows: List[Dict[str, object]] = []
    for target_id in selected_targets:
        target = TARGETS[target_id]
        vcf_path = vcf_dir / str(target["vcf"])
        print(f"[INFO] Preparing {target_id} from {vcf_path}")
        try:
            db_dir = build_target_database(
                target_id=target_id,
                target=target,
                vcf_path=vcf_path,
                phenotype_df=phenotype_df,
                phenotype_xlsx=phenotype_xlsx,
                output_root=output_root,
                intermediate_root=intermediate_root,
                min_haplotype_count=args.min_haplotype_count,
                max_missing_rate=args.max_missing_rate,
            )
        except ValueError as e:
            status_rows.append({"target_id": target_id, "status": "blocked", "reason": str(e)})
            print(f"[BLOCKED] {target_id}: {e}")
            continue

        variant_count = len(pd.read_csv(db_dir / "variant_info.csv"))
        sample_count = len(pd.read_csv(db_dir / "phenotype_data.csv"))
        haplotype_count = len(pd.read_csv(db_dir / "haplotype_data.csv"))
        status_rows.append({
            "target_id": target_id,
            "status": "built",
            "variant_count": variant_count,
            "sample_count": sample_count,
            "haplotype_count": haplotype_count,
            "database": str(db_dir),
        })
        print(
            f"[INFO] Built {target_id}: {variant_count} variants, "
            f"{haplotype_count} haplotypes, {sample_count} samples -> {db_dir}"
        )

    intermediate_root.mkdir(parents=True, exist_ok=True)
    with (intermediate_root / "prepare_status.json").open("w", encoding="utf-8") as f:
        json.dump(status_rows, f, ensure_ascii=False, indent=2)
    print(
        "[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 "
        "--target VRN-A1-remoteSNP --target VRN-B1-remoteSNP --target VRN-D1-remoteSNP "
        "--score-mode robust_discovery"
    )
    return 0 if any(row["status"] == "built" for row in status_rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
