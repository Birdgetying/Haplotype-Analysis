#!/usr/bin/env python3
"""Prepare Rht-B1b/Rht-D1b functional-SNP positive-control databases."""

import argparse
import gzip
import io
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

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


REMOTE_WHEATOMICS_SNP_VCF = (
    "https://wheatomics.sdau.edu.cn/jbrowse-1.12.3-release/Chinese_Spring1.0/"
    "tracks/vcf/merge.SNP.Missing-unphasing.ID.ann.finalSID.allele2_retain."
    "hard_retain.InbreedingCoeff_retain.clean.anno.vcf.gz"
)
DEFAULT_PHENOTYPE_XLSX = Path(
    "external_data/wheat_nature_2024/watseq/"
    "Watkins_Collection_WGIN_WISP_DFW_watseq_phenotype_data_JIC.xlsx"
)
DEFAULT_OUTPUT_ROOT = Path("star_gene_database/wheat_nature_2024")
DEFAULT_INTERMEDIATE_ROOT = Path("external_data/wheat_nature_2024/rht1_functional")

PHENOTYPE_COLUMNS = ["PlantHeight_CFLN06", "GrowthHabit_CFLN06"]

TARGETS = {
    "Rht-B1b": {
        "target_id": "Rht-B1b",
        "gene_or_locus": "Rht-B1",
        "traes_id": "TraesCS4B02G043100",
        "wheatomics_id": "TraesCS4B01G043100",
        "chrom": "4B",
        "region_chrom": "chr4B",
        "region_start": 30_861_268,
        "region_end": 30_863_723,
        "strand": "+",
        "marker_id": "chr4B_30861571",
        "position": 30_861_571,
        "ref": "C",
        "alt": "T",
        "cdna_change": "c.190C>T",
        "protein_change": "p.Gln64*",
    },
    "Rht-D1b": {
        "target_id": "Rht-D1b",
        "gene_or_locus": "Rht-D1",
        "traes_id": "TraesCS4D02G040400",
        "wheatomics_id": "TraesCS4D01G040400",
        "chrom": "4D",
        "region_chrom": "chr4D",
        "region_start": 18_781_062,
        "region_end": 18_782_933,
        "strand": "+",
        "marker_id": "chr4D_18781242",
        "position": 18_781_242,
        "ref": "G",
        "alt": "T",
        "cdna_change": "c.181G>T",
        "protein_change": "p.Glu61*",
    },
}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build Rht-B1b/Rht-D1b functional SNP positive-control databases "
            "from WheatOmics SNP VCF calls and Watkins CFLN06 plant-height phenotypes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--vcf", default=REMOTE_WHEATOMICS_SNP_VCF,
                        help="Local/remote VCF containing Rht-B1b and Rht-D1b functional SNPs")
    parser.add_argument("--phenotype-xlsx", default=str(DEFAULT_PHENOTYPE_XLSX),
                        help="Watkins JIC phenotype workbook")
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT),
                        help="Output star-gene database root")
    parser.add_argument("--intermediate-root", default=str(DEFAULT_INTERMEDIATE_ROOT),
                        help="Directory for marker matrix, phenotype TSV, and micro-VCFs")
    parser.add_argument("--target", action="append", choices=sorted(TARGETS),
                        help="Target to prepare; repeatable. Defaults to both Rht-B1b and Rht-D1b.")
    parser.add_argument("--min-haplotype-count", type=int, default=1,
                        help="Minimum sample count for retained haplotypes")
    parser.add_argument("--bcftools", default="bcftools",
                        help="bcftools executable; on Windows this can be used through WSL.")
    return parser


def normalize_chrom(value: object) -> str:
    text = str(value).strip()
    return text[3:] if text.lower().startswith("chr") else text


def open_text_maybe_gzip(path_or_url: str):
    if str(path_or_url).startswith(("http://", "https://")):
        raise ValueError("remote VCF must be region-extracted with bcftools first")
    path = Path(path_or_url)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def load_watkins_height_phenotypes(phenotype_xlsx: Path) -> pd.DataFrame:
    if not phenotype_xlsx.exists():
        raise FileNotFoundError(f"phenotype workbook is missing: {phenotype_xlsx}")

    df = pd.read_excel(phenotype_xlsx, sheet_name="WGIN_Watkins_JIC_CFLN06")
    required = ["StoreCode", "PH_M_cm-CFLN06"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"WGIN_Watkins_JIC_CFLN06 missing columns: {', '.join(missing)}")

    out = df[["StoreCode", "PH_M_cm-CFLN06"]].copy()
    out = out.rename(columns={"StoreCode": "SampleID", "PH_M_cm-CFLN06": "PlantHeight_CFLN06"})
    if "GrwHabit_E_sw-CFLN06" in df.columns:
        out["GrowthHabit_CFLN06"] = df["GrwHabit_E_sw-CFLN06"]
    else:
        out["GrowthHabit_CFLN06"] = ""
    out["SampleID"] = out["SampleID"].astype(str).str.strip()
    out["PlantHeight_CFLN06"] = pd.to_numeric(out["PlantHeight_CFLN06"], errors="coerce")
    out["GrowthHabit_CFLN06"] = out["GrowthHabit_CFLN06"].astype(str).str.strip()
    out = out.dropna(subset=["SampleID", "PlantHeight_CFLN06"])
    out = out[out["SampleID"] != ""]
    out = out.groupby("SampleID", as_index=False).agg({
        "PlantHeight_CFLN06": "mean",
        "GrowthHabit_CFLN06": "first",
    })
    if out.empty:
        raise ValueError("no complete CFLN06 plant-height phenotypes after filtering")
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


def is_remote_vcf(vcf: str) -> bool:
    return str(vcf).startswith(("http://", "https://"))


def bcftools_command_for_current_platform(bcftools: str, region: str, vcf: str) -> List[str]:
    if sys.platform == "win32" and bcftools == "bcftools":
        return ["wsl", "-d", "Ubuntu", "--", "bcftools", "view", "-H", "-r", region, vcf]
    return [bcftools, "view", "-H", "-r", region, vcf]


def fetch_region_records(vcf: str, target: Dict[str, object], bcftools: str) -> List[str]:
    region = f"{target['region_chrom']}:{target['position']}-{target['position']}"
    if is_remote_vcf(vcf):
        cmd = bcftools_command_for_current_platform(bcftools, region, vcf)
        completed = subprocess.run(cmd, check=True, capture_output=True, text=True, encoding="utf-8")
        return [line for line in completed.stdout.splitlines() if line and not line.startswith("#")]

    records: List[str] = []
    target_chrom = normalize_chrom(target["region_chrom"])
    target_pos = int(target["position"])
    with open_text_maybe_gzip(vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 8:
                continue
            if normalize_chrom(parts[0]) == target_chrom and int(parts[1]) == target_pos:
                records.append(line.rstrip("\n\r"))
    return records


def read_vcf_header_samples(vcf: str, bcftools: str) -> List[str]:
    if is_remote_vcf(vcf):
        cmd = (["wsl", "-d", "Ubuntu", "--", "bcftools", "view", "-h", vcf]
               if sys.platform == "win32" and bcftools == "bcftools"
               else [bcftools, "view", "-h", vcf])
        completed = subprocess.run(cmd, check=True, capture_output=True, text=True, encoding="utf-8")
        for line in completed.stdout.splitlines():
            if line.startswith("#CHROM"):
                return line.rstrip("\n\r").split("\t")[9:]
        raise ValueError("VCF header contains no #CHROM line")

    with open_text_maybe_gzip(vcf) as f:
        for line in f:
            if line.startswith("#CHROM"):
                return line.rstrip("\n\r").split("\t")[9:]
    raise ValueError(f"VCF header contains no #CHROM line: {vcf}")


def build_marker_matrix_for_target(
    vcf: str,
    target: Dict[str, object],
    phenotype_samples: Sequence[str],
    bcftools: str,
) -> Tuple[pd.DataFrame, Dict[str, object]]:
    records = fetch_region_records(vcf, target, bcftools)
    if not records:
        raise ValueError(f"{target['target_id']}: functional SNP not found in VCF at {target['region_chrom']}:{target['position']}")

    sample_names = read_vcf_header_samples(vcf, bcftools)
    phenotype_sample_set = {str(sample) for sample in phenotype_samples}
    sample_indices = [
        i for i, sample in enumerate(sample_names, start=9)
        if sample in phenotype_sample_set
    ]
    if not sample_indices:
        raise ValueError(f"{target['target_id']}: no VCF samples overlap phenotype samples")

    record_parts = None
    for record in records:
        parts = record.split("\t")
        if len(parts) >= 10 and parts[2] == target["marker_id"]:
            record_parts = parts
            break
    if record_parts is None:
        record_parts = records[0].split("\t")

    ref = record_parts[3]
    alt = record_parts[4]
    fmt = record_parts[8]
    calls: Dict[str, str] = {}
    for sample_idx in sample_indices:
        sample = sample_names[sample_idx - 9]
        state = genotype_to_state(record_parts[sample_idx], fmt, ref, alt)
        if state is not None:
            calls[sample] = state
    if len(set(calls.values())) < 2:
        raise ValueError(f"{target['target_id']}: functional SNP is not segregating in phenotype-overlap samples")

    rows = [{"SampleID": sample, str(target["marker_id"]): calls[sample]} for sample in sorted(calls)]
    metadata = {
        "marker_id": str(target["marker_id"]),
        "position": int(record_parts[1]),
        "ref": ref,
        "alt": alt,
        "record_id": record_parts[2],
        "info": record_parts[7],
        "missing_rate": 1 - (len(calls) / len(sample_indices)),
    }
    return pd.DataFrame(rows), metadata


def update_metadata(
    db_dir: Path,
    target: Dict[str, object],
    source_vcf: str,
    phenotype_xlsx: Path,
    marker_metadata: Dict[str, object],
) -> None:
    gene_info_path = db_dir / "gene_info.json"
    with gene_info_path.open("r", encoding="utf-8") as f:
        gene_info = json.load(f)
    gene_info.update({
        "gene_id": target["target_id"],
        "gene_symbol": target["target_id"],
        "gene_or_locus": target["gene_or_locus"],
        "traes_id": target["traes_id"],
        "wheatomics_id": target["wheatomics_id"],
        "chrom": target["chrom"],
        "start": target["region_start"],
        "end": target["region_end"],
        "gene_start": target["region_start"],
        "gene_end": target["region_end"],
        "strand": target["strand"],
        "expected_direction": "decreases_trait",
        "source": "wheatomics_remote_rht1_functional_snp_vcf",
        "source_vcf": source_vcf,
        "source_phenotype": str(phenotype_xlsx),
        "source_note": (
            f"{target['target_id']} functional stop-gained SNP "
            f"{target['marker_id']} {target['ref']}>{target['alt']} "
            f"({target['cdna_change']}, {target['protein_change']}) from WheatOmics merged SNP VCF."
        ),
    })
    with gene_info_path.open("w", encoding="utf-8") as f:
        json.dump(gene_info, f, ensure_ascii=False, indent=2)

    variant_path = db_dir / "variant_info.csv"
    variant_df = pd.read_csv(variant_path)
    variant_df["position"] = int(marker_metadata["position"])
    variant_df["ref"] = marker_metadata["ref"]
    variant_df["alt"] = marker_metadata["alt"]
    variant_df["len_diff"] = 0
    variant_df["is_sv"] = False
    variant_df["missing_rate"] = float(marker_metadata["missing_rate"])
    variant_df["annotation"] = "stop_gained"
    variant_df["marker_id"] = target["marker_id"]
    variant_df["cdna_change"] = target["cdna_change"]
    variant_df["protein_change"] = target["protein_change"]
    variant_df["impact"] = "HIGH"
    variant_df.to_csv(variant_path, index=False)


def build_target_database(
    target_id: str,
    target: Dict[str, object],
    vcf: str,
    phenotype_df: pd.DataFrame,
    phenotype_xlsx: Path,
    output_root: Path,
    intermediate_root: Path,
    min_haplotype_count: int,
    bcftools: str,
) -> Path:
    marker_df, marker_metadata = build_marker_matrix_for_target(
        vcf=vcf,
        target=target,
        phenotype_samples=phenotype_df["SampleID"].tolist(),
        bcftools=bcftools,
    )
    phenotype_subset = phenotype_df[phenotype_df["SampleID"].isin(marker_df["SampleID"])].copy()
    phenotype_subset = phenotype_subset.sort_values("SampleID")
    marker_df = marker_df[marker_df["SampleID"].isin(phenotype_subset["SampleID"])].sort_values("SampleID")

    target_intermediate = intermediate_root / target_id
    target_intermediate.mkdir(parents=True, exist_ok=True)
    marker_output = target_intermediate / "marker_matrix.tsv"
    phenotype_output = target_intermediate / "phenotype.tsv"
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
        marker_columns=[str(target["marker_id"])],
        marker_positions={str(target["marker_id"]): int(target["position"])},
        expected_direction="decreases_trait",
        sample_column="SampleID",
        min_haplotype_count=min_haplotype_count,
    )
    update_metadata(
        db_dir=db_dir,
        target=target,
        source_vcf=vcf,
        phenotype_xlsx=phenotype_xlsx,
        marker_metadata=marker_metadata,
    )
    return db_dir


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    phenotype_xlsx = Path(args.phenotype_xlsx)
    output_root = Path(args.output_root)
    intermediate_root = Path(args.intermediate_root)
    phenotype_df = load_watkins_height_phenotypes(phenotype_xlsx)
    selected_targets = args.target or sorted(TARGETS)

    print(f"[INFO] Loaded CFLN06 plant-height phenotypes: {len(phenotype_df)} samples")
    status_rows: List[Dict[str, object]] = []
    for target_id in selected_targets:
        target = TARGETS[target_id]
        print(f"[INFO] Preparing {target_id} from {args.vcf}")
        try:
            db_dir = build_target_database(
                target_id=target_id,
                target=target,
                vcf=str(args.vcf),
                phenotype_df=phenotype_df,
                phenotype_xlsx=phenotype_xlsx,
                output_root=output_root,
                intermediate_root=intermediate_root,
                min_haplotype_count=args.min_haplotype_count,
                bcftools=args.bcftools,
            )
        except ValueError as e:
            status_rows.append({
                "target_id": target_id,
                "status": "blocked",
                "reason": str(e),
            })
            print(f"[BLOCKED] {target_id}: {e}")
            continue

        variant_count = len(pd.read_csv(db_dir / "variant_info.csv"))
        sample_count = len(pd.read_csv(db_dir / "phenotype_data.csv"))
        status_rows.append({
            "target_id": target_id,
            "status": "built",
            "variant_count": variant_count,
            "sample_count": sample_count,
            "database": str(db_dir),
        })
        print(f"[INFO] Built {target_id}: {variant_count} variants, {sample_count} samples -> {db_dir}")

    intermediate_root.mkdir(parents=True, exist_ok=True)
    with (intermediate_root / "prepare_status.json").open("w", encoding="utf-8") as f:
        json.dump(status_rows, f, ensure_ascii=False, indent=2)
    print("[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Rht-B1b --target Rht-D1b")
    return 0 if any(row["status"] == "built" for row in status_rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
