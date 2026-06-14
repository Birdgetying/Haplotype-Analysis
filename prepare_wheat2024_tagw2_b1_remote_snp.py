#!/usr/bin/env python3
"""Prepare a TaGW2-B1 promoter-SNP positive-control database."""

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

from prepare_wheat2024_tagw2_ad import load_watkins_tgw_phenotypes
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
DEFAULT_INTERMEDIATE_ROOT = Path("external_data/wheat_nature_2024/tagw2_b1_remote")

TARGET_ID = "TaGW2-B1-remoteSNP"
PHENOTYPE_COLUMNS = ["TGW_mean", "TGW_CFLN10", "TGW_CFLN14"]

TAGW2_B1 = {
    "gene_id": "TraesCS6B02G215300",
    "chrom": "6B",
    "region_chrom": "chr6B",
    "gene_start": 291_761_229,
    "gene_end": 291_778_752,
    "strand": "+",
    "cds_start_offset_1based": 169,
    "atg_position": 291_761_397,
    "literature_haplotype": "Hap-6B-1",
    "literature_source": "Qin et al. 2014 BMC Plant Biology Fig. 2B",
    "promoter_snps": [
        {
            "name": "TaGW2-6B_-1709",
            "offset": -1709,
            "position": 291_759_688,
            "marker_id": "chr6B_291759688",
            "expected_allele": "A",
        },
        {
            "name": "TaGW2-6B_-721",
            "offset": -721,
            "position": 291_760_676,
            "marker_id": "chr6B_291760676",
            "expected_allele": "G",
        },
        {
            "name": "TaGW2-6B_-83",
            "offset": -83,
            "position": 291_761_314,
            "marker_id": "chr6B_291761314",
            "expected_allele": "C",
        },
    ],
}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build a TaGW2-B1 promoter-SNP positive-control database from a local "
            "or remote WheatOmics merged SNP VCF and Watkins TGW phenotypes."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--vcf", default=REMOTE_WHEATOMICS_SNP_VCF,
                        help="Local or remote SNP VCF containing TaGW2-B1 promoter SNP calls")
    parser.add_argument("--phenotype-xlsx", default=str(DEFAULT_PHENOTYPE_XLSX),
                        help="Watkins JIC phenotype workbook")
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT),
                        help="Output star-gene database root")
    parser.add_argument("--intermediate-root", default=str(DEFAULT_INTERMEDIATE_ROOT),
                        help="Directory for marker matrix, phenotype TSV, and prepare status")
    parser.add_argument("--min-haplotype-count", type=int, default=3,
                        help="Minimum sample count for retained haplotypes")
    parser.add_argument("--bcftools", default="bcftools",
                        help="bcftools executable; on Windows this can be used through WSL")
    return parser


def is_remote_vcf(vcf: str) -> bool:
    return str(vcf).startswith(("http://", "https://"))


def open_text_maybe_gzip(path_or_url: str):
    if is_remote_vcf(path_or_url):
        raise ValueError("remote VCF records must be extracted with bcftools")
    path = Path(path_or_url)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


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


def bcftools_base_command(bcftools: str) -> List[str]:
    if sys.platform == "win32" and bcftools == "bcftools":
        return ["wsl", "-d", "Ubuntu", "--", "bcftools"]
    return [bcftools]


def read_vcf_header_samples(vcf: str, bcftools: str) -> List[str]:
    if is_remote_vcf(vcf):
        cmd = [*bcftools_base_command(bcftools), "view", "-h", vcf]
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


def fetch_vcf_records(vcf: str, region: str, bcftools: str) -> List[str]:
    if is_remote_vcf(vcf):
        cmd = [*bcftools_base_command(bcftools), "view", "-H", "-r", region, vcf]
        completed = subprocess.run(cmd, check=True, capture_output=True, text=True, encoding="utf-8")
        return [line for line in completed.stdout.splitlines() if line and not line.startswith("#")]

    records: List[str] = []
    chrom, span = region.split(":", 1)
    start_s, end_s = span.split("-", 1)
    start, end = int(start_s), int(end_s)
    with open_text_maybe_gzip(vcf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 10:
                continue
            if parts[0].replace("chr", "") != chrom.replace("chr", ""):
                continue
            pos = int(parts[1])
            if start <= pos <= end:
                records.append(line.rstrip("\n\r"))
    return records


def build_marker_matrix(
    vcf: str,
    phenotype_samples: Sequence[str],
    bcftools: str,
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, object]], List[str]]:
    region = f"{TAGW2_B1['region_chrom']}:{TAGW2_B1['promoter_snps'][0]['position']}-{TAGW2_B1['gene_end']}"
    records = fetch_vcf_records(vcf, region, bcftools)
    sample_names = read_vcf_header_samples(vcf, bcftools)
    phenotype_set = {str(sample) for sample in phenotype_samples}
    sample_indices = [
        i for i, sample in enumerate(sample_names, start=9)
        if sample in phenotype_set
    ]
    if not sample_indices:
        raise ValueError("TaGW2-B1: no VCF samples overlap Watkins TGW phenotypes")

    record_by_marker: Dict[str, List[str]] = {}
    for record in records:
        parts = record.split("\t")
        if len(parts) < 10:
            continue
        marker_id = parts[2] if parts[2] and parts[2] != "." else f"{parts[0]}_{parts[1]}"
        record_by_marker[marker_id] = parts

    marker_metadata: Dict[str, Dict[str, object]] = {}
    calls_by_marker: Dict[str, Dict[str, str]] = {}
    missing_markers: List[str] = []
    for snp in TAGW2_B1["promoter_snps"]:
        marker_id = str(snp["marker_id"])
        parts = record_by_marker.get(marker_id)
        if parts is None:
            missing_markers.append(marker_id)
            continue
        ref, alt, fmt = parts[3], parts[4], parts[8]
        calls: Dict[str, str] = {}
        for sample_idx in sample_indices:
            sample = sample_names[sample_idx - 9]
            state = genotype_to_state(parts[sample_idx], fmt, ref, alt)
            if state is not None:
                calls[sample] = state
        if len(set(calls.values())) < 2:
            missing_markers.append(f"{marker_id}:not_segregating")
            continue
        calls_by_marker[marker_id] = calls
        marker_metadata[marker_id] = {
            "position": int(parts[1]),
            "ref": ref,
            "alt": alt,
            "record_id": parts[2],
            "missing_rate": 1 - (len(calls) / len(sample_indices)),
        }

    if missing_markers:
        raise ValueError(
            "TaGW2-B1: literature promoter SNPs are missing or not segregating: "
            + ", ".join(missing_markers)
        )

    complete_samples = sorted(set.intersection(*(set(calls) for calls in calls_by_marker.values())))
    if len(complete_samples) < 2:
        raise ValueError("TaGW2-B1: fewer than two complete samples across literature promoter SNPs")

    rows = []
    marker_order = [str(snp["marker_id"]) for snp in TAGW2_B1["promoter_snps"]]
    for sample in complete_samples:
        row = {"SampleID": sample}
        for marker in marker_order:
            row[marker] = calls_by_marker[marker][sample]
        rows.append(row)
    return pd.DataFrame(rows), marker_metadata, marker_order


def update_gene_and_variant_metadata(
    db_dir: Path,
    source_vcf: str,
    phenotype_xlsx: Path,
    marker_metadata: Dict[str, Dict[str, object]],
) -> None:
    gene_info_path = db_dir / "gene_info.json"
    with gene_info_path.open("r", encoding="utf-8") as f:
        gene_info = json.load(f)
    gene_info.update({
        "gene_id": TARGET_ID,
        "gene_symbol": "TaGW2-B1",
        "traes_id": TAGW2_B1["gene_id"],
        "chrom": TAGW2_B1["chrom"],
        "start": TAGW2_B1["promoter_snps"][0]["position"],
        "end": TAGW2_B1["gene_end"],
        "gene_start": TAGW2_B1["gene_start"],
        "gene_end": TAGW2_B1["gene_end"],
        "strand": TAGW2_B1["strand"],
        "cds_start_offset_1based": TAGW2_B1["cds_start_offset_1based"],
        "atg_position": TAGW2_B1["atg_position"],
        "literature_haplotype": TAGW2_B1["literature_haplotype"],
        "source": "wheatomics_remote_tagw2_b1_promoter_snp_vcf",
        "source_vcf": source_vcf,
        "source_phenotype": str(phenotype_xlsx),
        "source_note": (
            "TaGW2-B1 natural promoter haplotype SNPs from Qin et al. 2014 Fig. 2B; "
            "expected Hap-6B-1 pattern is A/G/C at -1709/-721/-83."
        ),
    })
    with gene_info_path.open("w", encoding="utf-8") as f:
        json.dump(gene_info, f, ensure_ascii=False, indent=2)

    variant_path = db_dir / "variant_info.csv"
    variant_df = pd.read_csv(variant_path)
    for idx, row in variant_df.iterrows():
        marker_id = str(row["marker_id"])
        meta = marker_metadata[marker_id]
        variant_df.at[idx, "position"] = int(meta["position"])
        variant_df.at[idx, "ref"] = meta["ref"]
        variant_df.at[idx, "alt"] = meta["alt"]
        variant_df.at[idx, "missing_rate"] = float(meta["missing_rate"])
        variant_df.at[idx, "annotation"] = "promoter_snp"
    variant_df.to_csv(variant_path, index=False)


def build_database(
    vcf: str,
    phenotype_df: pd.DataFrame,
    phenotype_xlsx: Path,
    output_root: Path,
    intermediate_root: Path,
    min_haplotype_count: int,
    bcftools: str,
) -> Path:
    marker_df, marker_metadata, marker_order = build_marker_matrix(
        vcf=vcf,
        phenotype_samples=phenotype_df["SampleID"].tolist(),
        bcftools=bcftools,
    )
    phenotype_subset = phenotype_df[phenotype_df["SampleID"].isin(marker_df["SampleID"])].copy()
    phenotype_subset = phenotype_subset.sort_values("SampleID")
    marker_df = marker_df[marker_df["SampleID"].isin(phenotype_subset["SampleID"])].sort_values("SampleID")

    target_intermediate = intermediate_root / TARGET_ID
    target_intermediate.mkdir(parents=True, exist_ok=True)
    marker_output = target_intermediate / "marker_matrix.tsv"
    phenotype_output = target_intermediate / "phenotype.tsv"
    marker_df.to_csv(marker_output, sep="\t", index=False)
    phenotype_subset[["SampleID", *PHENOTYPE_COLUMNS]].to_csv(phenotype_output, sep="\t", index=False)

    marker_positions = {marker: int(marker_metadata[marker]["position"]) for marker in marker_order}
    db_dir = build_database_from_marker_matrix(
        marker_matrix=marker_output,
        phenotype_table=phenotype_output,
        output_root=output_root,
        target_id=TARGET_ID,
        chrom=str(TAGW2_B1["chrom"]),
        start=int(TAGW2_B1["promoter_snps"][0]["position"]),
        end=int(TAGW2_B1["gene_end"]),
        phenotype_columns=PHENOTYPE_COLUMNS,
        marker_columns=marker_order,
        marker_positions=marker_positions,
        expected_direction="increases_trait",
        sample_column="SampleID",
        min_haplotype_count=min_haplotype_count,
    )
    update_gene_and_variant_metadata(
        db_dir=db_dir,
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
    intermediate_root.mkdir(parents=True, exist_ok=True)

    phenotype_df = load_watkins_tgw_phenotypes(phenotype_xlsx)
    status: Dict[str, object]
    print(f"[INFO] Loaded grain-weight phenotypes: {len(phenotype_df)} samples")
    print(f"[INFO] Preparing {TARGET_ID} from {args.vcf}")
    try:
        db_dir = build_database(
            vcf=str(args.vcf),
            phenotype_df=phenotype_df,
            phenotype_xlsx=phenotype_xlsx,
            output_root=output_root,
            intermediate_root=intermediate_root,
            min_haplotype_count=args.min_haplotype_count,
            bcftools=args.bcftools,
        )
    except ValueError as e:
        status = {
            "target_id": TARGET_ID,
            "status": "blocked",
            "reason": str(e),
            "atg_position": TAGW2_B1["atg_position"],
            "literature_promoter_snps": TAGW2_B1["promoter_snps"],
            "source_vcf": str(args.vcf),
        }
        with (intermediate_root / "prepare_status.json").open("w", encoding="utf-8") as f:
            json.dump([status], f, ensure_ascii=False, indent=2)
        print(f"[BLOCKED] {status['reason']}")
        return 1

    status = {
        "target_id": TARGET_ID,
        "status": "built",
        "variant_count": len(pd.read_csv(db_dir / "variant_info.csv")),
        "sample_count": len(pd.read_csv(db_dir / "phenotype_data.csv")),
        "database": str(db_dir),
        "atg_position": TAGW2_B1["atg_position"],
        "literature_promoter_snps": TAGW2_B1["promoter_snps"],
        "source_vcf": str(args.vcf),
    }
    with (intermediate_root / "prepare_status.json").open("w", encoding="utf-8") as f:
        json.dump([status], f, ensure_ascii=False, indent=2)
    print(f"[INFO] Built {TARGET_ID}: {status['variant_count']} variants, {status['sample_count']} samples -> {db_dir}")
    print(f"[NEXT] python run_star_gene_validation.py --run-analysis --paper wheat2024 --target {TARGET_ID}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
