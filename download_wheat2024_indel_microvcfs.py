#!/usr/bin/env python3
"""Download target-region INDEL micro-VCFs from public WheatOmics tracks."""

import argparse
import csv
import io
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set

if sys.platform == "win32":
    try:
        if not isinstance(sys.stdout, io.TextIOWrapper) or sys.stdout.encoding != "utf-8":
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
        if not isinstance(sys.stderr, io.TextIOWrapper) or sys.stderr.encoding != "utf-8":
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
    except (ValueError, AttributeError):
        pass

import pandas as pd


WHEATOMICS_VCF_BASE = (
    "https://wheatomics.sdau.edu.cn/jbrowse-1.12.3-release/"
    "Chinese_Spring1.0/tracks/vcf"
)

DEFAULT_OUTPUT_DIR = Path("external_data/wheat_nature_2024/indel_microvcfs")
DEFAULT_PHENOTYPE_XLSX = Path(
    "external_data/wheat_nature_2024/watseq/"
    "Watkins_Collection_WGIN_WISP_DFW_watseq_phenotype_data_JIC.xlsx"
)

SOURCES = {
    "WEC_INDEL": "WEC_INDEL_IWGSCv1.0.eff.vcf.gz",
    "GBS_INDEL": "GBS_filtered_Indels_IWGSCv1.0.eff.vcf.gz",
    "WildEmmer10WGS_INDEL": "WildEmmer10WGS_INDEL_eff.vcf.gz",
}

TARGETS = {
    "VRN-A1": "chr5A:587409454-587425416",
    "VRN-B1": "chr5B:573800883-573818070",
    "VRN-D1": "chr5D:467174609-467186508",
    "Rht-B1": "chr4B:30861268-30863723",
    "Rht-D1": "chr4D:18781062-18782933",
    "TaGW2-A1": "chr6A:237732651-237760058",
    "TaGW2-B1": "chr6B:291759689-291778752",
    "TaGW2-D1": "chr6D:175710228-175721507",
}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract small target-region INDEL VCFs from WheatOmics public JBrowse tracks. "
            "These tracks are not the 1051-sample WatSeq SNP source, so the script also "
            "reports sample overlap with the Watkins phenotype workbook."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Output directory")
    parser.add_argument("--phenotype-xlsx", default=str(DEFAULT_PHENOTYPE_XLSX),
                        help="Watkins phenotype workbook for sample-overlap reporting")
    parser.add_argument("--source", action="append", choices=sorted(SOURCES),
                        help="Source to extract; repeatable. Defaults to all sources.")
    parser.add_argument("--target", action="append", choices=sorted(TARGETS),
                        help="Target region to extract; repeatable. Defaults to all targets.")
    parser.add_argument("--bcftools", default="bcftools",
                        help="bcftools executable. On Windows, default uses WSL Ubuntu.")
    parser.add_argument("--force", action="store_true", help="Overwrite existing micro-VCFs")
    return parser


def to_wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        return str(resolved).replace("\\", "/")
    rel = str(resolved)[3:].replace("\\", "/")
    return f"/mnt/{drive}/{rel}"


def bcftools_cmd(args: List[str], bcftools: str) -> List[str]:
    if sys.platform == "win32" and bcftools == "bcftools":
        return ["wsl", "-d", "Ubuntu", "--", "bcftools", *args]
    return [bcftools, *args]


def output_arg(path: Path, bcftools: str) -> str:
    if sys.platform == "win32" and bcftools == "bcftools":
        return to_wsl_path(path)
    return str(path)


def run_command(args: List[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, check=True, capture_output=True, text=True, encoding="utf-8")


def load_watkins_samples(phenotype_xlsx: Path) -> Set[str]:
    if not phenotype_xlsx.exists():
        return set()
    samples: Set[str] = set()
    workbook = pd.ExcelFile(phenotype_xlsx)
    for sheet in workbook.sheet_names:
        try:
            df = pd.read_excel(workbook, sheet_name=sheet)
        except Exception:
            continue
        for col in ("StoreCode", "SampleID", "Accession"):
            if col in df.columns:
                samples.update(df[col].dropna().astype(str).str.strip())
    return {sample for sample in samples if sample}


def query_samples(vcf_path: str, bcftools: str) -> Set[str]:
    completed = run_command(bcftools_cmd(["query", "-l", vcf_path], bcftools))
    return {line.strip() for line in completed.stdout.splitlines() if line.strip()}


def count_records(vcf_path: str, bcftools: str) -> int:
    completed = run_command(bcftools_cmd(["view", "-H", vcf_path], bcftools))
    return sum(1 for line in completed.stdout.splitlines() if line.strip())


def extract_microvcf(
    source_id: str,
    target_id: str,
    output_dir: Path,
    bcftools: str,
    force: bool,
) -> Dict[str, object]:
    source_name = SOURCES[source_id]
    region = TARGETS[target_id]
    source_url = f"{WHEATOMICS_VCF_BASE}/{source_name}"
    safe_region = region.replace(":", "_").replace("-", "_")
    output_path = output_dir / source_id / f"{target_id}.{safe_region}.vcf.gz"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if force or not output_path.exists() or output_path.stat().st_size == 0:
        print(f"[INFO] Extracting {source_id} {target_id} {region}")
        run_command(bcftools_cmd([
            "view", "-Oz", "-r", region, "-o", output_arg(output_path, bcftools), source_url
        ], bcftools))
        run_command(bcftools_cmd(["index", "-f", output_arg(output_path, bcftools)], bcftools))
    else:
        print(f"[INFO] Existing micro-VCF: {output_path}")

    vcf_for_bcftools = output_arg(output_path, bcftools)
    records = count_records(vcf_for_bcftools, bcftools)
    samples = query_samples(vcf_for_bcftools, bcftools)
    return {
        "source_id": source_id,
        "source_url": source_url,
        "target_id": target_id,
        "region": region,
        "output_vcf": str(output_path),
        "record_count": records,
        "sample_count": len(samples),
        "samples": sorted(samples),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sources = args.source or sorted(SOURCES)
    targets = args.target or sorted(TARGETS)
    watkins_samples = load_watkins_samples(Path(args.phenotype_xlsx))

    rows: List[Dict[str, object]] = []
    for source_id in sources:
        for target_id in targets:
            row = extract_microvcf(source_id, target_id, output_dir, args.bcftools, args.force)
            overlap = sorted(set(row["samples"]) & watkins_samples)
            row["watkins_overlap_count"] = len(overlap)
            row["watkins_overlap_samples"] = ";".join(overlap[:20])
            rows.append(row)
            print(
                f"[INFO] {source_id}/{target_id}: records={row['record_count']} "
                f"samples={row['sample_count']} watkins_overlap={len(overlap)}"
            )

    summary_json = output_dir / "microvcf_status.json"
    summary_tsv = output_dir / "microvcf_status.tsv"
    summary_json.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    with summary_tsv.open("w", encoding="utf-8", newline="") as f:
        fieldnames = [
            "source_id", "target_id", "region", "record_count", "sample_count",
            "watkins_overlap_count", "watkins_overlap_samples", "source_url", "output_vcf",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})
    print(f"[INFO] Wrote {summary_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
