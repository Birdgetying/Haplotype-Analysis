#!/usr/bin/env python3
"""Download WheatOmics SNP micro-VCFs for VRN-A1/VRN-B1/VRN-D1 regions."""

import argparse
import io
import subprocess
import sys
from pathlib import Path

if sys.platform == "win32":
    try:
        if not isinstance(sys.stdout, io.TextIOWrapper) or sys.stdout.encoding != "utf-8":
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
        if not isinstance(sys.stderr, io.TextIOWrapper) or sys.stderr.encoding != "utf-8":
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
    except (ValueError, AttributeError):
        pass


REMOTE_WHEATOMICS_SNP_VCF = (
    "https://wheatomics.sdau.edu.cn/jbrowse-1.12.3-release/Chinese_Spring1.0/"
    "tracks/vcf/merge.SNP.Missing-unphasing.ID.ann.finalSID.allele2_retain."
    "hard_retain.InbreedingCoeff_retain.clean.anno.vcf.gz"
)

DEFAULT_OUTPUT_DIR = Path("external_data/wheat_nature_2024/vrn_remote")

TARGETS = {
    "VRN-A1": {
        "region": "chr5A:587409454-587425416",
        "output": "VRN-A1.wheatomics_snp.vcf.gz",
        "gene_id": "TraesCS5A01G391700",
    },
    "VRN-B1": {
        "region": "chr5B:573800883-573818070",
        "output": "VRN-B1.wheatomics_snp.vcf.gz",
        "gene_id": "TraesCS5B01G396600",
    },
    "VRN-D1": {
        "region": "chr5D:467174609-467186508",
        "output": "VRN-D1.wheatomics_snp.vcf.gz",
        "gene_id": "TraesCS5D01G401500",
    },
}


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Extract VRN-A1/VRN-B1/VRN-D1 SNP-only micro-VCFs from the WheatOmics "
            "1051-accession merged SNP VCF. These files do not include VRN INDEL/SV/CNV "
            "causal alleles unless those alleles are represented in the source SNP VCF."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--source-vcf", default=REMOTE_WHEATOMICS_SNP_VCF,
                        help="Remote or local WheatOmics merged SNP VCF")
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR),
                        help="Directory for downloaded micro-VCFs")
    parser.add_argument("--target", action="append", choices=sorted(TARGETS),
                        help="Target to download; repeatable. Defaults to all VRN targets.")
    parser.add_argument("--bcftools", default="bcftools",
                        help="bcftools executable. On Windows, default uses WSL Ubuntu.")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing micro-VCFs")
    return parser


def to_wsl_path(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        return str(resolved).replace("\\", "/")
    rel = str(resolved)[3:].replace("\\", "/")
    return f"/mnt/{drive}/{rel}"


def bcftools_cmd(args: list[str], bcftools: str) -> list[str]:
    if sys.platform == "win32" and bcftools == "bcftools":
        return ["wsl", "-d", "Ubuntu", "--", "bcftools", *args]
    return [bcftools, *args]


def run_command(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(args, check=True, capture_output=True, text=True, encoding="utf-8")


def output_arg(path: Path, bcftools: str) -> str:
    if sys.platform == "win32" and bcftools == "bcftools":
        return to_wsl_path(path)
    return str(path)


def download_target(target_id: str, source_vcf: str, output_dir: Path, bcftools: str, force: bool) -> None:
    target = TARGETS[target_id]
    output_path = output_dir / str(target["output"])
    output_dir.mkdir(parents=True, exist_ok=True)

    if output_path.exists() and output_path.stat().st_size > 0 and not force:
        print(f"{target_id}: exists, skip {output_path}")
    else:
        region = str(target["region"])
        print(f"{target_id}: downloading {region} -> {output_path}")
        cmd = bcftools_cmd(
            ["view", "-Oz", "-r", region, "-o", output_arg(output_path, bcftools), source_vcf],
            bcftools,
        )
        run_command(cmd)

    index_path = output_path.with_suffix(output_path.suffix + ".csi")
    if force or not index_path.exists():
        print(f"{target_id}: indexing {output_path}")
        run_command(bcftools_cmd(["index", "-f", output_arg(output_path, bcftools)], bcftools))

    verify_path = output_arg(output_path, bcftools)
    records = run_command(bcftools_cmd(["view", "-H", verify_path], bcftools)).stdout.count("\n")
    samples = len([line for line in run_command(bcftools_cmd(["query", "-l", verify_path], bcftools)).stdout.splitlines() if line])
    print(f"{target_id}: records={records} samples={samples} gene_id={target['gene_id']}")


def main() -> int:
    args = build_arg_parser().parse_args()
    output_dir = Path(args.output_dir)
    targets = args.target or sorted(TARGETS)
    for target_id in targets:
        download_target(target_id, args.source_vcf, output_dir, args.bcftools, args.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
