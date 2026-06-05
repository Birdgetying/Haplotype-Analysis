#!/usr/bin/env python3
"""Lightweight data-source metadata for star-gene validation.

This module intentionally separates download instructions from execution. Large
paper datasets remain manual by default; small direct files can be downloaded by
copying the generated PowerShell commands.
"""

import csv
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence


PROJECT_ROOT = Path(__file__).resolve().parent


@dataclass(frozen=True)
class DataFile:
    paper_id: str
    short_name: str
    key: str
    label: str
    source: str
    local_path: str
    size_hint: str = "unknown"
    is_large: bool = False
    default_action: str = "instruction_only"
    notes: str = ""

    @property
    def absolute_path(self) -> Path:
        path = Path(self.local_path)
        if path.is_absolute():
            return path
        return PROJECT_ROOT / path

    @property
    def normalized_local_path(self) -> str:
        return self.local_path.replace("\\", "/")


DATA_FILES: List[DataFile] = [
    DataFile(
        paper_id="rice_science_2024",
        short_name="rice2024",
        key="rice2024_figshare_genotype_matrix",
        label="18K-rice NAM genotype matrix",
        source="https://doi.org/10.6084/m9.figshare.19166475",
        local_path="external_data/rice_science_2024/figshare_19166475/NAM_variations",
        size_hint="12.64 GB",
        is_large=True,
        default_action="manual_large_download",
        notes="Figshare dataset NAM_variations. Inspect file names and matrix format after download.",
    ),
    DataFile(
        paper_id="rice_science_2024",
        short_name="rice2024",
        key="rice2024_parent_genomes",
        label="Rice parent genomes",
        source="https://doi.org/10.6084/m9.figshare.18972851",
        local_path="external_data/rice_science_2024/figshare_18972851",
        size_hint="large",
        is_large=True,
        default_action="manual_large_download",
        notes="Optional context; not required for the lightweight positive-control pass.",
    ),
    DataFile(
        paper_id="maize_natgenet_2019",
        short_name="maize2019",
        key="maize2019_sv_382254",
        label="MaizeGo SV.382254 structural variation map",
        source="http://www.maizego.org/download/all.sv/SV.382254.zip",
        local_path="external_data/maize_natgenet_2019/maizego/SV.382254.zip",
        size_hint="31.37 MB",
        is_large=False,
        default_action="small_direct_download",
        notes="Candidate source for qHKW1/ZmBAM1d SV marker lookup; inspect contents before analysis.",
    ),
    DataFile(
        paper_id="maize_natgenet_2019",
        short_name="maize2019",
        key="maize2019_gwas_psv_21081",
        label="MaizeGo GWAS.pSV.21081",
        source="http://www.maizego.org/download/all.sv/GWAS.pSV.21081.zip",
        local_path="external_data/maize_natgenet_2019/maizego/GWAS.pSV.21081.zip",
        size_hint="2.96 MB",
        is_large=False,
        default_action="small_direct_download",
        notes="Paper-related pSV/GWAS resource; may identify qHKW1 association records.",
    ),
    DataFile(
        paper_id="maize_natgenet_2019",
        short_name="maize2019",
        key="maize2019_psv_80614",
        label="MaizeGo pSV.80614",
        source="http://www.maizego.org/download/all.sv/pSV.80614.zip",
        local_path="external_data/maize_natgenet_2019/maizego/pSV.80614.zip",
        size_hint="8.63 MB",
        is_large=False,
        default_action="small_direct_download",
        notes="Paper-related pSV resource; inspect whether it contains sample-level states.",
    ),
    DataFile(
        paper_id="maize_natgenet_2019",
        short_name="maize2019",
        key="maize2019_agronomic_blup",
        label="MaizeGo agronomic BLUP traits",
        source="http://www.maizego.org/download/blup_traits_final.csv",
        local_path="external_data/maize_natgenet_2019/maizego/blup_traits_final.csv",
        size_hint="117 KB",
        is_large=False,
        default_action="small_direct_download",
        notes="Phenotype table for the MaizeGo association panel; confirm HKW column availability after download.",
    ),
    DataFile(
        paper_id="wheat_nature_2024",
        short_name="wheat2024",
        key="wheat2024_wwwg2b_portal",
        label="WWWG2B breeding portal",
        source="https://wwwg2b.com/dataAvailable",
        local_path="external_data/wheat_nature_2024/watseq",
        size_hint="portal",
        is_large=False,
        default_action="instruction_only",
        notes="Use portal metadata to obtain exact phenotype and target-region VCF object names.",
    ),
    DataFile(
        paper_id="wheat_nature_2024",
        short_name="wheat2024",
        key="wheat2024_earlham_watseq",
        label="Earlham WatSeq OpenData path",
        source="https://opendata.earlham.ac.uk/wheat/under_license/toronto/WatSeq_2023-09-15_landrace_modern_Variation_Data/",
        local_path="external_data/wheat_nature_2024/watseq",
        size_hint="mixed",
        is_large=False,
        default_action="instruction_only",
        notes="Paper confirms this path, but directory listing/object names were not available from local requests.",
    ),
    DataFile(
        paper_id="wheat_nature_2024",
        short_name="wheat2024",
        key="wheat2024_ngdc_gsa",
        label="NGDC/GSA raw WGS projects",
        source="PRJCA019636; CRA012590",
        local_path="external_data/wheat_nature_2024/raw",
        size_hint="large",
        is_large=True,
        default_action="manual_large_download",
        notes="Raw WGS is not needed for the lightweight positive-control pass.",
    ),
]


def _matches_paper(data_file: DataFile, paper: Optional[str]) -> bool:
    if not paper:
        return True
    value = paper.lower()
    return value in {
        data_file.paper_id.lower(),
        data_file.short_name.lower(),
    }


def iter_data_files(paper: Optional[str] = None) -> Iterable[DataFile]:
    matches = [data_file for data_file in DATA_FILES if _matches_paper(data_file, paper)]
    if paper and not matches:
        raise ValueError(f"No data files matched paper filter: {paper}")
    return iter(matches)


def _powershell_download_command(data_file: DataFile) -> str:
    local_path = data_file.normalized_local_path
    parent = str(Path(local_path).parent).replace("\\", "/")
    return (
        f"New-Item -ItemType Directory -Force -Path '{parent}' | Out-Null; "
        f"Invoke-WebRequest -Uri '{data_file.source}' -OutFile '{local_path}'"
    )


def build_download_commands(paper: Optional[str] = None, include_large: bool = False) -> List[str]:
    commands: List[str] = []
    for data_file in iter_data_files(paper=paper):
        if data_file.default_action != "small_direct_download":
            continue
        if data_file.is_large and not include_large:
            continue
        commands.append(_powershell_download_command(data_file))
    return commands


def summarize_downloads(paper: Optional[str] = None) -> str:
    lines = []
    for data_file in iter_data_files(paper=paper):
        lines.append(f"- {data_file.key}")
        lines.append(f"  paper: {data_file.short_name}")
        lines.append(f"  label: {data_file.label}")
        lines.append(f"  source: {data_file.source}")
        lines.append(f"  local_path: {data_file.normalized_local_path}")
        lines.append(f"  size_hint: {data_file.size_hint}")
        lines.append(f"  action: {data_file.default_action}")
        if data_file.notes:
            lines.append(f"  notes: {data_file.notes}")
    return "\n".join(lines)


def _read_table(path: Path) -> List[Dict[str, str]]:
    with open(path, "r", encoding="utf-8", newline="") as f:
        sample = f.read(4096)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        reader = csv.DictReader(f, dialect=dialect)
        return [dict(row) for row in reader]


def _write_csv(path: Path, fieldnames: Sequence[str], rows: Sequence[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _is_configured_missing_value(value: object, missing_values: Sequence[str]) -> bool:
    text = str(value).strip()
    if text in missing_values:
        return True
    try:
        numeric_value = float(text)
    except ValueError:
        return False

    for missing in missing_values:
        try:
            if numeric_value == float(str(missing).strip()):
                return True
        except ValueError:
            continue
    return False


def _infer_marker_len_diff(marker: str, ref: str, alt: str) -> int:
    match = re.search(r'_(insertion|deletion)_(\d+)(?:$|_)', marker)
    if match:
        length = int(match.group(2))
        return length if match.group(1) == "insertion" else -length
    return len(alt) - len(ref) if alt else 0


def build_database_from_marker_matrix(
    marker_matrix: Path,
    phenotype_table: Path,
    output_root: Path,
    target_id: str,
    chrom: str,
    start: int,
    end: int,
    phenotype_columns: Sequence[str],
    marker_columns: Optional[Sequence[str]] = None,
    marker_positions: Optional[Dict[str, int]] = None,
    expected_direction: str = "unknown",
    sample_column: str = "SampleID",
    min_haplotype_count: int = 1,
    phenotype_missing_values: Optional[Sequence[str]] = None,
) -> Path:
    """Convert a small sample-by-marker matrix into the precomputed DB format.

    The marker matrix must be a CSV/TSV with one sample row and one marker column
    per variant. This is intended for target-region or marker-level validation,
    not whole-genome matrices.
    """
    marker_matrix = Path(marker_matrix)
    phenotype_table = Path(phenotype_table)
    output_root = Path(output_root)

    marker_rows = _read_table(marker_matrix)
    phenotype_rows = _read_table(phenotype_table)
    if not marker_rows:
        raise ValueError(f"marker matrix is empty: {marker_matrix}")
    if not phenotype_rows:
        raise ValueError(f"phenotype table is empty: {phenotype_table}")
    if sample_column not in marker_rows[0]:
        raise ValueError(f"marker matrix missing sample column: {sample_column}")
    if sample_column not in phenotype_rows[0]:
        raise ValueError(f"phenotype table missing sample column: {sample_column}")

    if marker_columns is None:
        marker_columns = [c for c in marker_rows[0].keys() if c != sample_column]
    marker_columns = list(marker_columns)
    if not marker_columns:
        raise ValueError("no marker columns selected")

    missing_markers = [c for c in marker_columns if c not in marker_rows[0]]
    if missing_markers:
        raise ValueError(f"marker matrix missing columns: {', '.join(missing_markers)}")

    marker_positions = marker_positions or {
        marker: start + i for i, marker in enumerate(marker_columns)
    }
    positions = [int(marker_positions[marker]) for marker in marker_columns]

    missing_values = tuple(str(value).strip() for value in phenotype_missing_values or [])
    phenotype_by_sample = {}
    for row in phenotype_rows:
        has_missing_phenotype = False
        for col in phenotype_columns:
            if col not in row:
                raise ValueError(f"phenotype table missing selected column: {col}")
            if _is_configured_missing_value(row.get(col, ""), missing_values):
                has_missing_phenotype = True
                break
        if not has_missing_phenotype:
            phenotype_by_sample[row[sample_column]] = row

    sample_hap_sequences: Dict[str, str] = {}
    alleles_by_marker = {marker: [] for marker in marker_columns}

    for row in marker_rows:
        sample_id = row[sample_column]
        if sample_id not in phenotype_by_sample:
            continue
        alleles = []
        skip = False
        for marker in marker_columns:
            allele = str(row.get(marker, "")).strip()
            if not allele or allele.upper() in {"N", "NA", "NAN", "."}:
                skip = True
                break
            alleles.append(allele)
            alleles_by_marker[marker].append(allele)
        if not skip:
            sample_hap_sequences[sample_id] = "|".join(alleles)

    if not sample_hap_sequences:
        raise ValueError("no overlapping samples with complete marker genotypes")

    seq_counts = Counter(sample_hap_sequences.values())
    retained_seqs = {
        seq for seq, count in seq_counts.items()
        if count >= min_haplotype_count
    }
    seq_to_hap = {
        seq: f"Hap{i + 1}"
        for i, (seq, _) in enumerate(seq_counts.most_common())
        if seq in retained_seqs
    }

    hap_sample_rows = []
    phenotype_data_rows = []
    for sample_id, seq in sample_hap_sequences.items():
        if seq not in seq_to_hap:
            continue
        hap_name = seq_to_hap[seq]
        hap_sample_rows.append({
            "SampleID": sample_id,
            "Haplotype_Seq": seq,
            "Hap_Name": hap_name,
        })
        pheno_row = phenotype_by_sample[sample_id]
        merged = {
            "SampleID": sample_id,
            "Haplotype_Seq": seq,
            "Hap_Name": hap_name,
        }
        for col in phenotype_columns:
            if col not in pheno_row:
                raise ValueError(f"phenotype table missing selected column: {col}")
            merged[col] = pheno_row[col]
        phenotype_data_rows.append(merged)

    haplotype_rows = []
    for seq, hap_name in seq_to_hap.items():
        haplotype_rows.append({
            "Haplotype_Seq": seq,
            "Count": seq_counts[seq],
            "Hap_Name": hap_name,
            "Alleles": "|".join(seq.split("|")),
        })

    variant_rows = []
    for marker, pos in zip(marker_columns, positions):
        allele_counts = Counter(alleles_by_marker[marker])
        sorted_alleles = [a for a, _ in allele_counts.most_common()]
        ref = sorted_alleles[0] if sorted_alleles else ""
        alt = sorted_alleles[1] if len(sorted_alleles) > 1 else ""
        sorted_counts = [count for _, count in allele_counts.most_common()]
        maf = sorted_counts[1] / sum(sorted_counts) if len(sorted_counts) > 1 else 0.0
        len_diff = _infer_marker_len_diff(marker, ref, alt)
        variant_rows.append({
            "position": pos,
            "ref": ref,
            "alt": alt,
            "len_diff": len_diff,
            "is_sv": abs(len_diff) >= 50,
            "maf": round(maf, 6),
            "missing_rate": 0.0,
            "annotation": "SV" if abs(len_diff) >= 50 else "other",
            "marker_id": marker,
        })

    db_dir = output_root / target_id
    db_dir.mkdir(parents=True, exist_ok=True)
    gene_info = {
        "gene_id": target_id,
        "chrom": chrom,
        "start": int(start),
        "end": int(end),
        "gene_start": int(start),
        "gene_end": int(end),
        "strand": "+",
        "length": int(end) - int(start) + 1,
        "promoter_length": 0,
        "exons": [],
        "cds": [],
        "expected_direction": expected_direction,
        "source": "marker_matrix",
    }
    with open(db_dir / "gene_info.json", "w", encoding="utf-8") as f:
        json.dump(gene_info, f, ensure_ascii=False, indent=2)

    _write_csv(db_dir / "haplotype_data.csv",
               ["Haplotype_Seq", "Count", "Hap_Name", "Alleles"], haplotype_rows)
    _write_csv(db_dir / "haplotype_samples.csv",
               ["SampleID", "Haplotype_Seq", "Hap_Name"], hap_sample_rows)
    _write_csv(db_dir / "phenotype_data.csv",
               ["SampleID", "Haplotype_Seq", "Hap_Name", *phenotype_columns], phenotype_data_rows)
    _write_csv(db_dir / "variant_info.csv",
               ["position", "ref", "alt", "len_diff", "is_sv", "maf", "missing_rate", "annotation", "marker_id"],
               variant_rows)

    return db_dir
