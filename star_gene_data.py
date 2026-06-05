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
        key="maize2019_sv_386014_full",
        label="MaizeGo SV.386014 full structural variation genotype package",
        source="https://pan.baidu.com/s/10ieQpWGTEC805K4sI4RHOg",
        local_path="external_data/maize_natgenet_2019/maizego/SV.386014.zip",
        size_hint="102.42 MB page listing; 107,394,066 bytes in Baidu metadata",
        is_large=False,
        default_action="manual_baidu_download",
        notes=(
            "Required for the paper-genotype qHKW1/ZmBAM1d 8.9-kb indel check if "
            "the small direct MaizeGo matrices do not contain the exact marker. "
            "Baidu share metadata is public, but direct download can require web verification/login."
        ),
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


def _normal_chrom(value: object) -> str:
    text = str(value).strip()
    return text[3:] if text.lower().startswith("chr") else text


def _parse_maizego_sv_marker(marker: str) -> Optional[Dict[str, object]]:
    match = re.search(
        r'^(?P<chrom>chr[^_]+)_(?P<start>\d+)_(?P<end>\d+)_'
        r'(?P<variant_type>insertion|deletion)_(?P<length>\d+)(?:$|_)',
        marker,
    )
    if not match:
        return None
    start = int(match.group("start"))
    end = int(match.group("end"))
    return {
        "chrom": _normal_chrom(match.group("chrom")),
        "marker_start": min(start, end),
        "marker_end": max(start, end),
        "variant_type": match.group("variant_type"),
        "sv_length": int(match.group("length")),
    }


def _hmp_sample_start(header: Sequence[str]) -> int:
    if "QCcode" in header:
        return list(header).index("QCcode") + 1
    if "pos" in header:
        return list(header).index("pos") + 1
    raise ValueError("MaizeGo/HapMap matrix header must contain QCcode or pos")


def _valid_maizego_state(value: object) -> bool:
    text = str(value).strip()
    return bool(text) and text.upper() not in {"NA", "NAN", "."}


def _format_counts(values: Sequence[str]) -> str:
    counts = Counter(values)
    return ";".join(f"{state}:{count}" for state, count in counts.most_common())


def scan_maizego_sv_candidates(
    matrix_paths: Sequence[Path],
    chrom: str,
    window_start: int,
    window_end: int,
    length_min: int = 8500,
    length_max: int = 9500,
    variant_types: Optional[Sequence[str]] = None,
) -> List[Dict[str, object]]:
    """Scan MaizeGo/HapMap-style SV matrices for candidate paper markers.

    MaizeGo structural-variation files store marker IDs such as
    ``chr1_30450000_30458900_insertion_8900`` in the ``rs#`` column and one
    accession per following sample column. The scanner is intentionally strict:
    it only reports markers whose parsed chromosome, marker interval, and SV
    length match the requested paper window.
    """
    accepted_types = {str(v).lower() for v in variant_types or []}
    target_chrom = _normal_chrom(chrom)
    candidates: List[Dict[str, object]] = []

    for matrix_path in matrix_paths:
        matrix_path = Path(matrix_path)
        with open(matrix_path, "r", encoding="utf-8", newline="") as f:
            header_line = f.readline()
            if not header_line:
                continue
            header = header_line.rstrip("\n\r").split("\t")
            try:
                marker_idx = header.index("rs#")
                alleles_idx = header.index("alleles") if "alleles" in header else None
                chrom_idx = header.index("chrom") if "chrom" in header else None
                pos_idx = header.index("pos") if "pos" in header else None
                sample_start = _hmp_sample_start(header)
            except ValueError as e:
                raise ValueError(f"unsupported MaizeGo matrix header in {matrix_path}: {e}") from e

            for line_number, line in enumerate(f, start=2):
                parts = line.rstrip("\n\r").split("\t")
                if len(parts) <= max(marker_idx, sample_start - 1):
                    continue
                marker = parts[marker_idx]
                parsed = _parse_maizego_sv_marker(marker)
                if not parsed:
                    continue
                if _normal_chrom(parsed["chrom"]) != target_chrom:
                    continue
                if accepted_types and str(parsed["variant_type"]).lower() not in accepted_types:
                    continue
                if not (length_min <= int(parsed["sv_length"]) <= length_max):
                    continue
                marker_start = int(parsed["marker_start"])
                marker_end = int(parsed["marker_end"])
                if marker_end < window_start or marker_start > window_end:
                    continue

                states = [value.strip() for value in parts[sample_start:] if _valid_maizego_state(value)]
                pos = ""
                if pos_idx is not None and pos_idx < len(parts):
                    pos = parts[pos_idx]
                alleles = ""
                if alleles_idx is not None and alleles_idx < len(parts):
                    alleles = parts[alleles_idx]
                chrom_value = parsed["chrom"]
                if chrom_idx is not None and chrom_idx < len(parts) and parts[chrom_idx]:
                    chrom_value = _normal_chrom(parts[chrom_idx])
                candidates.append({
                    "source_path": str(matrix_path),
                    "line_number": line_number,
                    "marker": marker,
                    "alleles": alleles,
                    "chrom": chrom_value,
                    "pos": pos,
                    "marker_start": marker_start,
                    "marker_end": marker_end,
                    "variant_type": parsed["variant_type"],
                    "sv_length": parsed["sv_length"],
                    "valid_sample_count": len(states),
                    "state_count": len(set(states)),
                    "counts": _format_counts(states),
                })

    return sorted(candidates, key=lambda row: (
        str(row["source_path"]),
        int(row["line_number"]),
    ))


def extract_maizego_marker_matrix(
    matrix_path: Path,
    marker_id: str,
    output_path: Path,
) -> Path:
    """Transpose one MaizeGo/HapMap marker row into a sample-by-marker table."""
    matrix_path = Path(matrix_path)
    output_path = Path(output_path)

    with open(matrix_path, "r", encoding="utf-8", newline="") as f:
        header_line = f.readline()
        if not header_line:
            raise ValueError(f"matrix is empty: {matrix_path}")
        header = header_line.rstrip("\n\r").split("\t")
        try:
            marker_idx = header.index("rs#")
            sample_start = _hmp_sample_start(header)
        except ValueError as e:
            raise ValueError(f"unsupported MaizeGo matrix header in {matrix_path}: {e}") from e

        samples = header[sample_start:]
        selected_parts = None
        for line in f:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) > marker_idx and parts[marker_idx] == marker_id:
                selected_parts = parts
                break

    if selected_parts is None:
        raise ValueError(f"marker not found in {matrix_path}: {marker_id}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["SampleID", marker_id])
        for i, sample in enumerate(samples, start=sample_start):
            allele = selected_parts[i] if i < len(selected_parts) else ""
            writer.writerow([sample, allele])

    return output_path


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
