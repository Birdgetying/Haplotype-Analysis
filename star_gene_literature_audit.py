#!/usr/bin/env python3
"""Audit whether scored haplotypes contain literature functional variants."""

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple


AUDIT_COLUMNS = [
    "record_type",
    "paper",
    "target_id",
    "variant_name",
    "chrom",
    "position",
    "marker_id",
    "expected_allele",
    "expected_markers",
    "variant_class",
    "source_note",
    "covered_in_database",
    "covered_in_marker_matrix",
    "segregating_in_current_samples",
    "allele_counts",
    "carrier_count",
    "top_scored_haplotype",
    "top_scored_haplotype_score",
    "top_haplotype_sample_count",
    "top_haplotype_carrier_count",
    "top_haplotype_contains_expected",
    "validation_status",
    "covered_marker_count",
    "expected_marker_count",
    "missing_markers",
]


def _read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def _json_dump(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def _safe_float(value: Any) -> Optional[float]:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _marker_id(row: Dict[str, Any]) -> str:
    marker = row.get("marker_id")
    if marker is not None and str(marker).strip():
        return str(marker).strip()
    chrom = str(row.get("chrom") or row.get("chr") or "").strip()
    pos = row.get("position")
    if chrom and pos not in (None, ""):
        return f"{chrom}_{pos}"
    return ""


def _split_haplotype(seq: Any) -> List[str]:
    text = str(seq or "").strip()
    if not text:
        return []
    return text.split("|")


def _allele_matches(observed: Any, expected: Any) -> bool:
    if expected is None or expected == "":
        return False
    observed_text = str(observed or "").strip()
    expected_text = str(expected).strip()
    if not observed_text:
        return False
    if observed_text == expected_text:
        return True
    observed_parts = {part.strip() for part in observed_text.replace("|", "/").split("/") if part.strip()}
    return expected_text in observed_parts


def _format_counts(counter: Counter) -> str:
    return ";".join(f"{allele}:{count}" for allele, count in sorted(counter.items()))


def _load_marker_context(database_path: Path) -> Dict[str, Any]:
    variant_rows = _read_csv_rows(database_path / "variant_info.csv")
    haplotype_rows = _read_csv_rows(database_path / "haplotype_data.csv")
    sample_rows = _read_csv_rows(database_path / "haplotype_samples.csv")

    marker_order = [_marker_id(row) for row in variant_rows]
    marker_order = [marker for marker in marker_order if marker]
    marker_to_index = {marker: i for i, marker in enumerate(marker_order)}

    hap_alleles: Dict[str, List[str]] = {}
    hap_counts: Counter = Counter()
    for row in haplotype_rows:
        hap_name = str(row.get("Hap_Name") or "").strip()
        if not hap_name:
            continue
        hap_alleles[hap_name] = _split_haplotype(row.get("Haplotype_Seq"))
        count_value = row.get("Count")
        try:
            hap_counts[hap_name] = int(float(count_value)) if count_value not in (None, "") else 0
        except (TypeError, ValueError):
            hap_counts[hap_name] = 0

    sample_haps: Dict[str, str] = {}
    sample_alleles: Dict[str, List[str]] = {}
    sample_hap_counts: Counter = Counter()
    for row in sample_rows:
        sample_id = str(row.get("SampleID") or "").strip()
        hap_name = str(row.get("Hap_Name") or "").strip()
        if not sample_id or not hap_name:
            continue
        sample_haps[sample_id] = hap_name
        sample_hap_counts[hap_name] += 1
        sample_alleles[sample_id] = _split_haplotype(row.get("Haplotype_Seq"))

    if sample_hap_counts:
        for hap_name, count in sample_hap_counts.items():
            hap_counts[hap_name] = count

    return {
        "variant_rows": variant_rows,
        "marker_order": marker_order,
        "marker_to_index": marker_to_index,
        "hap_alleles": hap_alleles,
        "hap_counts": hap_counts,
        "sample_haps": sample_haps,
        "sample_alleles": sample_alleles,
    }


def _pick_top_haplotype(result_path: Path) -> Tuple[Optional[str], Optional[float]]:
    score_path = result_path / "haplotype_scores.json"
    if not score_path.exists():
        return None, None
    try:
        with score_path.open("r", encoding="utf-8") as f:
            score_data = json.load(f)
    except Exception:
        return None, None

    best_hap = None
    best_score = None
    for trait_scores in score_data.values():
        if not isinstance(trait_scores, dict):
            continue
        per_haplotype = trait_scores.get("per_haplotype") or {}
        if not isinstance(per_haplotype, dict):
            continue
        for hap_name, values in per_haplotype.items():
            score = _safe_float((values or {}).get("total"))
            if score is None:
                continue
            if best_score is None or score > best_score:
                best_hap = str(hap_name)
                best_score = score
    return best_hap, best_score


def _allele_for_marker(alleles: Sequence[str], marker_to_index: Dict[str, int], marker_id: str) -> Optional[str]:
    index = marker_to_index.get(marker_id)
    if index is None or index >= len(alleles):
        return None
    return alleles[index]


def _variant_record(
    paper: Dict[str, Any],
    target: Dict[str, Any],
    variant: Dict[str, Any],
    context: Dict[str, Any],
    top_hap: Optional[str],
    top_score: Optional[float],
) -> Dict[str, Any]:
    paper_id = paper.get("paper_id") or paper.get("short_name") or ""
    target_id = target.get("target_id") or target.get("gene_or_locus") or ""
    marker_id = str(variant.get("marker_id") or "").strip()
    if not marker_id and variant.get("chrom") and variant.get("position") is not None:
        marker_id = f"{variant.get('chrom')}_{variant.get('position')}"
    expected = variant.get("expected_allele")
    covered = marker_id in context["marker_to_index"]

    allele_counts: Counter = Counter()
    carrier_count = 0
    top_carrier_count = 0
    if covered and expected not in (None, ""):
        for sample_id, alleles in context["sample_alleles"].items():
            observed = _allele_for_marker(alleles, context["marker_to_index"], marker_id)
            if observed is None:
                continue
            allele_counts[str(observed)] += 1
            if _allele_matches(observed, expected):
                carrier_count += 1
                if context["sample_haps"].get(sample_id) == top_hap:
                    top_carrier_count += 1

    top_hap_count = int(context["hap_counts"].get(top_hap, 0)) if top_hap else 0
    top_contains = False
    if covered and expected not in (None, "") and top_hap:
        top_allele = _allele_for_marker(
            context["hap_alleles"].get(top_hap, []),
            context["marker_to_index"],
            marker_id,
        )
        top_contains = _allele_matches(top_allele, expected)

    segregating = len(allele_counts) > 1
    if not covered:
        status = "missing_from_database"
    elif expected in (None, ""):
        status = "not_testable_no_expected_allele"
    elif not segregating:
        status = "present_not_segregating"
    elif top_contains:
        status = "matched_top_haplotype"
    else:
        status = "present_but_not_top"

    return {
        "record_type": "variant",
        "paper": paper_id,
        "target_id": target_id,
        "variant_name": variant.get("variant_name") or marker_id,
        "chrom": variant.get("chrom") or "",
        "position": variant.get("position") or "",
        "marker_id": marker_id,
        "expected_allele": expected or "",
        "expected_markers": "",
        "variant_class": variant.get("variant_class") or "",
        "source_note": variant.get("source_note") or "",
        "covered_in_database": covered,
        "covered_in_marker_matrix": covered,
        "segregating_in_current_samples": segregating,
        "allele_counts": _format_counts(allele_counts),
        "carrier_count": carrier_count,
        "top_scored_haplotype": top_hap or "",
        "top_scored_haplotype_score": top_score if top_score is not None else "",
        "top_haplotype_sample_count": top_hap_count,
        "top_haplotype_carrier_count": top_carrier_count,
        "top_haplotype_contains_expected": top_contains,
        "validation_status": status,
        "covered_marker_count": 1 if covered else 0,
        "expected_marker_count": 1,
        "missing_markers": "" if covered else marker_id,
    }


def _haplotype_record(
    paper: Dict[str, Any],
    target: Dict[str, Any],
    haplotype: Dict[str, Any],
    context: Dict[str, Any],
    top_hap: Optional[str],
    top_score: Optional[float],
) -> Dict[str, Any]:
    paper_id = paper.get("paper_id") or paper.get("short_name") or ""
    target_id = target.get("target_id") or target.get("gene_or_locus") or ""
    expected_markers = haplotype.get("expected_markers") or {}
    expected_markers = {str(marker): str(allele) for marker, allele in expected_markers.items()}
    missing = [marker for marker in expected_markers if marker not in context["marker_to_index"]]
    covered = not missing and bool(expected_markers)

    carrier_count = 0
    top_carrier_count = 0
    hap_counts: Counter = Counter()
    if covered:
        for sample_id, alleles in context["sample_alleles"].items():
            is_carrier = all(
                _allele_matches(
                    _allele_for_marker(alleles, context["marker_to_index"], marker),
                    expected,
                )
                for marker, expected in expected_markers.items()
            )
            if is_carrier:
                carrier_count += 1
                sample_hap = context["sample_haps"].get(sample_id)
                hap_counts[sample_hap] += 1
                if sample_hap == top_hap:
                    top_carrier_count += 1

    top_hap_count = int(context["hap_counts"].get(top_hap, 0)) if top_hap else 0
    top_contains = False
    if covered and top_hap:
        top_alleles = context["hap_alleles"].get(top_hap, [])
        top_contains = all(
            _allele_matches(
                _allele_for_marker(top_alleles, context["marker_to_index"], marker),
                expected,
            )
            for marker, expected in expected_markers.items()
        )

    if not expected_markers:
        status = "not_testable_no_expected_allele"
    elif not covered:
        status = "missing_from_database"
    elif carrier_count == 0:
        status = "haplotype_not_observed"
    elif top_contains:
        status = "matched_top_haplotype"
    else:
        status = "present_but_not_top"

    return {
        "record_type": "haplotype",
        "paper": paper_id,
        "target_id": target_id,
        "variant_name": haplotype.get("haplotype_name") or "",
        "chrom": haplotype.get("chrom") or "",
        "position": "",
        "marker_id": "",
        "expected_allele": "",
        "expected_markers": ";".join(f"{marker}={allele}" for marker, allele in expected_markers.items()),
        "variant_class": haplotype.get("variant_class") or "literature_haplotype",
        "source_note": haplotype.get("source_note") or "",
        "covered_in_database": covered,
        "covered_in_marker_matrix": covered,
        "segregating_in_current_samples": carrier_count > 0 and carrier_count < len(context["sample_alleles"]),
        "allele_counts": _format_counts(hap_counts),
        "carrier_count": carrier_count,
        "top_scored_haplotype": top_hap or "",
        "top_scored_haplotype_score": top_score if top_score is not None else "",
        "top_haplotype_sample_count": top_hap_count,
        "top_haplotype_carrier_count": top_carrier_count,
        "top_haplotype_contains_expected": top_contains,
        "validation_status": status,
        "covered_marker_count": len(expected_markers) - len(missing),
        "expected_marker_count": len(expected_markers),
        "missing_markers": ";".join(missing),
    }


def run_literature_audit(
    paper: Dict[str, Any],
    target: Dict[str, Any],
    database_path: Path,
    result_path: Path,
) -> List[Dict[str, Any]]:
    """Write literature functional-variant audit files for one target result."""
    variants = list(target.get("literature_variants") or [])
    haplotypes = list(target.get("literature_haplotypes") or [])
    if not variants and not haplotypes:
        return []

    database_path = Path(database_path)
    result_path = Path(result_path)
    result_path.mkdir(parents=True, exist_ok=True)

    context = _load_marker_context(database_path)
    top_hap, top_score = _pick_top_haplotype(result_path)
    records: List[Dict[str, Any]] = []
    for variant in variants:
        records.append(_variant_record(paper, target, variant, context, top_hap, top_score))
    for haplotype in haplotypes:
        records.append(_haplotype_record(paper, target, haplotype, context, top_hap, top_score))

    csv_path = result_path / "literature_variant_audit.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=AUDIT_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow(record)
    _json_dump(result_path / "literature_variant_audit.json", records)
    return records
