#!/usr/bin/env python3
"""Download target-region marker annotations from WheatOmics Zhai INDEL track."""

import argparse
import csv
import io
import json
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib import request

if sys.platform == "win32":
    try:
        if not isinstance(sys.stdout, io.TextIOWrapper) or sys.stdout.encoding != "utf-8":
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
        if not isinstance(sys.stderr, io.TextIOWrapper) or sys.stderr.encoding != "utf-8":
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
    except (ValueError, AttributeError):
        pass


TRACK_BASE = (
    "https://wheatomics.sdau.edu.cn/jbrowse-1.12.3-release/"
    "Chinese_Spring1.0/tracks/Indel_marker_from_zhai"
)
DEFAULT_OUTPUT_DIR = Path("external_data/wheat_nature_2024/zhai_indel_markers")
DEFAULT_PADDING_BP = 1_000_000

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
            "Fetch target-region marker annotations from the small WheatOmics "
            "Indel_marker_from_zhai JBrowse track. This is marker annotation, "
            "not per-sample genotype data."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR), help="Output directory")
    parser.add_argument("--target", action="append", choices=sorted(TARGETS),
                        help="Target region to extract; repeatable. Defaults to all targets.")
    parser.add_argument("--padding-bp", type=int, default=DEFAULT_PADDING_BP,
                        help="Flanking sequence to include on each side of the target interval")
    return parser


def parse_region(region: str) -> Tuple[str, int, int]:
    chrom, rest = region.split(":", 1)
    start_text, end_text = rest.split("-", 1)
    return chrom, int(start_text), int(end_text)


def fetch_json(url: str) -> object:
    req = request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with request.urlopen(req, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def decode_feature(raw: List[object], classes: List[Dict[str, object]]) -> Dict[str, object]:
    class_id = int(raw[0])
    attrs = classes[class_id].get("attributes", []) if class_id < len(classes) else []
    feature: Dict[str, object] = {"class_id": class_id}
    for attr, value in zip(attrs, raw[1:]):
        feature[str(attr)] = value
    return feature


def flatten_features(features: Iterable[List[object]]) -> Iterable[List[object]]:
    for feature in features:
        yield feature
        for value in feature:
            if isinstance(value, list):
                for child in value:
                    if isinstance(child, list) and child and isinstance(child[0], int):
                        yield child
            elif isinstance(value, dict):
                for sublist in value.values():
                    if isinstance(sublist, list):
                        for child in sublist:
                            if isinstance(child, list) and child and isinstance(child[0], int):
                                yield child


def target_rows(target_id: str, output_dir: Path, padding_bp: int) -> List[Dict[str, object]]:
    chrom, target_start, target_end = parse_region(TARGETS[target_id])
    start = max(1, target_start - padding_bp)
    end = target_end + padding_bp
    track_url = f"{TRACK_BASE}/{chrom}/trackData.json"
    track_data = fetch_json(track_url)
    intervals = track_data.get("intervals", {})
    classes = intervals.get("classes", [])
    chunk_template = intervals.get("urlTemplate", "lf-{Chunk}.json")
    rows: List[Dict[str, object]] = []

    raw_dir = output_dir / "raw_track_chunks" / chrom
    raw_dir.mkdir(parents=True, exist_ok=True)
    track_copy = raw_dir / "trackData.json"
    if not track_copy.exists():
        track_copy.write_text(json.dumps(track_data, ensure_ascii=False, indent=2), encoding="utf-8")

    for lazy_feature in intervals.get("nclist", []):
        if len(lazy_feature) < 4 or int(lazy_feature[0]) != int(intervals.get("lazyClass", 128)):
            continue
        chunk_start = int(lazy_feature[1])
        chunk_end = int(lazy_feature[2])
        chunk_id = str(lazy_feature[3])
        if chunk_end < start or chunk_start > end:
            continue

        chunk_name = chunk_template.replace("{Chunk}", chunk_id)
        chunk_url = f"{TRACK_BASE}/{chrom}/{chunk_name}"
        chunk_data = fetch_json(chunk_url)
        chunk_copy = raw_dir / chunk_name
        if not chunk_copy.exists():
            chunk_copy.write_text(json.dumps(chunk_data, ensure_ascii=False, indent=2), encoding="utf-8")

        for raw_feature in flatten_features(chunk_data):
            if not raw_feature or int(raw_feature[0]) == int(intervals.get("lazyClass", 128)):
                continue
            feature = decode_feature(raw_feature, classes)
            try:
                marker_start = int(feature.get("Start", -1))
                marker_end = int(feature.get("End", -1))
            except (TypeError, ValueError):
                continue
            if marker_end < start or marker_start > end:
                continue
            rows.append({
                "target_id": target_id,
                "region": TARGETS[target_id],
                "query_region": f"{chrom}:{start}-{end}",
                "padding_bp": padding_bp,
                "chrom": chrom,
                "start": marker_start,
                "end": marker_end,
                "name": feature.get("Name", ""),
                "id": feature.get("Id", ""),
                "type": feature.get("Type", ""),
                "source": feature.get("Source", ""),
                "polymorphism_ratio": feature.get("Polymorphism_ratio", ""),
                "indel_size_bp": feature.get("Indel_size_bp", ""),
                "amplicon_size_cs_bp": feature.get("Amplicon_size_in_chinese_spring_bp", ""),
                "amplicon_size_e6015_3s_bp": feature.get("Amplicon_size_in_e6015_3s_bp", ""),
                "forward_primer": feature.get("Forward.primer", ""),
                "reverse_primer": feature.get("Reverse.primer", ""),
                "seq_id": feature.get("Seq_id", ""),
                "track_url": track_url,
                "chunk_url": chunk_url,
            })

    # The track encodes transcript and exon subfeatures for the same marker.
    # Keep the first row for each coordinate/name/type tuple.
    deduped: List[Dict[str, object]] = []
    seen = set()
    for row in rows:
        key = (row["target_id"], row["start"], row["end"], row["name"], row["type"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)
    return deduped


def write_tsv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_arg_parser().parse_args(argv)
    output_dir = Path(args.output_dir)
    targets = args.target or sorted(TARGETS)
    output_dir.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "target_id", "region", "query_region", "padding_bp", "chrom", "start", "end",
        "name", "id", "type", "source", "polymorphism_ratio", "indel_size_bp",
        "amplicon_size_cs_bp", "amplicon_size_e6015_3s_bp", "forward_primer",
        "reverse_primer", "seq_id", "track_url", "chunk_url",
    ]
    all_rows: List[Dict[str, object]] = []
    summary_rows: List[Dict[str, object]] = []
    for target_id in targets:
        rows = target_rows(target_id, output_dir, args.padding_bp)
        all_rows.extend(rows)
        chrom, target_start, target_end = parse_region(TARGETS[target_id])
        summary_rows.append({
            "target_id": target_id,
            "region": TARGETS[target_id],
            "query_region": f"{chrom}:{max(1, target_start - args.padding_bp)}-{target_end + args.padding_bp}",
            "padding_bp": args.padding_bp,
            "marker_rows": len(rows),
        })
        print(f"[INFO] {target_id}: marker_rows={len(rows)}")

    write_tsv(output_dir / "zhai_indel_marker_annotations.tsv", all_rows, fieldnames)
    write_tsv(output_dir / "zhai_indel_marker_status.tsv", summary_rows,
              ["target_id", "region", "query_region", "padding_bp", "marker_rows"])
    (output_dir / "zhai_indel_marker_annotations.json").write_text(
        json.dumps(all_rows, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"[INFO] Wrote {output_dir / 'zhai_indel_marker_annotations.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
