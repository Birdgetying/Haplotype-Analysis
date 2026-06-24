#!/usr/bin/env python3
"""
Star-gene positive-control validation framework for HaplotypeScorer.

The default mode is intentionally conservative: inspect the manifest, local files,
format hints, sample overlap, and coordinate readiness without downloading data or
running large analyses. Real paper-data analysis only runs when --run-analysis is
explicitly requested and the target has resolved coordinates plus local inputs.
"""

import argparse
import csv
import gzip
import io
import json
import os
import sys
import traceback
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

if sys.platform == 'win32':
    try:
        if not isinstance(sys.stdout, io.TextIOWrapper) or sys.stdout.encoding != 'utf-8':
            sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        if not isinstance(sys.stderr, io.TextIOWrapper) or sys.stderr.encoding != 'utf-8':
            sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    except (ValueError, AttributeError):
        pass

try:
    import pandas as pd
except Exception:  # pragma: no cover - pandas is a project dependency, fallback keeps check-only importable
    pd = None


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_MANIFEST = PROJECT_ROOT / 'star_gene_manifest.json'
DEFAULT_DATABASE_ROOT = PROJECT_ROOT / 'star_gene_database'
DEFAULT_RESULTS_ROOT = PROJECT_ROOT / 'star_gene_results'

SUMMARY_COLUMNS = [
    'paper', 'doi', 'species', 'assembly', 'gene_or_locus', 'target_id',
    'target_type', 'trait', 'n_samples', 'n_haplotypes', 'n_variants',
    'association_pvalue', 'score_r_squared', 'score_regression_pvalue',
    'slope', 'pve', 'confidence_level', 'low_confidence',
    'circularity_warning', 'top_scored_haplotype',
    'top_scored_haplotype_score', 'top_haplotype_sample_count',
    'directional_top_haplotype', 'directional_top_haplotype_score',
    'directional_top_haplotype_mean', 'directional_top_haplotype_sample_count',
    'top_core_group', 'top_core_group_score', 'top_core_group_mean',
    'top_core_group_sample_count', 'directional_top_core_group',
    'directional_top_core_group_score', 'directional_top_core_group_mean',
    'directional_top_core_group_sample_count',
    'top_functional_group', 'top_functional_group_score',
    'top_functional_group_mean', 'top_functional_group_sample_count',
    'directional_top_functional_group',
    'directional_top_functional_group_score',
    'directional_top_functional_group_mean',
    'directional_top_functional_group_sample_count', 'expected_direction',
    'direction_consistency', 'data_status',
    'database_path', 'result_path', 'integrated_html', 'haplotype_score_html',
    'analysis_report_html', 'haplotype_score_json_path', 'score_mode', 'notes',
]

DB_REQUIRED_FILES = [
    'gene_info.json',
    'haplotype_data.csv',
    'haplotype_samples.csv',
    'phenotype_data.csv',
]

VCF_INDEX_SUFFIXES = ['.tbi', '.csi']


def _to_path(path_value: Optional[str], base_dir: Path = PROJECT_ROOT) -> Optional[Path]:
    if not path_value:
        return None
    path = Path(path_value)
    if path.is_absolute():
        return path
    return base_dir / path


def _json_dump(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def _safe_float(value: Any) -> Optional[float]:
    if value is None or value == '':
        return None
    try:
        if pd is not None and pd.isna(value):
            return None
    except Exception:
        pass
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _safe_int(value: Any) -> Optional[int]:
    if value is None or value == '':
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return None


def _first_present(mapping: Dict[str, Any], keys: Sequence[str]) -> Any:
    for key in keys:
        if key in mapping and mapping[key] is not None:
            return mapping[key]
    return None


def _open_text_maybe_gzip(path: Path):
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return open(path, 'r', encoding='utf-8', errors='replace')


def read_database_gene_info(database_path: Path) -> Tuple[Dict[str, Any], Optional[str]]:
    gene_info_path = database_path / 'gene_info.json'
    if not gene_info_path.exists():
        return {}, 'gene_info_missing'
    try:
        with open(gene_info_path, 'r', encoding='utf-8') as f:
            return json.load(f), None
    except Exception as e:
        return {}, f'gene_info_read_error:{e}'


def read_database_phenotype_columns(database_path: Path) -> Tuple[List[str], Optional[str]]:
    phenotype_path = database_path / 'phenotype_data.csv'
    if not phenotype_path.exists():
        return [], 'phenotype_data_missing'
    try:
        with open(phenotype_path, 'r', encoding='utf-8', newline='') as f:
            reader = csv.reader(f)
            header = next(reader, [])
    except Exception as e:
        return [], f'phenotype_data_read_error:{e}'

    reserved = {'SampleID', 'Hap_Name', 'Haplotype_Seq'}
    return [col for col in header if col not in reserved], None


def read_vcf_samples(vcf_path: Path, max_header_lines: int = 200000) -> Tuple[Set[str], Optional[str]]:
    """Read VCF sample IDs from the #CHROM header without scanning variants."""
    try:
        with _open_text_maybe_gzip(vcf_path) as f:
            for i, line in enumerate(f):
                if i > max_header_lines:
                    return set(), f'VCF header exceeded {max_header_lines} lines before #CHROM'
                if line.startswith('#CHROM'):
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) <= 9:
                        return set(), 'VCF header contains no sample columns'
                    return set(parts[9:]), None
        return set(), 'VCF #CHROM header not found'
    except Exception as e:
        return set(), f'failed_to_read_vcf_header: {e}'


def read_phenotype_samples(path: Path) -> Tuple[Set[str], Optional[str], List[str]]:
    """Read phenotype sample IDs and columns using a small set of separator guesses."""
    if pd is None:
        return set(), 'pandas_unavailable', []

    separators = ['\t', ',', r'\s+']
    errors: List[str] = []
    for sep in separators:
        try:
            df = pd.read_csv(path, sep=sep, engine='python' if sep == r'\s+' else 'c')
            if len(df.columns) < 2:
                errors.append(f'{sep}: only_one_column')
                continue
            columns = [str(c) for c in df.columns]
            sample_col = None
            for col in columns:
                if col.lower() in ['sampleid', 'sample_id', 'id', 'iid', 'fid', 'sample', 'vcfid']:
                    sample_col = col
                    break
            if sample_col is None:
                sample_col = columns[0]
            samples = set(df[sample_col].dropna().astype(str).str.strip())
            return samples, None, columns
        except Exception as e:
            errors.append(f'{sep}: {e}')
    return set(), '; '.join(errors), []


def has_vcf_index(path: Path) -> bool:
    return any(Path(str(path) + suffix).exists() for suffix in VCF_INDEX_SUFFIXES)


def classify_path_format(entry: Dict[str, Any], path: Path) -> str:
    fmt = str(entry.get('format') or '').lower()
    name = path.name.lower()
    if fmt:
        return fmt
    if name.endswith(('.vcf', '.vcf.gz')):
        return 'vcf'
    if name.endswith(('.tsv', '.csv', '.txt')):
        return 'table'
    if name.endswith(('.gff', '.gff3', '.gtf')):
        return 'annotation'
    return 'unknown'


def normalize_direction(direction: Optional[str]) -> str:
    if not direction:
        return 'unknown'
    value = str(direction).strip().lower()
    aliases = {
        'increase': 'increases_trait',
        'increases': 'increases_trait',
        'increases_trait': 'increases_trait',
        'higher': 'increases_trait',
        'positive': 'increases_trait',
        'decrease': 'decreases_trait',
        'decreases': 'decreases_trait',
        'decreases_trait': 'decreases_trait',
        'lower': 'decreases_trait',
        'negative': 'decreases_trait',
        'unknown': 'unknown',
        'na': 'unknown',
    }
    return aliases.get(value, value)


def pick_trait_columns(target: Dict[str, Any], paper_paths: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    traits = target.get('traits') or []
    if traits:
        return traits
    return [{
        'name': 'phenotype',
        'aliases': ['phenotype'],
        'local_columns': [],
        'expected_direction': target.get('expected_direction', 'unknown'),
    }]


def merge_path_entries(paper: Dict[str, Any], target: Dict[str, Any]) -> List[Dict[str, Any]]:
    entries: List[Dict[str, Any]] = []
    entries.extend(paper.get('local_expected_paths') or [])
    entries.extend(target.get('local_expected_paths') or [])
    return entries


def find_path_entry(entries: List[Dict[str, Any]], keys: Sequence[str] = (), formats: Sequence[str] = ()) -> Optional[Dict[str, Any]]:
    keys_l = {k.lower() for k in keys}
    formats_l = {f.lower() for f in formats}
    for entry in entries:
        entry_key = str(entry.get('key') or '').lower()
        entry_fmt = str(entry.get('format') or '').lower()
        if keys_l and entry_key in keys_l:
            return entry
        if formats_l and entry_fmt in formats_l:
            return entry
    return None


def load_manifest(path: Path) -> Dict[str, Any]:
    with open(path, 'r', encoding='utf-8') as f:
        manifest = json.load(f)
    if 'papers' not in manifest or not isinstance(manifest['papers'], list):
        raise ValueError('manifest must contain a list field named "papers"')
    return manifest


class StarGeneValidator:
    def __init__(self, manifest: Dict[str, Any], database_root: Path, results_root: Path,
                 check_only: bool = True, allow_large_download: bool = False,
                 accept_license: bool = False, score_mode: str = 'default'):
        self.manifest = manifest
        self.database_root = database_root
        self.results_root = results_root
        self.check_only = check_only
        self.allow_large_download = allow_large_download
        self.accept_license = accept_license
        self.score_mode = score_mode
        self.results_root.mkdir(parents=True, exist_ok=True)

    def iter_targets(self, papers: Optional[Sequence[str]] = None,
                     targets: Optional[Sequence[str]] = None) -> Iterable[Tuple[Dict[str, Any], Dict[str, Any]]]:
        paper_filter = {p.lower() for p in papers or []}
        target_filter = {t.lower() for t in targets or []}
        for paper in self.manifest.get('papers', []):
            paper_ids = {
                str(paper.get('paper_id', '')).lower(),
                str(paper.get('short_name', '')).lower(),
                str(paper.get('doi', '')).lower(),
            }
            if paper_filter and paper_filter.isdisjoint(paper_ids):
                continue
            for target in paper.get('targets', []) or []:
                aliases = target.get('aliases') or []
                target_ids = {
                    str(target.get('target_id', '')).lower(),
                    str(target.get('gene_or_locus', '')).lower(),
                    *{str(a).lower() for a in aliases},
                }
                if target_filter and target_filter.isdisjoint(target_ids):
                    continue
                yield paper, target

    def check_target(self, paper: Dict[str, Any], target: Dict[str, Any]) -> Dict[str, Any]:
        paper_id = paper.get('paper_id') or paper.get('short_name') or 'paper'
        target_id = target.get('target_id') or target.get('gene_or_locus') or 'target'
        entries = merge_path_entries(paper, target)
        database_path = self.database_root / str(paper_id) / str(target_id)
        result_path = self.results_root / str(paper_id) / str(target_id)
        if self.score_mode != 'default':
            result_path = result_path.with_name(f"{result_path.name}__{self.score_mode}")

        notes: List[str] = []
        missing_required: List[str] = []
        existing_paths: List[str] = []
        vcf_sample_sets: List[Set[str]] = []
        phenotype_samples: Set[str] = set()
        phenotype_columns: List[str] = []
        supported_vcf_path: Optional[Path] = None
        unsupported_vcf_reasons: List[str] = []
        unresolved_coordinates = bool(target.get('requires_coordinate_resolution')) or not target.get('coordinates')

        for entry in entries:
            path = _to_path(entry.get('path'))
            label = entry.get('key') or entry.get('path') or 'local_path'
            required = bool(entry.get('required'))
            fmt = classify_path_format(entry, path) if path else str(entry.get('format') or 'unknown')
            if not path or not path.exists():
                if required:
                    missing_required.append(str(label))
                notes.append(f'missing:{label}={path}')
                continue

            existing_paths.append(str(path))
            size_mb = path.stat().st_size / (1024 * 1024) if path.is_file() else 0
            notes.append(f'found:{label}({fmt},{size_mb:.2f}MB)')

            if fmt == 'vcf' or str(path).lower().endswith(('.vcf', '.vcf.gz')):
                vcf_index_ok = True
                if path.suffix == '.gz' or str(path).endswith('.vcf.gz'):
                    vcf_index_ok = has_vcf_index(path)
                    if vcf_index_ok:
                        notes.append(f'vcf_index_ok:{label}')
                    else:
                        notes.append(f'vcf_index_missing:{label}')
                        unsupported_vcf_reasons.append(f'{label}:missing_index')
                if vcf_index_ok:
                    supported_vcf_path = path
                samples, err = read_vcf_samples(path)
                if err:
                    notes.append(f'vcf_header_warning:{label}:{err}')
                elif samples:
                    vcf_sample_sets.append(samples)
                    notes.append(f'vcf_samples:{label}={len(samples)}')
            elif fmt in ['phenotype_table', 'phenotype', 'table'] or 'phenotype' in str(label).lower():
                samples, err, cols = read_phenotype_samples(path)
                phenotype_columns = cols
                if err:
                    notes.append(f'phenotype_warning:{label}:{err}')
                elif samples:
                    phenotype_samples = samples
                    notes.append(f'phenotype_samples:{label}={len(samples)}')
                    notes.append(f'phenotype_columns:{label}={",".join(cols[:20])}')

        database_complete = False
        database_source_ok = True
        required_database_source = target.get('required_database_source')
        if database_path.exists():
            missing_db = [name for name in DB_REQUIRED_FILES if not (database_path / name).exists()]
            if missing_db:
                notes.append(f'database_incomplete:missing={",".join(missing_db)}')
            else:
                database_complete = True
                notes.append('database_complete')
        else:
            notes.append('database_missing')

        database_gene_info: Dict[str, Any] = {}
        if database_complete:
            database_gene_info, gene_info_err = read_database_gene_info(database_path)
            if gene_info_err:
                notes.append(f'database_gene_info_warning:{gene_info_err}')
            else:
                if required_database_source:
                    actual_source = str(database_gene_info.get('source') or '')
                    if actual_source != str(required_database_source):
                        database_source_ok = False
                        database_complete = False
                        notes.append(
                            f'database_source_mismatch:required={required_database_source},actual={actual_source or "missing"}'
                        )
                db_chrom = database_gene_info.get('chrom') or database_gene_info.get('chr')
                db_start = _safe_int(database_gene_info.get('start') or database_gene_info.get('gene_start'))
                db_end = _safe_int(database_gene_info.get('end') or database_gene_info.get('gene_end'))
                if database_source_ok and db_chrom and db_start is not None and db_end is not None:
                    unresolved_coordinates = False
                    notes.append(f'database_coordinates:{target_id}={db_chrom}:{db_start}-{db_end}')
                elif database_source_ok:
                    notes.append('database_coordinates_missing')

            if database_source_ok:
                db_pheno_cols, db_pheno_err = read_database_phenotype_columns(database_path)
                if db_pheno_err:
                    notes.append(f'database_phenotype_warning:{db_pheno_err}')
                elif db_pheno_cols:
                    phenotype_columns = db_pheno_cols
                    notes.append(f'database_phenotype_columns:{",".join(db_pheno_cols[:20])}')

        if unresolved_coordinates:
            notes.append('coordinates_unresolved')

        if phenotype_samples and vcf_sample_sets:
            # Report overlap against the union because target projects may provide one VCF per locus.
            vcf_union = set().union(*vcf_sample_sets)
            overlap = phenotype_samples.intersection(vcf_union)
            notes.append(f'sample_overlap={len(overlap)}')
            if not overlap:
                notes.append('sample_overlap_empty')

        has_supported_genotype_source = database_complete or supported_vcf_path is not None
        if unsupported_vcf_reasons:
            notes.append('unsupported_vcf:' + ','.join(unsupported_vcf_reasons))

        if database_complete:
            missing_required = []

        if not database_source_ok:
            data_status = 'unsupported_input_format_for_analysis'
        elif missing_required:
            data_status = 'missing_required_files'
        elif unresolved_coordinates:
            data_status = 'requires_coordinate_resolution'
        elif not has_supported_genotype_source:
            data_status = 'unsupported_input_format_for_analysis'
        elif self.check_only:
            data_status = 'ready_for_analysis_check_only'
        else:
            data_status = 'ready_for_analysis'

        if self.allow_large_download:
            if not self.accept_license:
                notes.append('allow_large_download_ignored_without_accept_license')
            else:
                notes.append('download_not_implemented_in_minimal_framework')

        rows = []
        for trait in pick_trait_columns(target, entries):
            row = self.empty_summary_row(paper, target, trait, database_path, result_path)
            row.update({
                'data_status': data_status,
                'notes': '; '.join(notes),
            })
            rows.append(row)
        return {
            'row': rows[0] if rows else self.empty_summary_row(paper, target, {}, database_path, result_path),
            'rows': rows,
            'status': data_status,
            'notes': notes,
            'entries': entries,
            'database_path': database_path,
            'result_path': result_path,
            'phenotype_columns': phenotype_columns,
        }

    def empty_summary_row(self, paper: Dict[str, Any], target: Dict[str, Any], trait: Dict[str, Any],
                          database_path: Path, result_path: Path) -> Dict[str, Any]:
        paper_id = paper.get('paper_id') or paper.get('short_name') or ''
        target_id = target.get('target_id') or target.get('gene_or_locus') or ''
        expected_direction = normalize_direction(
            trait.get('expected_direction') or target.get('expected_direction') or 'unknown'
        )
        return {
            'paper': paper_id,
            'doi': paper.get('doi'),
            'species': paper.get('species'),
            'assembly': paper.get('assembly'),
            'gene_or_locus': target.get('gene_or_locus'),
            'target_id': target_id,
            'target_type': target.get('target_type'),
            'trait': trait.get('name') or 'phenotype',
            'n_samples': None,
            'n_haplotypes': None,
            'n_variants': None,
            'association_pvalue': None,
            'score_r_squared': None,
            'score_regression_pvalue': None,
            'slope': None,
            'pve': None,
            'confidence_level': None,
            'low_confidence': None,
            'circularity_warning': None,
            'top_scored_haplotype': None,
            'top_scored_haplotype_score': None,
            'top_haplotype_sample_count': None,
            'expected_direction': expected_direction,
            'direction_consistency': 'NA' if expected_direction == 'unknown' else None,
            'data_status': None,
            'database_path': str(database_path),
            'result_path': str(result_path),
            'integrated_html': None,
            'haplotype_score_html': None,
            'analysis_report_html': None,
            'haplotype_score_json_path': None,
            'score_mode': self.score_mode,
            'notes': '',
        }

    def resolve_target_coordinates(self, target: Dict[str, Any],
                                   database_path: Optional[Path] = None) -> Tuple[Optional[str], Optional[int], Optional[int]]:
        coords = target.get('coordinates') or {}
        chrom = coords.get('chrom') or coords.get('chr')
        start = _safe_int(coords.get('start'))
        end = _safe_int(coords.get('end'))
        if chrom and start is not None and end is not None:
            return str(chrom), start, end

        if database_path:
            gene_info, err = read_database_gene_info(database_path)
            if not err:
                chrom = gene_info.get('chrom') or gene_info.get('chr')
                start = _safe_int(gene_info.get('start') or gene_info.get('gene_start'))
                end = _safe_int(gene_info.get('end') or gene_info.get('gene_end'))
                if chrom and start is not None and end is not None:
                    return str(chrom), start, end

        return None, None, None

    def run_target(self, paper: Dict[str, Any], target: Dict[str, Any], check: Dict[str, Any]) -> List[Dict[str, Any]]:
        if self.check_only:
            return check.get('rows') or [check['row']]
        if self.score_mode not in ('default', 'robust_discovery'):
            rows = []
            for base_row in check.get('rows') or [check['row']]:
                row = dict(base_row)
                row['data_status'] = 'skipped_non_default_score_mode_not_wired'
                row['notes'] = (row.get('notes') or '') + '; score mode is not wired for analyzer execution'
                rows.append(row)
            return rows
        if check['status'] not in ['ready_for_analysis', 'ready_for_analysis_check_only']:
            rows = []
            for base_row in check.get('rows') or [check['row']]:
                row = dict(base_row)
                row['data_status'] = f'skipped_{check["status"]}'
                rows.append(row)
            return rows

        chrom, start, end = self.resolve_target_coordinates(target, check['database_path'])
        if not chrom or start is None or end is None:
            row = dict(check['row'])
            row['data_status'] = 'skipped_missing_coordinates'
            row['notes'] = (row.get('notes') or '') + '; coordinates missing chrom/start/end'
            return [row]

        entries = check['entries']
        vcf_entry = find_path_entry(entries, keys=['target_vcf', 'vcf'], formats=['vcf'])
        pheno_entry = find_path_entry(entries, keys=['phenotype_table', 'phenotype'], formats=['phenotype_table', 'phenotype'])
        gtf_entry = find_path_entry(entries, keys=['gtf_file', 'gff_file', 'annotation'], formats=['gtf', 'gff', 'gff3', 'annotation'])

        vcf_path = _to_path(vcf_entry.get('path')) if vcf_entry else None
        pheno_path = _to_path(pheno_entry.get('path')) if pheno_entry else None
        gtf_path = _to_path(gtf_entry.get('path')) if gtf_entry else None
        if (pheno_path is None or not pheno_path.exists()) and check['database_path'].exists():
            database_pheno_path = check['database_path'] / 'phenotype_data.csv'
            if database_pheno_path.exists():
                pheno_path = database_pheno_path
        database_root_for_paper = self.database_root / str(paper.get('paper_id') or paper.get('short_name') or 'paper')
        result_path = check['result_path']
        result_path.mkdir(parents=True, exist_ok=True)

        try:
            from haplotype_phenotype_analysis import HaplotypePhenotypeAnalyzer
        except Exception as e:
            row = dict(check['row'])
            row['data_status'] = 'skipped_import_error'
            row['notes'] = (row.get('notes') or '') + f'; import_error:{e}'
            return [row]

        gene_id = target.get('analyzer_gene_id') or target.get('target_id') or target.get('gene_or_locus')
        trait_rows: List[Dict[str, Any]] = []
        traits = pick_trait_columns(target, entries)
        pheno_columns = check.get('phenotype_columns') or []
        selected_traits: List[Tuple[Dict[str, Any], str]] = []

        for trait in traits:
            selected_col = self.resolve_trait_column(trait, pheno_columns)
            if selected_col:
                selected_traits.append((trait, selected_col))
            else:
                row = self.empty_summary_row(paper, target, trait, check['database_path'], result_path)
                row['data_status'] = 'skipped_trait_column_missing'
                row['notes'] = (check['row'].get('notes') or '') + '; trait column not found in phenotype table'
                trait_rows.append(row)

        if not selected_traits:
            return trait_rows

        try:
            analyzer = HaplotypePhenotypeAnalyzer(
                vcf_file=str(vcf_path) if vcf_path and vcf_path.exists() else None,
                phenotype_file=str(pheno_path) if pheno_path and pheno_path.exists() else None,
                output_dir=str(result_path),
                gtf_file=str(gtf_path) if gtf_path and gtf_path.exists() else None,
                score_mode=self.score_mode,
            )
            phenotype_cols = []
            for _, selected_col in selected_traits:
                if selected_col not in phenotype_cols:
                    phenotype_cols.append(selected_col)
            result = analyzer.analyze_gene(
                chrom=str(chrom),
                start=start,
                end=end,
                gene_id=str(gene_id),
                phenotype_cols=phenotype_cols,
                cluster_haplotypes=False,
                database_dir=str(database_root_for_paper),
            )
        except Exception as e:
            for trait, _ in selected_traits:
                row = self.empty_summary_row(paper, target, trait, check['database_path'], result_path)
                row['data_status'] = 'analysis_error'
                row['notes'] = (check['row'].get('notes') or '') + f'; analysis_error:{e}; traceback:{traceback.format_exc(limit=3)}'
                trait_rows.append(row)
            return trait_rows

        for trait, selected_col in selected_traits:
            row = self.empty_summary_row(paper, target, trait, check['database_path'], result_path)
            summary = self.summarize_analysis_result(
                paper, target, trait, result, result_path, selected_phenotype_col=selected_col
            )
            row.update(summary)
            row['data_status'] = 'analyzed' if summary.get('trait') else 'analysis_trait_missing'
            row['notes'] = check['row'].get('notes') or ''
            if row['data_status'] == 'analysis_trait_missing':
                row['notes'] += f'; analyzer result missing selected trait column: {selected_col}'
            trait_rows.append(row)

        self.write_literature_audit(paper, target, check['database_path'], result_path, trait_rows)
        return trait_rows

    def write_literature_audit(self, paper: Dict[str, Any], target: Dict[str, Any],
                               database_path: Path, result_path: Path,
                               rows: List[Dict[str, Any]]) -> None:
        if not (target.get('literature_variants') or target.get('literature_haplotypes')):
            return
        try:
            from star_gene_literature_audit import run_literature_audit

            records = run_literature_audit(paper, target, database_path, result_path, summary_rows=rows)
            note = f'literature_audit_records={len(records)}'
            if records:
                note += f'; literature_audit={result_path / "literature_variant_audit.csv"}'
        except Exception as e:
            note = f'literature_audit_error:{e}'
        for row in rows:
            existing = row.get('notes') or ''
            row['notes'] = f'{existing}; {note}' if existing else note

    def resolve_trait_column(self, trait: Dict[str, Any], phenotype_columns: Sequence[str]) -> Optional[str]:
        requested = []
        requested.extend(trait.get('local_columns') or [])
        requested.append(trait.get('name'))
        requested.extend(trait.get('aliases') or [])
        requested = [str(x) for x in requested if x]
        if not phenotype_columns:
            return requested[0] if requested else None
        exact = {c: c for c in phenotype_columns}
        lower = {c.lower(): c for c in phenotype_columns}
        for name in requested:
            if name in exact:
                return exact[name]
            if name.lower() in lower:
                return lower[name.lower()]
        return None

    def summarize_analysis_result(self, paper: Dict[str, Any], target: Dict[str, Any], trait: Dict[str, Any],
                                  result: Dict[str, Any], result_path: Path,
                                  selected_phenotype_col: Optional[str] = None) -> Dict[str, Any]:
        target_id = target.get('target_id') or target.get('gene_or_locus')
        trait_name = trait.get('name') or 'phenotype'
        phenotype_results = result.get('phenotype_results') or {}
        pheno_key = None
        candidates = []
        if selected_phenotype_col:
            candidates.append(selected_phenotype_col)
        candidates.extend([trait_name, *(trait.get('local_columns') or []), *(trait.get('aliases') or [])])
        for candidate in [str(c) for c in candidates if c]:
            if candidate in phenotype_results:
                pheno_key = candidate
                break
        if pheno_key is None:
            return {'trait': None}

        pheno_result = phenotype_results.get(pheno_key, {})
        assoc = pheno_result.get('association') or {}
        pve_result = pheno_result.get('pve') or {}
        score_results = pheno_result.get('haplotype_score') or {}
        if not score_results:
            all_scores = result.get('haplotype_scores') or {}
            if pheno_key in all_scores:
                score_results = all_scores[pheno_key]

        per_sample = score_results.get('per_sample') or []
        per_haplotype = score_results.get('per_haplotype') or {}
        top_hap, top_score = self.pick_top_haplotype(per_haplotype)
        top_count = None
        if top_hap is not None:
            top_count = sum(1 for sample in per_sample if str(sample.get('haplotype')) == str(top_hap))

        expected_direction = normalize_direction(trait.get('expected_direction') or target.get('expected_direction'))
        directional_top_hap, directional_top_score, directional_top_mean, directional_top_count = (
            self.pick_directional_top_haplotype(score_results, expected_direction)
        )
        top_core_group, top_core_score, top_core_mean, top_core_count = (
            self.pick_top_core_group(score_results)
        )
        directional_core_group, directional_core_score, directional_core_mean, directional_core_count = (
            self.pick_directional_top_core_group(score_results, expected_direction)
        )
        top_functional_group, top_functional_score, top_functional_mean, top_functional_count = (
            self.pick_top_functional_group(score_results)
        )
        directional_functional_group, directional_functional_score, directional_functional_mean, directional_functional_count = (
            self.pick_directional_top_functional_group(score_results, expected_direction)
        )
        direction_consistency = self.direction_consistency(expected_direction, score_results, top_hap)

        integrated_html = result_path / f'{target_id}.html'
        if not integrated_html.exists():
            alt = result_path / 'integrated_analysis.html'
            integrated_html = alt if alt.exists() else integrated_html
        hap_score_html = result_path / 'haplotype_score.html'
        analysis_report = result_path / 'analysis_report.html'

        return {
            'trait': pheno_key or trait_name,
            'n_samples': len(per_sample) if per_sample else None,
            'n_haplotypes': len(per_haplotype) if per_haplotype else result.get('n_haplotypes'),
            'n_variants': result.get('n_variants'),
            'association_pvalue': _first_present(assoc, ['p_value_corrected', 'p_value', 'anova_pvalue']),
            'score_r_squared': score_results.get('r_squared'),
            'score_regression_pvalue': score_results.get('regression_pvalue'),
            'slope': score_results.get('slope'),
            'pve': _first_present(pve_result, ['eta_squared', 'pve_percent', 'partial_r_squared']),
            'confidence_level': score_results.get('confidence_level'),
            'low_confidence': score_results.get('low_confidence'),
            'circularity_warning': score_results.get('circularity_warning'),
            'top_scored_haplotype': top_hap,
            'top_scored_haplotype_score': top_score,
            'top_haplotype_sample_count': top_count,
            'directional_top_haplotype': directional_top_hap,
            'directional_top_haplotype_score': directional_top_score,
            'directional_top_haplotype_mean': directional_top_mean,
            'directional_top_haplotype_sample_count': directional_top_count,
            'top_core_group': top_core_group,
            'top_core_group_score': top_core_score,
            'top_core_group_mean': top_core_mean,
            'top_core_group_sample_count': top_core_count,
            'directional_top_core_group': directional_core_group,
            'directional_top_core_group_score': directional_core_score,
            'directional_top_core_group_mean': directional_core_mean,
            'directional_top_core_group_sample_count': directional_core_count,
            'top_functional_group': top_functional_group,
            'top_functional_group_score': top_functional_score,
            'top_functional_group_mean': top_functional_mean,
            'top_functional_group_sample_count': top_functional_count,
            'directional_top_functional_group': directional_functional_group,
            'directional_top_functional_group_score': directional_functional_score,
            'directional_top_functional_group_mean': directional_functional_mean,
            'directional_top_functional_group_sample_count': directional_functional_count,
            'expected_direction': expected_direction,
            'direction_consistency': direction_consistency,
            'integrated_html': str(integrated_html) if integrated_html.exists() else None,
            'haplotype_score_html': str(hap_score_html) if hap_score_html.exists() else None,
            'analysis_report_html': str(analysis_report) if analysis_report.exists() else None,
            'haplotype_score_json_path': result.get('haplotype_score_json_path'),
        }

    def pick_top_haplotype(self, per_haplotype: Dict[str, Dict[str, Any]]) -> Tuple[Optional[str], Optional[float]]:
        best_hap = None
        best_score = None
        for hap, values in per_haplotype.items():
            score = _safe_float((values or {}).get('total'))
            if score is None:
                continue
            if best_score is None or score > best_score:
                best_hap = hap
                best_score = score
        return best_hap, best_score

    @staticmethod
    def _stable_support_candidates(candidates: Sequence[Tuple[Any, ...]], count_index: int) -> List[Tuple[Any, ...]]:
        """Filter low-support outliers when stable alternatives exist."""
        candidates = list(candidates)
        counts = [int(item[count_index] or 0) for item in candidates if int(item[count_index] or 0) > 0]
        if not counts:
            return candidates
        max_count = max(counts)
        min_count = min(max_count, max(3, min(10, int(max_count * 0.10))))
        stable_candidates = [item for item in candidates if int(item[count_index] or 0) >= min_count]
        return stable_candidates or candidates

    @staticmethod
    def _directional_candidate_value(
        candidate: Tuple[Any, ...],
        baseline: float,
        mean_index: int,
        count_index: int,
        reliability_k: float = 20.0,
    ) -> float:
        mean = float(candidate[mean_index])
        count = max(int(candidate[count_index] or 0), 0)
        reliability = count / (count + reliability_k) if count > 0 else 0.0
        return baseline + (mean - baseline) * reliability

    def pick_directional_top_haplotype(self, score_results: Dict[str, Any],
                                       expected_direction: str) -> Tuple[Optional[str], Optional[float], Optional[float], Optional[int]]:
        """Pick a phenotype-direction-aware top haplotype without using literature allele labels."""
        expected_direction = normalize_direction(expected_direction)
        per_haplotype = score_results.get('per_haplotype') or {}
        per_sample = score_results.get('per_sample') or []
        if not per_haplotype:
            return None, None, None, None
        if expected_direction == 'unknown':
            hap, score = self.pick_top_haplotype(per_haplotype)
            return hap, score, None, None

        stats_by_hap: Dict[str, Dict[str, float]] = {}
        for sample in per_sample:
            hap = str(sample.get('haplotype'))
            pheno = _safe_float(sample.get('phenotype'))
            if not hap or pheno is None:
                continue
            entry = stats_by_hap.setdefault(hap, {'n': 0.0, 'sum': 0.0})
            entry['n'] += 1.0
            entry['sum'] += pheno

        candidates: List[Tuple[str, float, float, int]] = []
        for hap, values in per_haplotype.items():
            score = _safe_float((values or {}).get('total'))
            if score is None:
                continue
            entry = stats_by_hap.get(str(hap), {})
            n = int(entry.get('n', 0.0) or 0)
            if n:
                mean = float(entry['sum'] / n)
            else:
                mean = _safe_float((values or {}).get('mean_phenotype'))
            if mean is None:
                continue
            candidates.append((str(hap), score, mean, n))

        if not candidates:
            hap, score = self.pick_top_haplotype(per_haplotype)
            return hap, score, None, None
        candidates = self._stable_support_candidates(candidates, count_index=3)
        total_n = sum(max(item[3], 0) for item in candidates)
        if total_n > 0:
            baseline = sum(item[2] * max(item[3], 0) for item in candidates) / total_n
        else:
            baseline = sum(item[2] for item in candidates) / len(candidates)
        if expected_direction == 'decreases_trait':
            return min(
                candidates,
                key=lambda item: (
                    self._directional_candidate_value(item, baseline, mean_index=2, count_index=3),
                    -item[1],
                    -item[3],
                )
            )
        return max(
            candidates,
            key=lambda item: (
                self._directional_candidate_value(item, baseline, mean_index=2, count_index=3),
                item[1],
                item[3],
            )
        )

    def pick_top_core_group(self, score_results: Dict[str, Any]) -> Tuple[Optional[str], Optional[float], Optional[float], Optional[int]]:
        """Pick the raw top robust core group by reliability-adjusted score."""
        core = score_results.get('core_haplotype_groups') or {}
        top = core.get('top_group') or {}
        group = str(top.get('core_sequence') or '').strip()
        if group:
            return (
                group,
                _safe_float(top.get('rank_score') or top.get('mean_score')),
                _safe_float(top.get('mean_phenotype')),
                int(_safe_float(top.get('sample_count')) or 0),
            )
        groups = core.get('groups') or {}
        if not isinstance(groups, dict) or not groups:
            return None, None, None, None
        candidates = []
        for group_name, values in groups.items():
            values = values or {}
            score = _safe_float(values.get('rank_score') or values.get('mean_score'))
            if score is None:
                continue
            candidates.append((
                str(values.get('core_sequence') or group_name),
                score,
                _safe_float(values.get('mean_phenotype')),
                int(_safe_float(values.get('sample_count')) or 0),
            ))
        if not candidates:
            return None, None, None, None
        return max(candidates, key=lambda item: (item[1], item[3]))

    def pick_directional_top_core_group(self, score_results: Dict[str, Any],
                                        expected_direction: str) -> Tuple[Optional[str], Optional[float], Optional[float], Optional[int]]:
        """Pick a phenotype-direction-aware robust core group without literature labels."""
        return self._pick_directional_top_position_group(
            score_results,
            expected_direction,
            group_key='core_haplotype_groups',
            sequence_key='core_sequence',
        )

    def pick_top_functional_group(self, score_results: Dict[str, Any]) -> Tuple[Optional[str], Optional[float], Optional[float], Optional[int]]:
        """Pick the raw top functional sub-haplotype group."""
        return self._pick_top_position_group(
            score_results,
            group_key='functional_haplotype_groups',
            sequence_key='functional_sequence',
        )

    def pick_directional_top_functional_group(self, score_results: Dict[str, Any],
                                             expected_direction: str) -> Tuple[Optional[str], Optional[float], Optional[float], Optional[int]]:
        """Pick a direction-aware functional sub-haplotype group."""
        return self._pick_directional_top_position_group(
            score_results,
            expected_direction,
            group_key='functional_haplotype_groups',
            sequence_key='functional_sequence',
        )

    def _pick_top_position_group(self, score_results: Dict[str, Any],
                                 group_key: str,
                                 sequence_key: str) -> Tuple[Optional[str], Optional[float], Optional[float], Optional[int]]:
        groups_payload = score_results.get(group_key) or {}
        top = groups_payload.get('top_group') or {}
        group = str(top.get(sequence_key) or '').strip()
        if group:
            return (
                group,
                _safe_float(top.get('rank_score') or top.get('mean_score')),
                _safe_float(top.get('mean_phenotype')),
                int(_safe_float(top.get('sample_count')) or 0),
            )
        groups = groups_payload.get('groups') or {}
        if not isinstance(groups, dict) or not groups:
            return None, None, None, None
        candidates = []
        for group_name, values in groups.items():
            values = values or {}
            score = _safe_float(values.get('rank_score') or values.get('mean_score'))
            if score is None:
                continue
            candidates.append((
                str(values.get(sequence_key) or group_name),
                score,
                _safe_float(values.get('mean_phenotype')),
                int(_safe_float(values.get('sample_count')) or 0),
            ))
        if not candidates:
            return None, None, None, None
        return max(candidates, key=lambda item: (item[1], item[3]))

    def _pick_directional_top_position_group(self, score_results: Dict[str, Any],
                                             expected_direction: str,
                                             group_key: str,
                                             sequence_key: str) -> Tuple[Optional[str], Optional[float], Optional[float], Optional[int]]:
        """Pick a phenotype-direction-aware position group without literature labels."""
        expected_direction = normalize_direction(expected_direction)
        if expected_direction == 'unknown':
            return self._pick_top_position_group(score_results, group_key, sequence_key)
        groups_payload = score_results.get(group_key) or {}
        groups = groups_payload.get('groups') or {}
        if not isinstance(groups, dict) or not groups:
            return None, None, None, None
        candidates: List[Tuple[str, float, float, int]] = []
        for group_name, values in groups.items():
            values = values or {}
            mean = _safe_float(values.get('mean_phenotype'))
            if mean is None:
                continue
            score = _safe_float(values.get('rank_score') or values.get('mean_score')) or 0.0
            count = int(_safe_float(values.get('sample_count')) or 0)
            candidates.append((str(values.get(sequence_key) or group_name), score, mean, count))
        if not candidates:
            return self._pick_top_position_group(score_results, group_key, sequence_key)
        candidates = self._stable_support_candidates(candidates, count_index=3)
        total_n = sum(max(item[3], 0) for item in candidates)
        if total_n > 0:
            baseline = sum(item[2] * max(item[3], 0) for item in candidates) / total_n
        else:
            baseline = sum(item[2] for item in candidates) / len(candidates)
        if expected_direction == 'decreases_trait':
            return min(
                candidates,
                key=lambda item: (
                    self._directional_candidate_value(item, baseline, mean_index=2, count_index=3),
                    -item[1],
                    -item[3],
                )
            )
        return max(
            candidates,
            key=lambda item: (
                self._directional_candidate_value(item, baseline, mean_index=2, count_index=3),
                item[1],
                item[3],
            )
        )

    def direction_consistency(self, expected_direction: str, score_results: Dict[str, Any],
                              top_hap: Optional[str]) -> str:
        expected_direction = normalize_direction(expected_direction)
        if expected_direction == 'unknown':
            return 'NA'
        slope = _safe_float(score_results.get('slope'))
        per_sample = score_results.get('per_sample') or []
        if slope is None or top_hap is None or not per_sample:
            return 'unknown'

        all_phenos = [_safe_float(s.get('phenotype')) for s in per_sample]
        top_phenos = [_safe_float(s.get('phenotype')) for s in per_sample if str(s.get('haplotype')) == str(top_hap)]
        all_phenos = [x for x in all_phenos if x is not None]
        top_phenos = [x for x in top_phenos if x is not None]
        if not all_phenos or not top_phenos:
            return 'unknown'
        overall_mean = sum(all_phenos) / len(all_phenos)
        top_mean = sum(top_phenos) / len(top_phenos)

        if expected_direction == 'increases_trait':
            return 'consistent' if slope > 0 and top_mean > overall_mean else 'inconsistent'
        if expected_direction == 'decreases_trait':
            return 'consistent' if slope < 0 and top_mean < overall_mean else 'inconsistent'
        return 'unknown'

    def run(self, papers: Optional[Sequence[str]] = None,
            targets: Optional[Sequence[str]] = None) -> List[Dict[str, Any]]:
        rows: List[Dict[str, Any]] = []
        selected = list(self.iter_targets(papers=papers, targets=targets))
        if not selected:
            raise ValueError('No manifest targets matched the requested --paper/--target filters')

        print(f"[INFO] Mode: {'check-only' if self.check_only else 'run-analysis'}")
        print(f"[INFO] Targets selected: {len(selected)}")
        print(f"[INFO] Results root: {self.results_root}")
        print(f"[INFO] Database root: {self.database_root}")

        for paper, target in selected:
            paper_id = paper.get('paper_id') or paper.get('short_name')
            target_id = target.get('target_id') or target.get('gene_or_locus')
            print(f"\n[CHECK] {paper_id} / {target_id}")
            check = self.check_target(paper, target)
            print(f"  status: {check['status']}")
            for note in check['notes'][:12]:
                print(f"  - {note}")
            if len(check['notes']) > 12:
                print(f"  - ... {len(check['notes']) - 12} more notes")
            rows.extend(self.run_target(paper, target, check))

        self.write_summary(rows)
        return rows

    def write_summary(self, rows: List[Dict[str, Any]]) -> None:
        csv_path = self.results_root / 'validation_summary.csv'
        json_path = self.results_root / 'validation_summary.json'
        with open(csv_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=SUMMARY_COLUMNS, extrasaction='ignore')
            writer.writeheader()
            for row in rows:
                writer.writerow({col: row.get(col) for col in SUMMARY_COLUMNS})
        _json_dump(json_path, rows)
        print(f"\n[INFO] Summary CSV: {csv_path}")
        print(f"[INFO] Summary JSON: {json_path}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description='Validate HaplotypeScorer against paper star genes/loci (safe check-only by default).',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--manifest', default=str(DEFAULT_MANIFEST), help='star-gene manifest JSON path')
    parser.add_argument('--database-root', default=str(DEFAULT_DATABASE_ROOT), help='root for per-paper star-gene databases')
    parser.add_argument('--results-root', default=str(DEFAULT_RESULTS_ROOT), help='root for validation summaries/results')
    parser.add_argument('--paper', action='append', default=[], help='paper_id, short_name, or DOI to include; repeatable')
    parser.add_argument('--target', action='append', default=[], help='target_id, gene/locus, or alias to include; repeatable')
    parser.add_argument('--all', action='store_true', help='explicitly include all manifest targets')

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--check-only', dest='check_only', action='store_true', default=True,
                      help='inspect manifest/local files only; do not run analyzer')
    mode.add_argument('--run-analysis', dest='check_only', action='store_false',
                      help='run analyzer for targets with resolved coordinates and local inputs')

    parser.add_argument('--no-download', action='store_true', default=True,
                        help='do not download data (default; retained for explicitness)')
    parser.add_argument('--allow-large-download', action='store_true',
                        help='placeholder flag for future large downloads; no downloads are implemented in this minimal framework')
    parser.add_argument('--accept-license', action='store_true',
                        help='confirm external data terms before any future download path is enabled')
    parser.add_argument('--score-mode', choices=['default', 'no_finemap', 'functional_only', 'robust_discovery'],
                        default='default',
                        help='haplotype scoring mode; robust_discovery penalizes tiny/ambiguous haplotypes for discovery ranking')
    parser.add_argument('--print-downloads', action='store_true',
                        help='print lightweight data-source notes and copyable download commands, then exit')
    parser.add_argument('--include-large-downloads', action='store_true',
                        help='include large/manual sources in printed notes; direct commands are still only emitted for safe small files')
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    manifest_path = _to_path(args.manifest) or DEFAULT_MANIFEST
    database_root = _to_path(args.database_root) or DEFAULT_DATABASE_ROOT
    results_root = _to_path(args.results_root) or DEFAULT_RESULTS_ROOT

    if args.print_downloads:
        from star_gene_data import build_download_commands, summarize_downloads

        paper_filter = args.paper[0] if len(args.paper) == 1 else None
        if len(args.paper) > 1:
            paper_filter = None
        if args.paper and len(args.paper) == 1:
            print(summarize_downloads(paper=paper_filter))
            commands = build_download_commands(
                paper=paper_filter,
                include_large=args.include_large_downloads,
            )
        else:
            print(summarize_downloads())
            commands = build_download_commands(include_large=args.include_large_downloads)

        if commands:
            print('\n[DOWNLOAD COMMANDS]')
            for command in commands:
                print(command)
        else:
            print('\n[DOWNLOAD COMMANDS]')
            print('No safe direct-download commands for the selected paper(s). See notes above.')
        return 0

    manifest = load_manifest(manifest_path)
    if not args.all and not args.paper and not args.target:
        print('[INFO] No --paper/--target filter supplied; checking all manifest targets. Use --all to make this explicit.')

    validator = StarGeneValidator(
        manifest=manifest,
        database_root=database_root,
        results_root=results_root,
        check_only=args.check_only,
        allow_large_download=args.allow_large_download,
        accept_license=args.accept_license,
        score_mode=args.score_mode,
    )
    validator.run(papers=args.paper, targets=args.target)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
