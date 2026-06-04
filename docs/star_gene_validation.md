# Star-gene positive-control validation

This repository now includes a minimal, conservative framework for validating `HaplotypeScorer` against high-confidence genes/loci from three crop genomics papers.

The first implementation is intentionally **check-only by default**. It records the target papers and loci, verifies whether local data files are present, checks obvious VCF/phenotype compatibility when files exist, and writes a uniform summary. It does **not** download large datasets or guess unresolved coordinates.

## Files

- [`star_gene_manifest.json`](../star_gene_manifest.json): paper/target manifest, expected local paths, trait aliases, coordinate-resolution status, and data-source notes.
- [`star_gene_validation.py`](../star_gene_validation.py): common validation library and CLI implementation.
- [`run_star_gene_validation.py`](../run_star_gene_validation.py): generic CLI entry point.
- [`run_rice2024_star_genes.py`](../run_rice2024_star_genes.py): thin wrapper for the Rice Science 2024 targets.
- [`run_maize2019_star_gene.py`](../run_maize2019_star_gene.py): thin wrapper for the Maize Nature Genetics 2019 qHKW1/ZmBAM1d targets.
- [`run_wheat2024_star_genes.py`](../run_wheat2024_star_genes.py): thin wrapper for Wheat Nature 2024 WatSeq targets.

## Default safe check

```bash
python run_star_gene_validation.py --check-only --all --no-download
```

Outputs are written under ignored directories:

- `star_gene_results/validation_summary.csv`
- `star_gene_results/validation_summary.json`

The summary uses one row per target check and includes status fields such as:

- `missing_required_files`
- `requires_coordinate_resolution`
- `ready_for_analysis_check_only`
- `skipped_*`
- `analyzed`

## Paper-specific checks

```bash
python run_rice2024_star_genes.py --check-only --no-download
python run_maize2019_star_gene.py --check-only --no-download
python run_wheat2024_star_genes.py --check-only --no-download
```

## Analysis mode

Analysis mode is available but guarded by the manifest:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target RHT8
```

A target will only run if the manifest has resolved `coordinates` (`chrom`, `start`, `end`) and the expected local files exist. For the current manifest, the planned targets intentionally keep `requires_coordinate_resolution: true` until paper-specific coordinates/builds are confirmed.

## Haplotype score export

[`haplotype_phenotype_analysis.py`](../haplotype_phenotype_analysis.py) now writes `haplotype_scores.json` in each per-gene result directory when the integrated report computes HaplotypeScorer outputs. The analyzer return value also includes:

- `haplotype_scores`
- `haplotype_score_json_path`

The JSON preserves the existing `HaplotypeScorer.score_all()` structure and does not alter the scoring formula.

## Data policy

- No large downloads are implemented in this minimal framework.
- `--allow-large-download` is only a future-facing placeholder and is ignored unless paired with `--accept-license`; even then, this version still reports that downloading is not implemented.
- `external_data/`, `star_gene_database/`, and `star_gene_results/` are ignored by git to avoid committing large or generated files.
