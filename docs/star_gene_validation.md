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

## Lightweight download instructions

Print data-source notes and copyable commands:

```bash
python run_star_gene_validation.py --print-downloads --paper maize2019
python run_star_gene_validation.py --print-downloads --paper rice2024 --include-large-downloads
python run_star_gene_validation.py --print-downloads --paper wheat2024
```

`--print-downloads` never downloads data. It prints:

- direct PowerShell commands only for small files with confirmed URLs, currently the MaizeGo SV/pSV and BLUP trait files;
- manual large-download notes for Rice Figshare `NAM_variations` (`10.6084/m9.figshare.19166475`, about 12.64 GB);
- instruction-only notes for Wheat WWWG2B/Earlham until the exact target-region VCF object names are confirmed.

The recommended first pass is:

1. Use `--print-downloads --paper maize2019`.
2. Download the small MaizeGo files into `external_data/maize_natgenet_2019/maizego/`.
3. Inspect whether the qHKW1/ZmBAM1d 8.9 kb indel has sample-level allele calls.
4. Convert the confirmed marker/phenotype table into the existing per-target database format before running `--run-analysis`.

## Build a precomputed database from marker tables

When a paper source provides a small sample-by-marker matrix, convert it with:

```bash
python build_star_gene_database.py \
  --marker-matrix external_data/maize_natgenet_2019/qHKW1_markers.tsv \
  --phenotype-table external_data/maize_natgenet_2019/hundred_kernel_weight.tsv \
  --output-root star_gene_database/maize_natgenet_2019 \
  --target-id qHKW1 \
  --chrom chr1 \
  --start 100 \
  --end 100 \
  --phenotype-column HKW \
  --marker-column qHKW1_8_9kb_indel \
  --marker-position qHKW1_8_9kb_indel=100
```

Use real marker names and coordinates from the paper data. The command writes:

- `gene_info.json`
- `haplotype_data.csv`
- `haplotype_samples.csv`
- `phenotype_data.csv`
- `variant_info.csv`

Those files are the same precomputed database format already consumed by `HaplotypePhenotypeAnalyzer`.

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
- Do not hard-code rice or wheat coordinates until the paper data files confirm assembly/build and coordinate system.
