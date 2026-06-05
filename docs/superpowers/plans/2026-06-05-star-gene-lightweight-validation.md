# Star Gene Lightweight Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a lightweight, reproducible path for validating HaplotypeScorer against star genes/loci from the rice, maize, and wheat papers without forcing large whole-dataset downloads.

**Architecture:** Keep `star_gene_validation.py` as the checker/runner and add a focused helper module for external data instructions and small-download metadata. Use the existing precomputed-database interface for actual analysis rather than changing the core analyzer.

**Tech Stack:** Python 3, pandas-compatible CSV/TSV inputs, existing `HaplotypePhenotypeAnalyzer`, JSON manifest metadata, standard-library download utilities.

---

## File Structure

- Modify `star_gene_manifest.json`: correct rice target traits and add explicit lightweight data-source metadata.
- Create `star_gene_data.py`: expose small-file download metadata, data-source reporting, and command generation without automatically downloading large files.
- Modify `star_gene_validation.py`: add CLI options to print download instructions and recognize lightweight database-ready inputs.
- Create `test_star_gene_data.py`: unit tests for metadata, large-download safeguards, and generated commands.
- Create or update `docs/star_gene_validation.md`: document the lightweight workflow and exact commands.
- Update `task_plan.md`, `findings.md`, and `progress.md`: record implementation progress.

## Tasks

### Task 1: Add Download Metadata Tests

**Files:**
- Create: `test_star_gene_data.py`

- [x] **Step 1: Write failing tests**

Create tests that expect:
- rice Figshare source is marked large and no auto-download command is produced by default;
- maize MaizeGo SV and phenotype files are marked small and have direct URLs;
- wheat WWWG2B/Earlham entries are instruction-only until concrete object names are confirmed.

- [x] **Step 2: Run tests and verify RED**

Run: `python -m unittest test_star_gene_data.py -v`

Expected: import error for `star_gene_data`.

### Task 2: Implement Data Metadata Helper

**Files:**
- Create: `star_gene_data.py`

- [x] **Step 1: Implement metadata structures and safe command generation**

Add `DataFile` and functions:
- `iter_data_files(paper=None)`
- `summarize_downloads(paper=None)`
- `build_download_commands(paper=None, include_large=False)`

- [x] **Step 2: Run tests and verify GREEN**

Run: `python -m unittest test_star_gene_data.py -v`

Expected: all tests pass.

### Task 3: Wire CLI Download Instructions

**Files:**
- Modify: `star_gene_validation.py`
- Test: `test_star_gene_data.py`

- [x] **Step 1: Add CLI tests**

Add tests for `--print-downloads` behavior by calling helper functions directly, avoiding network.

- [x] **Step 2: Implement CLI option**

Add `--print-downloads` and `--include-large-downloads` to print reproducible instructions and exit before normal validation.

- [x] **Step 3: Run focused tests**

Run: `python -m unittest test_star_gene_data.py -v`

Expected: all tests pass.

### Task 4: Correct Manifest Positive Controls

**Files:**
- Modify: `star_gene_manifest.json`

- [x] **Step 1: Correct rice traits**

Set `OsMADS22` primary trait to panicle number aliases and keep `OsFTL1` heading date. Do not invent coordinates.

- [x] **Step 2: Add lightweight path notes**

Add notes that maize has direct small MaizeGo files and wheat requires concrete WWWG2B/Earlham object names before hard-coded downloads.

- [x] **Step 3: Run validation check**

Run: `python run_star_gene_validation.py --check-only --all --no-download`

Expected: still safely reports missing local data and unresolved coordinates.

### Task 5: Documentation and Verification

**Files:**
- Modify: `docs/star_gene_validation.md`
- Modify: `task_plan.md`, `findings.md`, `progress.md`

- [x] **Step 1: Document commands**

Include:
- `python run_star_gene_validation.py --print-downloads --paper maize2019`
- `python run_star_gene_validation.py --print-downloads --paper rice2024 --include-large-downloads`
- `python run_star_gene_validation.py --check-only --all --no-download`

- [x] **Step 2: Run required tests**

Run:
- `python -m unittest test_star_gene_data.py -v`
- `python run_star_gene_validation.py --check-only --all --no-download`
- `python run_rice_test.py`

Observed: focused tests passed, star-gene check remained safe, and rice quick test exited 0 with one existing no-phenotype-overlap status.

- [x] **Step 3: Git commit and push**

Stage only files changed for this task, commit, and push. If push fails due to network, record that instead of retrying indefinitely.
