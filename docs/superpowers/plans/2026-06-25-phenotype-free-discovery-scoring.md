# Phenotype-Free Discovery Scoring Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove validation phenotype leakage from robust haplotype discovery scoring.

**Architecture:** Keep `HaplotypeScorer` as the scoring boundary. Add local site-weighting helpers for annotation, distance, quality, and external evidence; use phenotype only in validation/report outputs.

**Tech Stack:** Python, pandas, numpy, scipy, unittest.

---

### Task 1: Lock Phenotype-Free Robust Behavior

**Files:**
- Modify: `test_star_gene_data.py`
- Modify: `haplotype_phenotype_analysis.py`

- [ ] Add a regression test showing that robust discovery totals and selected functional positions do not change when the phenotype column is inverted.
- [ ] Run the targeted test and confirm it fails against the current phenotype-leaking implementation.
- [ ] Remove phenotype-derived robust scoring inputs from `HaplotypeScorer`.
- [ ] Run the targeted test and confirm it passes.

### Task 2: Preserve Post Hoc Validation

**Files:**
- Modify: `test_star_gene_data.py`
- Modify: `haplotype_phenotype_analysis.py`

- [ ] Add a regression test showing that `expected_direction` no longer changes the discovery `score_axis`.
- [ ] Keep `direction_score` and `directional_total` in output as validation/audit fields, but set active `score_axis` to `total`.
- [ ] Run targeted tests for robust scoring and existing literature audit behavior.

### Task 3: Document Star-Gene Validation Change

**Files:**
- Modify: `star_gene_validation_record.md`

- [ ] Add an entry explaining that robust discovery is now phenotype-free and PostGWAS-like site evidence is local code, not a live web call.
- [ ] Run quick syntax/tests after code changes.
- [ ] Stage only relevant files and commit.
