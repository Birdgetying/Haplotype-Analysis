# Functional Sub-Haplotype Robust Scoring Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a generic functional sub-haplotype discovery layer and validate it on TaGW2-B1 full-region, VRN-D1-Kiss2014, and Rht-Zanke2014.

**Architecture:** Extend `HaplotypeScorer` to emit `functional_haplotype_groups` beside `core_haplotype_groups`. Extend validation summary and literature audit to read functional-group fields without feeding literature variants into scoring.

**Tech Stack:** Python, pandas, scipy, existing unittest suite, existing star-gene validation CLI.

---

### Task 1: Scorer Functional Groups

**Files:**
- Modify: `haplotype_phenotype_analysis.py`
- Test: `test_star_gene_data.py`

- [ ] Write a failing test that constructs regional haplotypes split by background SNPs and asserts `functional_haplotype_groups` groups by selected functional positions.
- [ ] Run the targeted test and confirm it fails because the field is absent.
- [ ] Add functional-position selection using local annotation, phenotype signal, EB effect, MAF stability, missingness, and LD pruning.
- [ ] Emit group sequence, selected positions, mean score, rank score, mean phenotype, sample count, representative haplotype, and member haplotypes.
- [ ] Run the targeted test and confirm it passes.

### Task 2: Summary and Audit

**Files:**
- Modify: `star_gene_validation.py`
- Modify: `star_gene_literature_audit.py`
- Test: `test_star_gene_data.py`

- [ ] Write failing tests for top/directional functional group selection and literature audit matching through functional groups.
- [ ] Add summary columns for raw and direction-aware functional groups.
- [ ] Add audit columns/status for functional group post hoc matching.
- [ ] Run targeted tests and confirm they pass.

### Task 3: Three-Target Validation

**Files:**
- Modify: `star_gene_validation_record.md`
- Modify: `findings.md`
- Modify: `progress.md`

- [ ] Run `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --target VRN-D1-Kiss2014 --target Rht-Zanke2014 --score-mode robust_discovery`.
- [ ] Extract summary and audit results for all three targets.
- [ ] Update validation record with target, source, score mode, command, output, top functional group, literature match status, and limitations.
- [ ] Run `python -m py_compile haplotype_phenotype_analysis.py star_gene_validation.py star_gene_literature_audit.py`.
- [ ] Run `python -m unittest test_star_gene_data.py -v`.
- [ ] Run `python run_star_gene_validation.py --check-only --paper wheat2024 --target TaGW2-B1-remoteSNP --target VRN-D1-Kiss2014 --target Rht-Zanke2014 --no-download`.
- [ ] Commit and push only relevant source, test, plan, and documentation files.
