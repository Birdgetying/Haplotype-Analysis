# Discovery Key Sites Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a literature-free report panel that ranks likely driver sites distinguishing the current top-scored haplotype from other haplotypes.

**Architecture:** Reuse the existing discovery candidate panel data path. Compute key-site rows from `hap_sample_df`, `Haplotype_Seq`, display-position mappings, score results, phenotype values, and local variant metadata; pre-render HTML per score mode and phenotype so Original/Robust switching stays consistent.

**Tech Stack:** Python, pandas, unittest, generated standalone HTML/JavaScript.

---

### Task 1: Add Red Tests For Key-Site Ranking

**Files:**
- Modify: `D:/Desktop/project1/test_star_gene_data.py`

- [ ] Add a test that imports `_build_top_haplotype_key_site_rows` and `_render_top_haplotype_key_sites`, builds three haplotypes where `Hap1` is top-scored, and expects only sites where `Hap1` differs from other haplotypes to be ranked.
- [ ] Add assertions for top allele, other allele distribution, phenotype contrast, annotation/type fields, reliability flag, and no literature wording.
- [ ] Add integrated-report source assertions for `top_haplotype_key_sites_html`, `top-hap-key-site`, and `highlightKeySiteColumns`.
- [ ] Run `python -m unittest test_star_gene_data.StarGeneDataTests.test_render_discovery_candidate_panels_are_literature_free -v` and confirm failure caused by the missing helper/panel.

### Task 2: Implement Key-Site Helper And Renderer

**Files:**
- Modify: `D:/Desktop/project1/haplotype_phenotype_analysis.py`

- [ ] Implement `_build_top_haplotype_key_site_rows(...)`.
- [ ] Rank positions by top-hap specificity, phenotype contrast magnitude, GWAS support, annotation priority, missing-rate/MAF reliability, and top-hap sample count.
- [ ] Implement `_render_top_haplotype_key_sites(...)` as a compact table with data-position attributes.
- [ ] Extend `_build_discovery_candidate_panel_data(...)` and `_select_discovery_candidate_panel(...)` to include `top_haplotype_key_sites_html`.
- [ ] Run the red test again and confirm it passes.

### Task 3: Wire HTML Panel And Sequence Highlighting

**Files:**
- Modify: `D:/Desktop/project1/haplotype_phenotype_analysis.py`

- [ ] Add the new panel below the discovery candidate list in the Local Candidate Evidence block.
- [ ] Add CSS for `.top-hap-key-site-table`, `.key-site-flag`, and `.key-site-highlight`.
- [ ] Add `data-pos` and `data-key-site` attributes to sequence header/cells.
- [ ] Add JavaScript `highlightKeySiteColumns(panel)` and call it from `updateDiscoveryCandidatePanels()`.
- [ ] Run `python -m unittest test_star_gene_data.py -v`.

### Task 4: Regenerate Reports And Verify

**Files:**
- Modify generated outputs under `D:/Desktop/project1/star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP*`
- Modify: `D:/Desktop/project1/star_gene_validation_record.md`

- [ ] Run `python -m py_compile haplotype_phenotype_analysis.py star_gene_validation.py`.
- [ ] Run `python run_rice_test.py`.
- [ ] Run `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP`.
- [ ] Run `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --score-mode robust_discovery`.
- [ ] Check the current browser report shows the new Candidate Key Sites panel and that Original/Robust switching updates it.
- [ ] Update `star_gene_validation_record.md` with the feature/rerun note.
- [ ] Stage only relevant files, commit, and try one push.
