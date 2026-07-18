# Compact Initial Display Range Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Keep automatically selected report windows at no more than 25 rendered variants and no more than 2 kb while avoiding empty midpoint windows in sparse genes.

**Architecture:** Extend `_select_initial_display_range()` with a physical-span cap after the existing neutral 25-variant selection. Preserve the old result when it is already compact; otherwise select a deterministic dense variant window using only coordinates and the gene/region midpoint. The generated report continues to pass the returned range into the existing applied/pending filter and unified gene/GWAS coordinate-domain paths.

**Tech Stack:** Python 3, `unittest`, generated HTML/JavaScript, D3.js, Microsoft Edge Chromium Playwright smoke testing, Git.

---

## File map

- Modify `haplotype_phenotype_analysis.py`: compact initial-range selection only; no scoring changes.
- Modify `test_star_gene_data.py`: Python regression coverage for span capping and deterministic compact-cluster selection.
- Modify `star_gene_validation_record.md`: record the VRN-B1 rerun, exact initial range, browser result, and unchanged scientific conclusion.
- Regenerate ignored/local output `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html` for user inspection.

### Task 1: Add failing compact-window tests

**Files:**
- Modify: `test_star_gene_data.py:29-145`
- Test: `test_star_gene_data.py`

- [ ] **Step 1: Add a test proving an already compact old window is unchanged**

Retain the existing `test_initial_display_range_centers_by_sorted_index_across_large_gaps` assertion:

```python
positions = list(range(1, 21)) + list(range(1000, 1020))
self.assertEqual(
    _select_initial_display_range(
        positions, 1, 2000,
        gene_start=999, gene_end=1001,
        max_variants=25, max_span_bp=2000,
    ),
    (9, 1012),
)
```

- [ ] **Step 2: Add a failing test for a 25-site window wider than 2 kb**

```python
def test_initial_display_range_falls_back_to_compact_cluster(self):
    from haplotype_phenotype_analysis import _select_initial_display_range

    positions = list(range(100, 4100, 100))
    start, end = _select_initial_display_range(
        positions, 1, 5000,
        gene_start=1500, gene_end=2500,
        max_variants=25, max_span_bp=2000,
    )

    selected = [pos for pos in positions if start <= pos <= end]
    self.assertEqual((start, end), (1000, 3000))
    self.assertEqual(len(selected), 21)
    self.assertLessEqual(end - start, 2000)
```

- [ ] **Step 3: Add failing tests for sparse and empty inputs**

```python
def test_initial_display_range_caps_sparse_small_variant_set(self):
    from haplotype_phenotype_analysis import _select_initial_display_range

    self.assertEqual(
        _select_initial_display_range(
            [100, 5000, 9000], 1, 10000,
            gene_start=4000, gene_end=6000,
            max_variants=25, max_span_bp=2000,
        ),
        (5000, 5000),
    )

def test_initial_display_range_uses_fixed_anchor_window_without_variants(self):
    from haplotype_phenotype_analysis import _select_initial_display_range

    self.assertEqual(
        _select_initial_display_range(
            [], 1, 10000,
            gene_start=4000, gene_end=6000,
            max_variants=25, max_span_bp=2000,
        ),
        (4000, 6000),
    )
```

- [ ] **Step 4: Add a failing VRN-B1-like clustered-distribution test**

```python
def test_initial_display_range_chooses_vrn_b1_like_dense_cluster(self):
    from haplotype_phenotype_analysis import _select_initial_display_range

    positions = (
        list(range(5447, 5459))
        + [7596, 7600, 7603, 7605]
        + list(range(12379, 12429))
    )
    start, end = _select_initial_display_range(
        positions, 1, 17867,
        gene_start=4673, gene_end=17867,
        max_variants=25, max_span_bp=2000,
    )

    selected = [pos for pos in positions if start <= pos <= end]
    self.assertEqual((start, end), (12379, 12403))
    self.assertEqual(len(selected), 25)
    self.assertLessEqual(end - start, 2000)
```

- [ ] **Step 5: Run the new tests and verify RED**

Run:

```powershell
python -m unittest `
  test_star_gene_data.StarGeneDataTests.test_initial_display_range_falls_back_to_compact_cluster `
  test_star_gene_data.StarGeneDataTests.test_initial_display_range_caps_sparse_small_variant_set `
  test_star_gene_data.StarGeneDataTests.test_initial_display_range_uses_fixed_anchor_window_without_variants `
  test_star_gene_data.StarGeneDataTests.test_initial_display_range_chooses_vrn_b1_like_dense_cluster
```

Expected: FAIL because `_select_initial_display_range()` does not accept `max_span_bp` and still returns the old broad/full windows.

### Task 2: Implement deterministic span capping

**Files:**
- Modify: `haplotype_phenotype_analysis.py:112-201`
- Test: `test_star_gene_data.py`

- [ ] **Step 1: Extend the helper signature and normalize the span limit**

```python
def _select_initial_display_range(
    variant_positions,
    region_start,
    region_end,
    gene_start=None,
    gene_end=None,
    max_variants=25,
    max_span_bp=2000,
):
    try:
        span_limit = max(1, int(max_span_bp))
    except (TypeError, ValueError):
        span_limit = 2000
```

- [ ] **Step 2: Preserve the old candidate selection but defer its return**

Store the existing output as `candidate_start` and `candidate_end`. Interpret `None` as `lower`/`upper` for the span check:

```python
effective_start = lower if candidate_start is None else candidate_start
effective_end = upper if candidate_end is None else candidate_end
if effective_end - effective_start <= span_limit:
    return candidate_start, candidate_end
```

For `len(positions) <= limit`, set the old candidate to `(None, None)` instead of returning immediately.

- [ ] **Step 3: Add the no-variant fixed anchor fallback**

```python
if not unique_positions:
    width = min(span_limit, upper - lower)
    fixed_start = int(round(anchor - width / 2.0))
    fixed_start = max(lower, min(upper - width, fixed_start))
    fixed_end = fixed_start + width
    return (
        None if fixed_start <= lower else fixed_start,
        None if fixed_end >= upper else fixed_end,
    )
```

- [ ] **Step 4: Enumerate compact windows with deterministic scoring**

```python
compact_best = None
for left_candidate in range(len(unique_positions)):
    for right_candidate in range(left_candidate, len(unique_positions)):
        selected_count = (
            prefix_counts[right_candidate + 1]
            - prefix_counts[left_candidate]
        )
        if selected_count > limit:
            break
        start_pos = unique_positions[left_candidate]
        end_pos = unique_positions[right_candidate]
        physical_span = end_pos - start_pos
        if physical_span > span_limit:
            break
        score = (
            -selected_count,
            abs((start_pos + end_pos) / 2.0 - anchor),
            physical_span,
            start_pos,
        )
        if compact_best is None or score < compact_best[0]:
            compact_best = (score, start_pos, end_pos)
```

Return the chosen coordinates using the existing `None` boundary convention. If no compact window is possible, use the same clamped anchor fallback as the no-variant branch.

- [ ] **Step 5: Run focused tests and verify GREEN**

Run:

```powershell
python -m unittest test_star_gene_data.StarGeneDataTests.test_initial_display_range_falls_back_to_compact_cluster test_star_gene_data.StarGeneDataTests.test_initial_display_range_caps_sparse_small_variant_set test_star_gene_data.StarGeneDataTests.test_initial_display_range_uses_fixed_anchor_window_without_variants test_star_gene_data.StarGeneDataTests.test_initial_display_range_chooses_vrn_b1_like_dense_cluster
```

Expected: `Ran 4 tests ... OK`.

- [ ] **Step 6: Run all initial-range tests**

Run:

```powershell
python -m unittest test_star_gene_data.StarGeneDataTests
```

Expected: all tests pass, including existing duplicate/boundary behavior.

- [ ] **Step 7: Commit the algorithm and tests**

```powershell
git add -- haplotype_phenotype_analysis.py test_star_gene_data.py
git diff --cached --check
git commit -m "Cap automatic report ranges by physical span"
```

### Task 3: Regenerate and verify VRN-B1

**Files:**
- Modify: `star_gene_validation_record.md`
- Regenerate: `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`

- [ ] **Step 1: Regenerate the robust-discovery report**

Run:

```powershell
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Expected: exit 0 and the main HTML timestamp updates.

- [ ] **Step 2: Verify the generated initial range from the HTML**

Read `initialDisplayRange` and assert:

```text
rangeEnd - rangeStart <= 2000
visible sequence columns <= 25
pending range == applied range on load
```

Also confirm the range is selected from the approximately 12.4 kb dense cluster and is not anchored to the literature marker at 13,077.

- [ ] **Step 3: Run the real Edge browser smoke test**

Open the generated `file:///` report with Playwright and verify:

```text
1. No page-level JavaScript errors.
2. Initial gene-structure and GWAS axes use the new compact range.
3. Apply is disabled on load because pending == applied.
4. Reset restores the same compact initial range.
5. Manual range edits remain pending until Apply.
```

- [ ] **Step 4: Verify scientific outputs are unchanged**

Read `haplotype_scores.json` and the existing post-hoc marker audit. Expected:

```text
top haplotype = Hap6
total score = 1.2016
sample count = 3
837 bp insertion post-hoc match = yes
```

- [ ] **Step 5: Update the validation record**

Append a dated section to `star_gene_validation_record.md` containing the target, literature variant, data source, `robust_discovery` mode, exact command/output directory, exact new initial range and visible count, unchanged Hap6 result, reliability limit, browser checks, and any blocked reason. State explicitly that the literature marker did not influence initial-range selection.

- [ ] **Step 6: Commit the regenerated-validation record**

```powershell
git add -- star_gene_validation_record.md
git diff --cached --check
git commit -m "Record compact VRN-B1 initial range validation"
```

### Task 4: Final regression, review, and push

**Files:**
- Verify: `haplotype_phenotype_analysis.py`
- Verify: `test_star_gene_data.py`
- Verify: `star_gene_validation_record.md`

- [ ] **Step 1: Run syntax and full unit tests**

```powershell
python -m py_compile haplotype_phenotype_analysis.py test_star_gene_data.py
python -m unittest test_star_gene_data
git diff --check
```

Expected: compilation succeeds, all tests pass, and no whitespace errors are reported.

- [ ] **Step 2: Run the project smoke test with a bounded timeout**

```powershell
python run_rice_test.py
```

Expected: pass when local indexed data/pysam permit it. If it remains silent on the known pure-Python large-VCF path, stop after the agreed bounded timeout and record the limitation rather than claiming success.

- [ ] **Step 3: Request independent code review**

Ask the reviewer to inspect only the compact-range commits for phenotype/literature leakage, deterministic tie-breaking, boundary behavior, duplicate handling, and test adequacy. Fix Important/Critical findings and rerun the focused/full tests.

- [ ] **Step 4: Confirm commit scope**

```powershell
git status --short --branch
git log -3 --oneline --decorate
```

Expected: only the intended source/test/record commits are new; unrelated historical dirty files remain untouched.

- [ ] **Step 5: Push once**

```powershell
git push origin work/star-gene-validation
```

Expected: branch updates on GitHub. If the network fails, report it once without repeated retries.
