# Explicit Range Apply Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make display-range edits cheap previews and update the report only after the user explicitly clicks Apply Range.

**Architecture:** Separate pending range state from `currentFilter`, which remains the applied range consumed by GWAS, gene structure, sequence table, and LD rendering. Input, slider, and Clear handlers update only pending state and dirty UI; one commit function canonicalizes the pending range, writes `currentFilter`, and invokes `applyFilters()` once.

**Tech Stack:** Python-generated HTML, vanilla JavaScript/CSS, Python `unittest`, Node.js JavaScript checks, Playwright browser verification.

---

### Task 1: Add failing explicit-apply regression tests

**Files:**
- Modify: `test_star_gene_data.py:1825-1935`

- [ ] **Step 1: Replace real-time range assertions with explicit-apply assertions**

Add static checks for `rangeApplyBtn`, `rangePendingStatus`, `pendingRangeFilter`, and `commitRangeFilter()`. Assert that `updateRangeFilterFromInputs`, `updateRangeFilterFromSliders`, and `clearRangeFilter` do not call `applyFilters()`, while `commitRangeFilter()` contains exactly one call path.

```python
self.assertIn('id="rangeApplyBtn"', integrated_block)
self.assertIn('id="rangePendingStatus"', integrated_block)
self.assertIn("var pendingRangeFilter = { start: null, end: null };", integrated_block)
self.assertIn("function commitRangeFilter()", integrated_block)
self.assertNotIn("function scheduleRangeFilterApply()", integrated_block)

for function_name in (
    "updateRangeFilterFromInputs",
    "updateRangeFilterFromSliders",
    "clearRangeFilter",
):
    function_block = extract_js_function(function_name)
    self.assertNotIn("applyFilters();", function_block)

commit_block = extract_js_function("commitRangeFilter")
self.assertEqual(commit_block.count("applyFilters();"), 1)
```

- [ ] **Step 2: Add executable JavaScript checks for pending and commit behavior**

Extract the pure range helpers and commit function into the existing Node test. Use minimal DOM elements and an `applyCount` counter to prove edits do not apply, one click applies once, and clicking again without changes does not apply again.

```javascript
var applyCount = 0;
function applyFilters() { applyCount += 1; }

updateRangeFilterFromInputs('start');
assertEqual(applyCount, 0, 'number edit stays pending');
clearRangeFilter();
assertEqual(applyCount, 0, 'clear stays pending');
commitRangeFilter();
assertEqual(applyCount, 1, 'commit applies once');
commitRangeFilter();
assertEqual(applyCount, 1, 'unchanged commit is skipped');
```

- [ ] **Step 3: Run the focused tests and verify RED**

Run:

```powershell
python -m unittest test_star_gene_data.StarGeneDataTests.test_integrated_report_has_position_range_filter_controls test_star_gene_data.StarGeneDataTests.test_display_range_slider_executes_overlap_and_timer_logic -v
```

Expected: FAIL because `rangeApplyBtn`, `pendingRangeFilter`, and `commitRangeFilter()` do not exist and the existing handlers still schedule or invoke filtering.

### Task 2: Implement pending range state and explicit Apply Range UI

**Files:**
- Modify: `haplotype_phenotype_analysis.py:9932-9951`
- Modify: `haplotype_phenotype_analysis.py:10077-10092`
- Modify: `haplotype_phenotype_analysis.py:11058-11415`
- Modify: `haplotype_phenotype_analysis.py:11474-11520`
- Modify: `haplotype_phenotype_analysis.py:11910-11925`

- [ ] **Step 1: Add Apply Range and pending status controls**

Place the action next to Clear and provide a concise dirty-state label.

```html
<button class="range-clear-btn" type="button" onclick="clearRangeFilter()">Clear</button>
<button class="range-apply-btn" id="rangeApplyBtn" type="button" onclick="commitRangeFilter()" disabled>Apply Range</button>
<span class="range-pending-status" id="rangePendingStatus" hidden>Unapplied changes</span>
```

Add light-theme styles for enabled, disabled, and pending states. Preserve the existing responsive wrapping behavior.

- [ ] **Step 2: Introduce independent pending state and pure range reads**

Replace `rangeFilterFrame` with pending state.

```javascript
var pendingRangeFilter = { start: null, end: null };
```

Make `readRangeFilter(changedHandle)` return normalized control values without modifying `currentFilter`. Add canonicalization and equality helpers so full-boundary values are equivalent to `null`.

```javascript
function canonicalizeRangeFilter(range) {
    var bounds = getRangeBounds();
    return {
        start: range.start === null || range.start <= bounds.min ? null : range.start,
        end: range.end === null || range.end >= bounds.max ? null : range.end,
        min: bounds.min,
        max: bounds.max
    };
}

function rangeFiltersEqual(a, b) {
    return a.start === b.start && a.end === b.end;
}
```

- [ ] **Step 3: Make all edit paths pending-only**

Add a helper that stores pending values, synchronizes controls, and updates the dirty indicator without filtering.

```javascript
function setPendingRangeFilter(range) {
    pendingRangeFilter = canonicalizeRangeFilter(range);
    syncRangeControls(pendingRangeFilter);
    updateRangePendingState();
}
```

Use it from pointer dragging, native slider input, number input, and Clear. Remove `scheduleRangeFilterApply()` entirely. Pointer up/cancel only releases capture.

- [ ] **Step 4: Add the single explicit commit path**

```javascript
function commitRangeFilter() {
    var nextRange = canonicalizeRangeFilter(readRangeFilter());
    var appliedRange = canonicalizeRangeFilter(currentFilter);
    setPendingRangeFilter(nextRange);
    if (rangeFiltersEqual(nextRange, appliedRange)) return;
    currentFilter.rangeStart = nextRange.start;
    currentFilter.rangeEnd = nextRange.end;
    updateRangePendingState();
    applyFilters();
}
```

The dirty-state helper enables the button and reveals the status only when pending and applied ranges differ.

- [ ] **Step 5: Remove hidden implicit range commits**

Delete `syncRangeControls(readRangeFilter())` from the beginning of `applyFilters()` so changing another filter or exporting cannot silently commit pending range edits. In reset and cross-window message handlers, explicitly update both `currentFilter` and `pendingRangeFilter` before syncing controls, because those paths intentionally apply complete filter state.

- [ ] **Step 6: Run focused tests and verify GREEN**

Run the same two focused tests from Task 1.

Expected: both PASS; the Node counter confirms no apply during edits and exactly one apply on a changed commit.

- [ ] **Step 7: Run the full unit suite**

```powershell
python -m py_compile haplotype_phenotype_analysis.py test_star_gene_data.py star_gene_validation.py
python -m unittest test_star_gene_data -v
```

Expected: compilation exits 0 and all tests pass.

### Task 3: Regenerate and verify the VRN-B1 report

**Files:**
- Modify: `star_gene_validation_record.md`
- Regenerate: `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`

- [ ] **Step 1: Rerun the existing robust discovery target**

```powershell
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Expected: exit 0 and the HTML contains `rangeApplyBtn`, `rangePendingStatus`, and `commitRangeFilter()`.

- [ ] **Step 2: Run Chromium and Firefox interaction checks**

Wrap `window.applyFilters` with a counter after load. Verify slider dragging, number entry, and Clear leave the count at zero; `Apply Range` increments it once; a second unchanged click does not increment it; the rendered range updates only after the first click.

- [ ] **Step 3: Update the validation record without changing scientific conclusions**

Append the rerun command, output path, unchanged robust top haplotype and literature-marker match, test evidence, and a note that this rerun changes only report interaction timing.

- [ ] **Step 4: Run final diff checks**

```powershell
git diff --check -- haplotype_phenotype_analysis.py test_star_gene_data.py star_gene_validation_record.md
git status --short -- haplotype_phenotype_analysis.py test_star_gene_data.py star_gene_validation_record.md
```

Expected: no whitespace errors and only the three intended source/record files are selected for the implementation commit; generated report output remains untracked or ignored according to existing repository rules.

### Task 4: Review, commit, and push

**Files:**
- Review: `haplotype_phenotype_analysis.py`
- Review: `test_star_gene_data.py`
- Review: `star_gene_validation_record.md`

- [ ] **Step 1: Review the implementation against the approved design**

Check that pending edits never affect `inDisplayRange()`, exports, other filter changes, or external filter messages until an intentional commit path runs.

- [ ] **Step 2: Stage only intended files and commit**

```powershell
git add -- haplotype_phenotype_analysis.py test_star_gene_data.py star_gene_validation_record.md
git diff --cached --check
git commit -m "Make display range apply explicit"
```

- [ ] **Step 3: Push once**

```powershell
git push origin work/star-gene-validation
```

If the network fails, record the local commit and error without repeated retries.
