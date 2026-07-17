# Initial Display Range and Zoom Layout Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Open long reports on a neutral 25-variant window, keep the main viewport full-width at every Zoom level, and make pending-range preview identical to the applied range.

**Architecture:** A pure Python helper selects the initial coordinate window from variant positions and gene midpoint without phenotype or literature input. Generated JavaScript stores this default separately from pending/applied range state, uses CSS `zoom` instead of `transform: scale()`, preserves the main-data scroll center, and redraws alignment-dependent overlays after zoom.

**Tech Stack:** Python, pandas-generated HTML, vanilla JavaScript/CSS, `unittest`, Node.js JavaScript execution tests, Playwright browser verification.

---

### Task 1: Select a neutral initial 25-variant window

**Files:**
- Modify: `test_star_gene_data.py`
- Modify: `haplotype_phenotype_analysis.py`

- [ ] **Step 1: Add failing pure-function tests**

Import `_select_initial_display_range` and add tests proving full-range fallback, midpoint selection, boundary clamping, and missing-gene fallback.

```python
def test_initial_display_range_uses_full_range_for_small_variant_sets(self):
    self.assertEqual(
        _select_initial_display_range([10, 20, 30], 1, 100, 20, 80),
        (None, None),
    )

def test_initial_display_range_caps_midpoint_window_at_25_variants(self):
    positions = list(range(100, 4100, 100))
    start, end = _select_initial_display_range(
        positions, 1, 5000, gene_start=1500, gene_end=2500, max_variants=25
    )
    selected = [pos for pos in positions if start <= pos <= end]
    self.assertEqual(len(selected), 25)
    self.assertLessEqual(start, 2000)
    self.assertGreaterEqual(end, 2000)

def test_initial_display_range_uses_region_midpoint_without_gene_coordinates(self):
    positions = list(range(1, 41))
    start, end = _select_initial_display_range(
        positions, 1, 40, gene_start=None, gene_end=None, max_variants=25
    )
    self.assertEqual(len([p for p in positions if start <= p <= end]), 25)
    self.assertLessEqual(start, 20)
    self.assertGreaterEqual(end, 20)

def test_initial_display_range_clamps_window_at_variant_boundaries(self):
    positions = list(range(1, 41))
    start, end = _select_initial_display_range(
        positions, 1, 40, gene_start=1, gene_end=1, max_variants=25
    )
    self.assertEqual((start, end), (1, 25))
```

- [ ] **Step 2: Run the focused tests and verify RED**

```powershell
python -m unittest test_star_gene_data.StarGeneDataTests.test_initial_display_range_uses_full_range_for_small_variant_sets test_star_gene_data.StarGeneDataTests.test_initial_display_range_caps_midpoint_window_at_25_variants test_star_gene_data.StarGeneDataTests.test_initial_display_range_uses_region_midpoint_without_gene_coordinates test_star_gene_data.StarGeneDataTests.test_initial_display_range_clamps_window_at_variant_boundaries -v
```

Expected: import or attribute failure because `_select_initial_display_range` does not exist.

- [ ] **Step 3: Implement the pure selector**

Add near the existing report helpers:

```python
def _select_initial_display_range(
    variant_positions,
    region_start,
    region_end,
    gene_start=None,
    gene_end=None,
    max_variants=25,
):
    try:
        lower = int(region_start)
        upper = int(region_end)
    except (TypeError, ValueError):
        return None, None
    if upper < lower:
        lower, upper = upper, lower
    normalized_positions = set()
    for pos in variant_positions or []:
        try:
            numeric_pos = int(pos)
        except (TypeError, ValueError):
            continue
        if lower <= numeric_pos <= upper:
            normalized_positions.add(numeric_pos)
    positions = sorted(normalized_positions)
    limit = max(1, int(max_variants or 25))
    if len(positions) <= limit:
        return None, None
    if gene_start is not None and gene_end is not None:
        anchor = (int(gene_start) + int(gene_end)) / 2.0
    else:
        anchor = (lower + upper) / 2.0
    center_index = min(range(len(positions)), key=lambda i: abs(positions[i] - anchor))
    first_index = max(0, min(len(positions) - limit, center_index - limit // 2))
    selected = positions[first_index:first_index + limit]
    return selected[0], selected[-1]
```

- [ ] **Step 4: Run the focused tests and verify GREEN**

Run the command from Step 2. Expected: all three tests PASS.

- [ ] **Step 5: Commit the selector**

```powershell
git add -- haplotype_phenotype_analysis.py test_star_gene_data.py
git commit -m "Add neutral initial display range selector"
```

### Task 2: Generate and restore the initial range consistently

**Files:**
- Modify: `test_star_gene_data.py:1820-2040`
- Modify: `haplotype_phenotype_analysis.py:9500-9605`
- Modify: `haplotype_phenotype_analysis.py:10900-11570`
- Modify: `haplotype_phenotype_analysis.py:12870-12935`

- [ ] **Step 1: Add failing generated-JavaScript tests**

Extend the integrated-report checks with:

```python
self.assertIn("var initialDisplayRange = {", integrated_block)
self.assertIn("rangeStart: initialDisplayRange.start", integrated_block)
self.assertIn("rangeEnd: initialDisplayRange.end", integrated_block)
self.assertIn("setPendingRangeFilter(initialDisplayRange);", integrated_block)
self.assertIn("canonicalizeRangeFilter(pendingRangeFilter)", commit_block)
self.assertNotIn("canonicalizeRangeFilter(readRangeFilter())", commit_block)
```

Extend the Node script to reproduce the reviewed crossing case:

```javascript
elements.rangeEndInput.value = '80';
updateRangeFilterFromInputs('end');
elements.rangeStartInput.value = '90';
updateRangeFilterFromInputs('start');
assertEqual(pendingRangeFilter, {start: 80, end: 80}, 'crossed input preview');
commitRangeFilter();
assertEqual(
    {start: currentFilter.rangeStart, end: currentFilter.rangeEnd},
    {start: 80, end: 80},
    'commit uses pending preview'
);
```

- [ ] **Step 2: Run the two focused integrated-report tests and verify RED**

```powershell
python -m unittest test_star_gene_data.StarGeneDataTests.test_integrated_report_has_position_range_filter_controls test_star_gene_data.StarGeneDataTests.test_display_range_slider_executes_overlap_and_timer_logic -v
```

Expected: FAIL because initial state remains full range and `commitRangeFilter()` rereads crossed inputs.

- [ ] **Step 3: Inject initial-range constants into generated HTML**

After `display_positions` is finalized, calculate:

```python
initial_range_start, initial_range_end = _select_initial_display_range(
    display_positions,
    region_start,
    region_end,
    g_start,
    g_end,
    max_variants=25,
)
initial_display_range_json = _script_json_dumps({
    "start": initial_range_start,
    "end": initial_range_end,
})
```

Add `var initialDisplayRange = {initial_display_range_json};`, initialize both range states from it, and replace the placeholder with `initial_display_range_json` during final HTML substitution.

- [ ] **Step 4: Make Apply and Reset use the intended state**

Use pending state as the single Apply source:

```javascript
function commitRangeFilter() {
    var nextRange = canonicalizeRangeFilter(pendingRangeFilter);
    var appliedRange = canonicalizeRangeFilter(currentFilter);
    setPendingRangeFilter(nextRange);
    if (rangeFiltersEqual(nextRange, appliedRange)) return;
    currentFilter.rangeStart = nextRange.start;
    currentFilter.rangeEnd = nextRange.end;
    updateRangePendingState();
    applyFilters();
    scheduleAppliedRangeCentering();
}
```

Reset restores the generated default:

```javascript
currentFilter = {
    maf: 0.05,
    missingRate: 0.2,
    rangeStart: initialDisplayRange.start,
    rangeEnd: initialDisplayRange.end
};
setPendingRangeFilter(initialDisplayRange);
```

On DOMContentLoaded, call `setPendingRangeFilter(initialDisplayRange)` before the first `applyFilters()`.

- [ ] **Step 5: Run focused tests and verify GREEN**

Run the command from Step 2. Expected: both tests PASS, including the `90/80 -> 80/80` commit assertion.

- [ ] **Step 6: Commit initial-state behavior**

```powershell
git add -- haplotype_phenotype_analysis.py test_star_gene_data.py
git commit -m "Initialize reports with a bounded display range"
```

### Task 3: Keep the viewport full-width during Zoom

**Files:**
- Modify: `test_star_gene_data.py:1788-1865`
- Modify: `haplotype_phenotype_analysis.py:9685-9998`
- Modify: `haplotype_phenotype_analysis.py:10713-10732`
- Modify: `haplotype_phenotype_analysis.py:11568-11960`

- [ ] **Step 1: Add failing zoom implementation assertions**

```python
self.assertIn("zc.style.zoom = String(cz / 100)", integrated_block)
self.assertIn("zc.style.transform = 'none'", integrated_block)
self.assertNotIn("zc.style.transform = 'scale('", integrated_block)
self.assertIn("zoom: 1 !important", integrated_block)
self.assertIn("function captureMainDataCenter()", integrated_block)
self.assertIn("function restoreMainDataCenter(centerRatio)", integrated_block)
self.assertIn("scheduleConnectorRedraw();", integrated_block)
self.assertIn("scheduleLDTriangleRedraw();", integrated_block)
```

- [ ] **Step 2: Run the layout test and verify RED**

```powershell
python -m unittest test_star_gene_data.StarGeneDataTests.test_integrated_report_uses_right_sidebar_workbench_without_losing_controls -v
```

Expected: FAIL because `applyZoom()` still uses `transform: scale()`.

- [ ] **Step 3: Replace report-level transform scaling with CSS zoom**

```javascript
function captureMainDataCenter() {
    var section = document.querySelector('.main-data-section');
    if (!section || section.scrollWidth <= 0) return 0.5;
    return (section.scrollLeft + section.clientWidth / 2) / section.scrollWidth;
}

function restoreMainDataCenter(centerRatio) {
    var section = document.querySelector('.main-data-section');
    if (!section) return;
    section.scrollLeft = Math.max(
        0,
        centerRatio * section.scrollWidth - section.clientWidth / 2
    );
}

function applyZoom() {
    var centerRatio = captureMainDataCenter();
    if (zc) {
        zc.style.transform = 'none';
        zc.style.zoom = String(cz / 100);
    }
    if (zs) zs.value = cz;
    if (zl) zl.innerText = cz + '%';
    requestAnimationFrame(function() {
        restoreMainDataCenter(centerRatio);
        scheduleConnectorRedraw();
        scheduleLDTriangleRedraw();
    });
}
```

Remove the obsolete transform transition and add `zoom: 1 !important` to print CSS.

For SVG export, measure at intrinsic 100% without changing the user's stored `cz` value:

```javascript
function withIntrinsicReportZoom(callback) {
    if (!zc) return callback();
    var previousZoom = zc.style.zoom;
    zc.style.zoom = '1';
    try {
        return callback();
    } finally {
        zc.style.zoom = previousZoom;
    }
}
```

Wrap `prepareReportExportVisuals()`, SVG collection, sizing, cloning, and serialization inside `withIntrinsicReportZoom(...)`. Keep `beforeprint`/`afterprint` handlers that save the interactive zoom, set `zc.style.zoom='1'` before printing, and call `applyZoom()` after printing.

- [ ] **Step 4: Make Fit independent of the current zoom**

```javascript
function getIntrinsicReportWidth() {
    var candidates = document.querySelectorAll(
        '.gene-gwas-panel, #gene-structure-svg, .data-table, .ld-right-panel'
    );
    var width = 0;
    candidates.forEach(function(node) { width = Math.max(width, node.offsetWidth || 0); });
    return Math.max(width, zc ? zc.scrollWidth : 0, 1);
}

function fitToWindow() {
    var wrapper = document.querySelector('.content-wrapper');
    if (!wrapper) return;
    cz = Math.max(20, Math.min(150, Math.floor(
        wrapper.clientWidth / getIntrinsicReportWidth() * 100
    )));
    applyZoom();
}
```

- [ ] **Step 5: Add range-centering helpers**

Inject `regionStart` and `regionEnd` JS constants and center `.main-data-section` on the current coordinate range after initial load and Apply.

```javascript
function scrollToAppliedRange() {
    var section = document.querySelector('.main-data-section');
    if (!section) return;
    var start = currentFilter.rangeStart === null ? regionStart : currentFilter.rangeStart;
    var end = currentFilter.rangeEnd === null ? regionEnd : currentFilter.rangeEnd;
    var center = (start + end) / 2;
    var ratio = (center - regionStart) / Math.max(1, regionEnd - regionStart);
    var logicalX = gwasLeftMargin + ratio * geneAreaWidth;
    section.scrollLeft = Math.max(0, logicalX - section.clientWidth / 2);
}

function scheduleAppliedRangeCentering() {
    requestAnimationFrame(scrollToAppliedRange);
}
```

- [ ] **Step 6: Run focused and full tests**

```powershell
python -m py_compile haplotype_phenotype_analysis.py test_star_gene_data.py star_gene_validation.py
python -m unittest test_star_gene_data -v
```

Expected: compilation exits 0 and all tests pass.

- [ ] **Step 7: Commit zoom behavior**

```powershell
git add -- haplotype_phenotype_analysis.py test_star_gene_data.py
git commit -m "Keep report viewport stable during zoom"
```

### Task 4: Regenerate VRN-B1 and verify production behavior

**Files:**
- Modify: `star_gene_validation_record.md`
- Regenerate: `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`

- [ ] **Step 1: Rerun the robust target**

```powershell
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Expected: exit 0; generated HTML contains initial range `5447–12446`, CSS zoom, and explicit Apply Range.

- [ ] **Step 2: Run Edge and Firefox interaction verification**

Verify in each browser:

```text
initial currentFilter == pendingRangeFilter == {start: 5447, end: 12446}
visible sequence columns <= 25
Apply button disabled on load
50% Zoom leaves currentFilter unchanged
50% Zoom keeps zoomContent visual width approximately equal to content-wrapper width
Start=90, End=80 previews and applies 80–80
Clear causes no apply; Clear + Apply switches to full range
Reset returns 5447–12446
no pageerror or console error
```

- [ ] **Step 3: Confirm scientific outputs did not change**

Read `validation_summary.json` and the marker database directly. Expected: top `Hap6`, total `1.2016`, n=3, and the 837 bp inserted allele remains present at marker index 251.

- [ ] **Step 4: Update the validation record**

Append target, literature variant, data source, robust mode, command, output directory, unchanged top/match/sample limitation, initial range, zoom browser evidence, and any blocked reason. State explicitly that the initial window uses only positions and gene coordinates.

- [ ] **Step 5: Run final checks, review, commit, and push once**

```powershell
git diff --check -- haplotype_phenotype_analysis.py test_star_gene_data.py star_gene_validation_record.md
git add -- haplotype_phenotype_analysis.py test_star_gene_data.py star_gene_validation_record.md
git diff --cached --check
git commit -m "Verify bounded initial range and zoom layout"
git push origin work/star-gene-validation
```

If GitHub resets the connection, report the local commit and stop retrying.
