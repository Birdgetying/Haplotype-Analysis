# Dynamic Visible Report Width Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Resize the integrated report's GWAS and gene-structure canvases from the currently applied visible-site set so compact ranges no longer retain the full-report width.

**Architecture:** Keep `applyFilters()` as the only applied-state redraw boundary. Add one pure JavaScript width calculator and one DOM adapter; after range, quality, annotation, type, synonymous, and manual filters determine the unique visible positions, update the shared width globals and DOM dimensions before redrawing the gene model, GWAS plot, connectors, and LD triangle.

**Tech Stack:** Python report generator, embedded JavaScript, D3.js, Node.js regression checks, Python `unittest`, Chromium through `agent-browser`.

---

## File map

- Modify `haplotype_phenotype_analysis.py`: embedded integrated-report JavaScript width model, DOM resizing, and `applyFilters()` redraw ordering.
- Modify `test_star_gene_data.py`: Node-backed width behavior tests and source-order regression assertions.
- Modify `star_gene_validation_record.md`: append the real VRN-B1 rerun, layout measurements, unchanged scientific result, and limitations.
- Regenerate `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`: runtime artifact used for browser verification; do not stage it unless it is tracked and intentionally part of the existing repository policy.

### Task 1: Add failing width-model and redraw-order regressions

**Files:**
- Modify: `test_star_gene_data.py:2336-2573`
- Test: `test_star_gene_data.py`

- [ ] **Step 1: Add a source-order regression test**

Add this test near the existing applied-range report tests:

```python
def test_integrated_report_resizes_shared_canvases_from_applied_visible_sites(self):
    source = Path("haplotype_phenotype_analysis.py").read_text(encoding="utf-8")
    integrated_start = source.index("def generate_integrated_html")
    integrated_end = source.index("def generate_haplotype_network_html", integrated_start)
    integrated_block = source[integrated_start:integrated_end]

    self.assertIn("function calculateVisibleLayoutDimensions(", integrated_block)
    self.assertIn("function applyVisibleLayoutDimensions(", integrated_block)
    apply_start = integrated_block.index("function applyFilters()")
    apply_end = integrated_block.index("function updateTableColumns", apply_start)
    apply_block = integrated_block[apply_start:apply_end]
    layout_idx = apply_block.index("applyVisibleLayoutDimensions(visiblePositions.length)")
    gene_idx = apply_block.index("updateGeneStructureDomain()")
    gwas_idx = apply_block.index("drawGWASPlot(visibleGwasData, coordinateDomain)")
    self.assertLess(layout_idx, gene_idx)
    self.assertLess(layout_idx, gwas_idx)
    self.assertIn("new Set(varPositions)", apply_block)
    self.assertIn("manualBlacklist", apply_block)
```

- [ ] **Step 2: Extend the Node-backed JavaScript test with expected widths**

Include `calculateVisibleLayoutDimensions` in the extracted function list and add:

```javascript
assertEqual(
    calculateVisibleLayoutDimensions(0),
    {geneAreaWidth: 320, svgTotalWidth: 990, gwasPlotWidth: 540},
    'zero sites use readable minimum'
);
assertEqual(calculateVisibleLayoutDimensions(1).geneAreaWidth, 320, 'one site uses minimum');
assertEqual(
    calculateVisibleLayoutDimensions(25),
    {geneAreaWidth: 500, svgTotalWidth: 1170, gwasPlotWidth: 720},
    'twenty-five sites fit the report viewport'
);
assertEqual(calculateVisibleLayoutDimensions(100).geneAreaWidth, 2000, 'large sets remain scrollable');
```

Define the test harness constants before these assertions:

```javascript
var gwasLeftMargin = 450;
var sequenceColumnWidth = 20;
var minimumGeneAreaWidth = 320;
var legendAreaWidth = 220;
```

- [ ] **Step 3: Run the targeted tests and verify RED**

Run:

```powershell
python -m unittest test_star_gene_data.StarGeneDataTests.test_integrated_report_resizes_shared_canvases_from_applied_visible_sites test_star_gene_data.StarGeneDataTests.test_display_range_slider_executes_overlap_and_timer_logic
```

Expected: FAIL because `calculateVisibleLayoutDimensions` and `applyVisibleLayoutDimensions` do not exist.

- [ ] **Step 4: Commit the failing tests**

```powershell
git add -- test_star_gene_data.py
git diff --cached --check
git commit -m "Test dynamic visible report widths"
```

### Task 2: Implement applied-visible-site layout resizing

**Files:**
- Modify: `haplotype_phenotype_analysis.py:11144-11153`
- Modify: `haplotype_phenotype_analysis.py:11648-11886`
- Modify: `haplotype_phenotype_analysis.py:12445-12534`
- Test: `test_star_gene_data.py`

- [ ] **Step 1: Add explicit shared layout constants**

Replace the implicit fixed-width assumptions beside the existing width globals with:

```javascript
var sequenceColumnWidth = 20;
var minimumGeneAreaWidth = 320;
var legendAreaWidth = 220;
var svgTotalWidth = {svg_total_width};
var geneAreaWidth = {gene_area_width_js};
var gwasPlotWidth = {gwas_plot_width};
var gwasLeftMargin = {gwas_left_margin};
```

Do not change the fixed left-column widths or the 20 px sequence-column width.

- [ ] **Step 2: Add the pure width calculator**

Place this function before coordinate mapping functions:

```javascript
function calculateVisibleLayoutDimensions(visibleVariantCount) {
    var count = Math.max(0, Math.floor(Number(visibleVariantCount) || 0));
    var nextGeneAreaWidth = Math.max(minimumGeneAreaWidth, count * sequenceColumnWidth);
    return {
        geneAreaWidth: nextGeneAreaWidth,
        svgTotalWidth: gwasLeftMargin + nextGeneAreaWidth + legendAreaWidth,
        gwasPlotWidth: nextGeneAreaWidth + legendAreaWidth
    };
}
```

- [ ] **Step 3: Add the DOM adapter**

Add one function that mutates the shared runtime dimensions and the existing panel/SVG elements:

```javascript
function applyVisibleLayoutDimensions(visibleVariantCount) {
    var dimensions = calculateVisibleLayoutDimensions(visibleVariantCount);
    geneAreaWidth = dimensions.geneAreaWidth;
    svgTotalWidth = dimensions.svgTotalWidth;
    gwasPlotWidth = dimensions.gwasPlotWidth;

    var panel = document.getElementById('gene-gwas-panel-workbench');
    if (panel) {
        panel.style.width = svgTotalWidth + 'px';
        panel.style.minWidth = svgTotalWidth + 'px';
    }
    var gwasContainer = document.getElementById('gwas-gene-viz');
    if (gwasContainer) gwasContainer.style.width = svgTotalWidth + 'px';

    var geneSvg = document.getElementById('gene-structure-svg');
    if (geneSvg) {
        geneSvg.setAttribute('width', String(svgTotalWidth));
        geneSvg.setAttribute('data-gene-width', String(geneAreaWidth));
        geneSvg.style.width = svgTotalWidth + 'px';
    }
    return dimensions;
}
```

Do not use viewport width, physical base-pair span, phenotype data, score data, or literature markers in this function.

- [ ] **Step 4: Reorder `applyFilters()` around the final visible-site set**

Keep the current filter predicates, then use the manual-filtered position set to define unique visible positions and GWAS points:

```javascript
var posSet = {};
filtered.forEach(function(d) { posSet[d.pos] = true; });
manualBlacklist.forEach(function(pos) { delete posSet[pos]; });

// Existing table-header scan fills varIndices and varPositions.
var visiblePositions = Array.from(new Set(varPositions));
var visibleGwasData = filtered.filter(function(d) {
    return posSet[d.pos] !== undefined;
});

updateTableColumns(varIndices, varPositions);
applyVisibleLayoutDimensions(visiblePositions.length);
var coordinateDomain = getAppliedCoordinateDomain();
updateGeneStructureDomain();
drawGWASPlot(visibleGwasData, coordinateDomain);
```

Move the existing gene/GWAS redraw from before the position-set calculation to this location. Keep variant visibility, connector scheduling, LD redraw scheduling, and manual-filter visuals after these calls.

- [ ] **Step 5: Run targeted tests and verify GREEN**

Run:

```powershell
python -m unittest test_star_gene_data.StarGeneDataTests.test_integrated_report_resizes_shared_canvases_from_applied_visible_sites test_star_gene_data.StarGeneDataTests.test_display_range_slider_executes_overlap_and_timer_logic test_star_gene_data.StarGeneDataTests.test_applied_display_range_drives_gene_structure_and_gwas_domains
```

Expected: `Ran 3 tests` and `OK`.

- [ ] **Step 6: Run the complete unit suite**

```powershell
python -m py_compile haplotype_phenotype_analysis.py test_star_gene_data.py
python -m unittest test_star_gene_data
git diff --check
```

Expected: compile exit 0, at least 155 tests with `OK`, and no whitespace errors.

- [ ] **Step 7: Commit the implementation**

```powershell
git add -- haplotype_phenotype_analysis.py test_star_gene_data.py
git diff --cached --check
git commit -m "Resize reports from visible variants"
```

### Task 3: Regenerate and verify the real VRN-B1 report

**Files:**
- Modify: `star_gene_validation_record.md`
- Regenerate: `D:\Desktop\project1\star_gene_results\wheat_nature_2024\VRN-B1-fullSequence-IJMS2021__robust_discovery\VRN-B1-fullSequence-IJMS2021.html`

- [ ] **Step 1: Run the real target from the feature worktree**

```powershell
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery --database-root D:\Desktop\project1\star_gene_database --results-root D:\Desktop\project1\star_gene_results --manifest D:\Desktop\project1\star_gene_manifest.json
```

Expected: exit 0 and the integrated HTML is regenerated at the path above.

- [ ] **Step 2: Verify embedded scientific and range data before opening the browser**

Use PowerShell `Get-Content -Raw` plus regular expressions/`ConvertFrom-Json` to verify:

```text
initialDisplayRange = 12379-12466
visible range positions = 25
score_mode = robust_discovery
top total haplotype = Hap6
Hap6 total = 1.2016
Hap6 sample count = 3
literature marker 13077 is outside the neutral initial window
```

If the scientific values differ, stop and investigate rather than recording an unchanged result.

- [ ] **Step 3: Run Chromium layout verification**

Open the file in a named browser session:

```powershell
npx --yes agent-browser --session dynamic-width-check --allow-file-access open "file:///D:/Desktop/project1/star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html"
npx --yes agent-browser --session dynamic-width-check wait --load networkidle
```

Evaluate DOM metrics and require:

```text
currentFilter = pendingRangeFilter = 12379-12466
visible unique variant positions = 25
geneAreaWidth = 500
svgTotalWidth = 1170
GWAS SVG width = gene-structure SVG width = 1170
main-data scroll width is no longer the old 5390 px
pageerror count = 0
console error count = 0
```

Then inspect browser diagnostics explicitly:

```powershell
npx --yes agent-browser --session dynamic-width-check console
npx --yes agent-browser --session dynamic-width-check errors
```

Expected: neither command reports an error entry.

- [ ] **Step 4: Verify explicit-Apply behavior**

Record the current SVG width, edit the numeric range inputs to `12381-12465`, and confirm the width and coordinate axis do not change while the state is pending. Click the real Apply button and verify all of the following change together:

```text
currentFilter
visible unique position count
geneAreaWidth/svgTotalWidth
GWAS x-axis domain
gene-structure x-axis domain
sequence columns/connectors
```

Click Reset and confirm `12379-12466`, 25 positions, `500/1170 px`, and disabled Apply are restored.

- [ ] **Step 5: Capture evidence outside the repository and close the browser**

Write any screenshot to `%TEMP%\VRN-B1-dynamic-visible-width.png`, not the repository, then run:

```powershell
npx --yes agent-browser --session dynamic-width-check close
```

- [ ] **Step 6: Update the required star-gene validation record**

Append a dated entry containing:

```text
target: VRN-B1-fullSequence-IJMS2021
literature variant: Vrn-B1f 837 bp promoter insertion, coordinate 13077 in report audit data
data source: existing IJMS2021 full-sequence database
score mode: robust_discovery
exact run command
output directory
top: Hap6, total 1.2016, n=3
post-hoc literature match: yes
layout: initial 12379-12466, 25 sites, gene area 500 px, total SVG 1170 px
scientific result changed: no
limitations: Hap6 has only 3 samples; range selection and layout do not use the literature marker or phenotype
```

- [ ] **Step 7: Commit the validation record**

```powershell
git add -- star_gene_validation_record.md
git diff --cached --check
git commit -m "Record dynamic VRN-B1 report layout"
```

### Task 4: Final regression, review, and handoff

**Files:**
- Verify: `haplotype_phenotype_analysis.py`
- Verify: `test_star_gene_data.py`
- Verify: `star_gene_validation_record.md`

- [ ] **Step 1: Run fresh final verification**

```powershell
python -m py_compile haplotype_phenotype_analysis.py test_star_gene_data.py
python -m unittest test_star_gene_data
git diff --check 9eb1911..HEAD
git status --short --branch
```

Expected: compile exit 0, all tests `OK`, commit-range diff check clean, and no uncommitted feature files.

- [ ] **Step 2: Attempt the project smoke test with a bounded runtime**

Run `python run_rice_test.py` with a 180-second timeout. Record `PASS` only on exit 0. If the no-`pysam` linear scan times out, record `TIMEOUT / not verified`, terminate only the process started by this command, and confirm no child remains.

- [ ] **Step 3: Request code review**

Review the exact range `9eb1911..HEAD` against the approved design. Require explicit checks for width calculation, duplicated/hidden DOM columns, manual blacklist behavior, Apply-only redraw, export dimensions, scientific-data leakage, tests, and changed-file scope. Fix every Critical or Important issue and rerun Step 1.

- [ ] **Step 4: Push the feature branch once**

```powershell
git push -u origin codex/dynamic-visible-width
```

If GitHub connectivity fails, report the exact error and do not loop indefinitely.

- [ ] **Step 5: Use the finishing workflow**

Use `finishing-a-development-branch` to offer local integration, pull request, keep-as-is, or discard. Do not remove the worktree or delete the feature branch without explicit confirmation because repository instructions prohibit unconfirmed deletion.
