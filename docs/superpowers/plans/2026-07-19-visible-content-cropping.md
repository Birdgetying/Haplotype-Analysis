# Visible Content Cropping Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the empty initial viewport and multi-thousand-pixel blank horizontal area from the integrated report while preserving explicit Apply behavior and all coordinate-aligned views.

**Architecture:** Keep every sequence column and allele cell in the DOM so range changes remain reversible. At the existing `applyFilters()` commit boundary, set the table to the exact width implied by the retained sequence columns, resize the GWAS/gene canvases with the existing visible-count model, and center scrolling from real pixel offsets instead of multiplying an SVG ratio by the full table scroll width.

**Tech Stack:** Python 3, generated HTML/CSS/JavaScript, Node-backed JavaScript unit harness, `unittest`, D3, Chromium via `agent-browser`.

---

## File map

- Modify `haplotype_phenotype_analysis.py`: add visible-table width helpers, call them from `updateTableColumns()`, and replace proportional auto-centering with real-pixel centering.
- Modify `test_star_gene_data.py`: add RED/GREEN regression coverage for 0/2/23/25/large visible sets, exact table sizing, and scroll behavior.
- Modify `star_gene_validation_record.md`: append the required real VRN-B1 rerun and browser-layout audit without changing its negative scientific conclusion.
- Regenerate `D:\Desktop\project1\star_gene_results\wheat_nature_2024\VRN-B1-fullSequence-IJMS2021__robust_discovery\VRN-B1-fullSequence-IJMS2021.html`: deliver the corrected report artifact outside Git.

### Task 1: Add failing layout regression tests

**Files:**
- Modify: `test_star_gene_data.py:2730-3130`
- Test: `test_star_gene_data.py`

- [ ] **Step 1: Extract the new helpers into the existing Node harness**

Add `calculateVisibleTableWidth`, `applyVisibleTableWidth`, `scrollToAppliedRange`, and the existing range/layout helpers to the function list used by the Node-backed integrated-report test.

- [ ] **Step 2: Add exact-width RED assertions**

Extend the fake DOM with a `.data-table`, `.main-data-section`, and `gene-gwas-panel-workbench`, then assert:

```javascript
assertEqual(calculateVisibleTableWidth(0), 510, 'zero-site table width');
assertEqual(calculateVisibleTableWidth(2), 550, 'two-site table width');
assertEqual(calculateVisibleTableWidth(23), 970, 'twenty-three-site table width');
assertEqual(calculateVisibleTableWidth(25), 1010, 'twenty-five-site table width');
assertEqual(calculateVisibleTableWidth(100), 2510, 'large table remains scrollable');

applyVisibleTableWidth(25);
assertEqual(elements.dataTable.style, {
    width: '1010px', minWidth: '1010px', maxWidth: '1010px'
}, 'table uses exact visible width');
```

- [ ] **Step 3: Add real-pixel scroll RED assertions**

Use a fake section with `clientWidth=1200`, `scrollWidth=1200`, and a panel whose content center is inside the viewport:

```javascript
elements.mainDataSection.scrollLeft = 900;
scrollToAppliedRange();
assertEqual(elements.mainDataSection.scrollLeft, 0, 'narrow content resets blank scroll');

elements.mainDataSection.clientWidth = 800;
elements.mainDataSection.scrollWidth = 2510;
elements.genePanel.offsetLeft = 20;
geneAreaWidth = 2000;
svgTotalWidth = 2670;
scrollToAppliedRange();
assertEqual(
    elements.mainDataSection.scrollLeft,
    Math.min(1710, 20 + 450 + 1000 - 400),
    'wide content centers from the panel pixel origin'
);
```

Also add source-order assertions that `applyVisibleTableWidth(keepPositions.length)` occurs inside `updateTableColumns()` before the forced table reflow and before connector/LD scheduling.

- [ ] **Step 4: Run the targeted test and confirm RED**

Run:

```powershell
python -m unittest test_star_gene_data.StarGeneDataTests.test_display_range_slider_executes_overlap_and_timer_logic
```

Expected: FAIL because `calculateVisibleTableWidth` and `applyVisibleTableWidth` do not exist and `scrollToAppliedRange()` still uses proportional `scrollWidth` mapping.

- [ ] **Step 5: Commit the RED test**

```powershell
git add -- test_star_gene_data.py
git commit -m "test: reproduce integrated report blank scrolling"
```

### Task 2: Implement exact cropping and real-pixel centering

**Files:**
- Modify: `haplotype_phenotype_analysis.py:12020-12040`
- Modify: `haplotype_phenotype_analysis.py:12603-12780`
- Test: `test_star_gene_data.py`

- [ ] **Step 1: Add table-width constants and pure calculation helper**

Near the existing visible-layout constants, define:

```javascript
var haplotypeColumnWidth = 90;
var effectColumnWidth = 180;
var phenotypeColumnWidth = 180;
var sampleCountColumnWidth = 60;

function calculateVisibleTableWidth(visibleCount) {
    var count = Math.max(0, parseInt(visibleCount, 10) || 0);
    return haplotypeColumnWidth + effectColumnWidth + phenotypeColumnWidth +
        count * sequenceColumnWidth + sampleCountColumnWidth;
}
```

- [ ] **Step 2: Add the DOM width adapter**

Implement:

```javascript
function applyVisibleTableWidth(visibleCount) {
    var table = document.querySelector('.data-table');
    var width = calculateVisibleTableWidth(visibleCount);
    if (table) {
        table.style.width = width + 'px';
        table.style.minWidth = width + 'px';
        table.style.maxWidth = width + 'px';
    }
    return width;
}
```

Call `applyVisibleTableWidth(keepPositions.length)` at the end of `updateTableColumns()` after `th/td/col` visibility updates and before `void table.offsetHeight`.

- [ ] **Step 3: Replace proportional scrolling**

Replace the current `contentRatio * section.scrollWidth` calculation with:

```javascript
function scrollToAppliedRange() {
    var section = document.querySelector('.main-data-section');
    var panel = document.getElementById('gene-gwas-panel-workbench');
    if (!section || !panel) return;
    var maxScroll = Math.max(0, section.scrollWidth - section.clientWidth);
    if (maxScroll <= 0) {
        section.scrollLeft = 0;
        return;
    }
    var panelOffset = panel.offsetLeft || 0;
    var geneCenter = panelOffset + gwasLeftMargin + geneAreaWidth / 2;
    section.scrollLeft = Math.max(
        0,
        Math.min(maxScroll, geneCenter - section.clientWidth / 2)
    );
}
```

Do not call this function from pending input/slider handlers. Keep the current double-`requestAnimationFrame` scheduling after initial load, Apply, Reset, and external committed filters.

- [ ] **Step 4: Run targeted tests and confirm GREEN**

Run:

```powershell
python -m unittest test_star_gene_data.StarGeneDataTests.test_display_range_slider_executes_overlap_and_timer_logic
python -m unittest test_star_gene_data.StarGeneDataTests.test_integrated_report_auto_center_scrolls_content_wrapper
```

Expected: both PASS.

- [ ] **Step 5: Run syntax and complete unit verification**

Run:

```powershell
python -m py_compile haplotype_phenotype_analysis.py test_star_gene_data.py
python -m unittest test_star_gene_data
git diff --check
```

Expected: compilation succeeds, all tests pass with zero failures, and `git diff --check` prints no errors.

- [ ] **Step 6: Commit implementation**

```powershell
git add -- haplotype_phenotype_analysis.py test_star_gene_data.py
git commit -m "Fix integrated report visible content cropping"
```

### Task 3: Regenerate and verify the real VRN-B1 report

**Files:**
- Modify: `star_gene_validation_record.md`
- Regenerate: `D:\Desktop\project1\star_gene_results\wheat_nature_2024\VRN-B1-fullSequence-IJMS2021__robust_discovery\VRN-B1-fullSequence-IJMS2021.html`

- [ ] **Step 1: Rerun the real analysis**

Run:

```powershell
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery --database-root D:\Desktop\project1\star_gene_database --results-root D:\Desktop\project1\star_gene_results --manifest D:\Desktop\project1\star_gene_manifest.json
```

Expected: exit code 0 and updated integrated HTML. The scientific result must remain the corrected negative validation: Hap11 raw top, Hap4 displayed top, Hap5 anchor top, and Hap6 insertion carrier at rank 15.

- [ ] **Step 2: Verify initial browser state**

Open the generated HTML in a named Chromium session at 1440×900. Assert:

```text
range = 12379-12466
visible sequence columns = 25
table width = 1010 px
gene/GWAS panel width = 1170 px
main-data-section.scrollLeft = 0
main-data-section.scrollWidth <= main-data-section.clientWidth (at this viewport)
```

Capture a screenshot showing data in the initial viewport rather than a blank page.

- [ ] **Step 3: Verify Apply behavior at 23 and 2 sites**

Edit numeric inputs without Apply and confirm widths/scroll position do not change. Click Apply for `12381-12465` and assert 23 columns, 970 px table width, 1130 px gene/GWAS width, and aligned connector positions. Apply `13070-13080` and assert two columns, 550 px table width, 990 px gene/GWAS width, `scrollLeft=0`, and one visible sequence/gene/GWAS/connector representation for position 13077.

- [ ] **Step 4: Verify browser health and close the session**

Confirm no page errors and no console error/warning entries. Check LD redraw and Reset, capture the final screenshot, then close the named browser session.

- [ ] **Step 5: Append the required validation record**

Add a dated entry to `star_gene_validation_record.md` containing the target, literature 837 bp insertion, data source, `robust_discovery`, exact run command, output directory, raw/displayed/anchor tops, negative match, n limitations, layout metrics for 25/23/2 sites, browser diagnostics, and any blocked reason. Explicitly state that this is a layout-only rerun and that the literature marker remains post-hoc only.

- [ ] **Step 6: Run final verification, commit, and push**

Run:

```powershell
python -m py_compile haplotype_phenotype_analysis.py test_star_gene_data.py
python -m unittest test_star_gene_data
python run_rice_test.py
git diff --check
git status --short --branch
```

If `run_rice_test.py` exceeds the local timeout without a traceback, record it as timeout rather than pass/fail. Commit only the record change after verification:

```powershell
git add -- star_gene_validation_record.md
git commit -m "Record visible content cropping verification"
git push -u origin codex/dynamic-visible-width
```
