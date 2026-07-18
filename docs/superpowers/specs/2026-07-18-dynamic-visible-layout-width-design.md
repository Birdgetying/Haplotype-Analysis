# Dynamic visible-range report width design

## Problem

The integrated report correctly changes the coordinate domain after a display
range is applied, but the gene-structure and GWAS canvases keep the width that
was calculated from every variant in the report. In the current VRN-B1 report,
25 applied-range variants are drawn in a gene area sized for 235 variants
(`4700 px`), producing a `5370 px` SVG inside an approximately `1221 px`
viewport.

The coordinate domain and the canvas width therefore describe different
states: the domain represents the applied range, while the width represents
the unfiltered report.

## Desired behavior

1. The horizontal layout width follows the variants that are actually visible
   after the applied coordinate range and the existing variant filters.
2. Pending edits to numeric range inputs or slider handles do not resize or
   redraw the report. Width changes only when the user clicks **Apply** or when
   another existing action explicitly applies filters.
3. The GWAS plot, gene-structure SVG, sequence table, variant connectors, and
   LD panel remain synchronized after every applied filter change.
4. A small visible set should fit comfortably in the report viewport. A large
   visible set may retain horizontal scrolling so sequence columns remain
   readable.
5. Width selection uses presentation state only. It must not depend on
   phenotype values, haplotype scores, literature markers, validation markers,
   or post-hoc audit outcomes.

## Width model

The visible gene-area width is calculated from the number of currently visible
variant columns:

```text
gene_area_width = max(320 px, visible_variant_count * 20 px)
svg_total_width = fixed_left_columns + gene_area_width + legend_area
```

- `20 px` preserves the existing sequence-column width.
- `320 px` keeps the gene model and axes readable when zero or only a few
  variants remain.
- No maximum is imposed: many visible variants continue to use the existing
  horizontal scroll container rather than compressing bases into unreadable
  columns.

For the current VRN-B1 initial range, 25 visible variants produce a `500 px`
gene area and an approximately `1170 px` total SVG instead of `5370 px`.

## Runtime data flow

`applyFilters()` remains the single applied-state redraw boundary:

1. Filter `gwasData` by the committed range and existing MAF, missingness,
   annotation, variant-type, synonymous, and manual-filter rules.
2. Determine the retained variant positions and table-column indices.
3. Recalculate the shared layout dimensions from the retained position count.
4. Apply the new width to the GWAS panel, GWAS SVG, gene-structure SVG, and
   shared width globals/data attributes.
5. Redraw the gene structure and GWAS plot using the committed coordinate
   domain and updated dimensions.
6. Collapse hidden table columns, then redraw connectors and the LD triangle
   after browser layout settles.

This ordering prevents the plots from using stale dimensions and preserves the
existing explicit-Apply performance rule.

## Edge cases

- **Zero visible variants:** use the `320 px` minimum gene area, keep the gene
  model/axis visible, render an empty GWAS layer, and hide sequence connectors.
- **One visible coordinate:** use the existing single-coordinate D3 domain and
  center it in the minimum-width gene area.
- **Many visible variants:** width grows by `20 px` per variant and remains
  horizontally scrollable.
- **Manual blacklist:** blacklisted variants are excluded from both the visible
  width count and plotted/connected visible-site set.
- **Reset:** restoring the generated initial range recalculates the same compact
  width deterministically.
- **Export/print:** dynamic SVG width attributes remain authoritative so the
  existing intrinsic-size export and print fallbacks measure the applied state.

## Verification

1. Add a Node-backed regression test for the width calculation, including 0,
   1, 25, and large visible-site counts.
2. Add source-level ordering assertions that layout resizing occurs after
   retained positions are known and before gene/GWAS redraw.
3. Run `py_compile` and the complete `test_star_gene_data` suite.
4. Regenerate the real VRN-B1 robust-discovery HTML and update
   `star_gene_validation_record.md` with the command, output path, layout-only
   change, and unchanged scientific result.
5. In Chromium, verify the initial range is still `12,379-12,466`, exactly 25
   variant columns are visible, the total SVG is approximately `1170 px`, the
   main pane no longer has unnecessary multi-thousand-pixel horizontal width,
   Apply/Reset resize only after application, and GWAS/gene/sequence positions
   remain connected.

## Out of scope

- Changing discovery scoring or haplotype ranking.
- Changing automatic range selection.
- Compressing the 20 px sequence columns.
- Reworking report colors, sidebar organization, or score panels.
