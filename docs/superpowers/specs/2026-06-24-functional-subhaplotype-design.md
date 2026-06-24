# Functional Sub-Haplotype Robust Scoring Design

## Goal

Make robust star-gene validation produce positive, biologically interpretable
evidence across TaGW2-B1 full-region, VRN-D1-Kiss2014, and Rht-Zanke2014
without using literature variants as scoring inputs or adding target-specific
rules.

## Problem

Full-region haplotypes can be split by neutral background SNPs. TaGW2-B1 now
has 71 regional SNPs and 112 full haplotypes, so the published promoter
haplotype is covered but no longer appears as the exact top full-region
haplotype. VRN-D1 and Rht-Zanke2014 still validate, which means the next
change must solve over-splitting while preserving the current positive controls.

## Approach

Add a generic `functional_haplotype_groups` layer under `robust_discovery`.
It groups full haplotypes by discovery-selected functional positions rather
than all regional SNPs. Candidate positions are selected from local evidence:
annotation/function weight, promoter/body context, phenotype marker signal,
EB marker effect, MAF stability, missingness, and LD pruning. Literature
markers are not read by this scorer.

Validation summary and literature audit will report raw top, direction-aware
top, core group, and functional group fields separately. A target can be
called positive when a stable discovery layer has significant score-phenotype
association and its post hoc audit matches the literature allele/haplotype.

## Success Criteria

- TaGW2-B1 full-region robust output reports a functional group that covers
  the Qin2014 promoter haplotype in post hoc audit.
- VRN-D1-Kiss2014 remains matched in top/directional validation.
- Rht-Zanke2014 remains matched in direction-aware validation.
- Tests cover background-SNP over-splitting, support filtering, and summary
  field extraction.
- `star_gene_validation_record.md` is updated with commands, outputs, top
  groups, and limitations.
