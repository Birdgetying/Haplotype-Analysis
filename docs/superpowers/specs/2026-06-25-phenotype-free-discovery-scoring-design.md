# Phenotype-Free Discovery Scoring Design

**Goal:** Make `robust_discovery` useful for prediction and discovery by removing validation phenotype leakage from discovery scores while retaining phenotype and literature checks as post hoc validation.

**Scope:** This change affects `HaplotypeScorer` in `haplotype_phenotype_analysis.py`, star-gene validation tests, and the validation record. It does not call PostGWAS at runtime. Instead, it implements a local PostGWAS-like site weighting model that can consume existing external evidence fields when present.

**Discovery Inputs:** variant annotation severity, promoter/gene boundary distance, variant type, MAF stability, missingness, LD redundancy, non-phenotype sample covariates, and externally supplied GWAS/eQTL-like evidence. Literature causal variants and the current validation phenotype are excluded from discovery scoring.

**Validation Inputs:** phenotype means, PVE, expected trait direction, and literature marker matches remain available only in summary/audit/report layers. These fields may explain whether a ranked haplotype is biologically consistent, but they must not change the discovery score or the score axis.

**Implementation Notes:** `robust_discovery` removes phenotype-derived components (`eb_effect`, `effect_size`, local phenotype p/effect, direction-adjusted total) from ranking. Functional/core position selection uses local annotation and external-evidence weights only. Existing direction-aware summary fields are kept as validation helpers and should be labeled as audit outputs, not discovery outputs.
