# Star Gene Validation Record

This file is the running record for star-gene validation of haplotype scoring.
Update it whenever a new star-gene target, score mode, dataset, or validation
result is tested.

## Update Rules

- Literature variants and literature haplotypes are used only for offline
  validation, not as inputs to discovery scoring.
- Record both the biological positive-control result and the score-rank result.
- Always record the score mode, command, output directory, top-scored haplotype,
  whether the top haplotype matches the literature allele/haplotype, and any
  data limitation.
- Keep blocked targets in this file with the reason, so they are not repeatedly
  treated as missing work.
- If an HTML/JSON output is ignored by git, record the local path and key
  numbers here instead of forcing large generated files into git.

## Current Summary

Current proof set for the teacher's "at least three" requirement:

1. `TaGW2-B1-remoteSNP`: strong regional proof after the 2026-06-23
   full-region rerun. The complete chr6B VCF gives a strong TaGW2-B1 TGW
   association (`P=1.98e-06`, PVE `27.17%`), and the Qin2014 `A/G/C`
   promoter haplotype is covered and segregating. However, full-region
   background SNPs split that promoter class, so the exact literature
   promoter haplotype is no longer the raw top full-region haplotype.
2. `VRN-D1-Kiss2014__robust_discovery`: stable diagnostic-marker proof.
   The top robust single-marker haplotype is `VRN-D1=1` (`Hap2`, n=38)
   and has significant DEV49/DEV59 heading-date score regressions.
3. `Rht-Zanke2014__robust_discovery`: marker-panel proof after direction-aware scoring.
   The retained raw `total` still explains the old high-plant-height Hap2
   artifact, but the active `directional_total` score ranks `Hap1=B1a|D1b`
   first (`n=214`), and the literature audit marks `Rht-D1b` as both
   `matched_top_haplotype` and `matched_directional_top_haplotype`.
   This is no longer treated as a strict single-gene/base-variant Rht proof
   because it combines `Rht-B1` and `Rht-D1` marker states.
4. `VRN-B1-fullSequence-IJMS2021__robust_discovery`: exact full-sequence
   proof after the 2026-06-29 sequence-order fix and 2026-06-30 rerun.
   The raw robust top and anchor top are both `Hap6`, and the rendered HTML
   shows `Hap6` carries the literature `Vrn-B1f` 837 bp insertion at
   `13077` (`+837bp`). Reliability caveat: n=3.

| Target | Trait | Score mode | Status | Key result | Conclusion |
|---|---|---|---|---|---|
| TaGW2-B1 / TaGW2-6B | TGW | default | tested | Full-region rerun uses 71 SNPs; top is rare full-region Hap43, not the exact Qin2014 promoter haplotype | Biological signal present, but full-region background SNPs split the literature promoter haplotype into many sub-haplotypes |
| TaGW2-B1 / TaGW2-6B | TGW | robust_discovery | tested | Full-region rerun uses 71 SNPs; top is Hap2, while Qin2014 `A/G/C` is present across many exact full-region haplotypes | Full-region discovery view is more realistic, but exact promoter-haplotype validation should be read from literature marker audit/grouping, not raw full-region top rank |
| TaGW2-A1 | TGW | default | tested | Gene-level signal exists, but `SNP-494 A` is absent in current complete Watkins samples | Not an exact Hap5 validation |
| TaGW2-D1 | TGW | default | tested | Weak low-confidence signal with only 2 retained SNPs | Limited marker coverage; not strong proof |
| VRN-Kiss2014 | heading date | default | tested | Literature all-spring `1\|1\|1` was top-scored but n=2 | Useful audit check, weak statistical proof |
| VRN-Kiss2014 | heading date | robust_discovery | tested | Top score became `Hap2=0\|0\|1`, n=31; all-spring `1\|1\|1` was no longer top | Robust mode favors stable VRN-D1 signal over rare all-spring haplotype |
| VRN-D1-Kiss2014 | heading date | robust_discovery | tested | Single-marker `VRN-D1=1` top haplotype `Hap2`, n=38; DEV49 P `2.98e-07`, DEV59 P `7.49e-05` | Counts as stable marker-level proof |
| Rht-D1b | plant height | default | tested | Top-scored haplotype contained functional T allele but only 5 carriers | Weak positive control |
| Rht-D1b | plant height | robust_discovery | tested | 2026-06-26 strict rerun: exact stop-gained `chr4D:18781242 G>T`; top score is functional T haplotype `Hap2`, n=3, mean `84.583 cm`, score `1.0938`; audit `matched_top_haplotype` | Weak exact single-variant positive control; carrier count too small and support-shrunk directional top is wild-type `Hap1` |
| Rht-Zanke2014 | plant height | robust_discovery + direction-aware score axis | tested | Active top-scored haplotype is `Hap1=B1a\|D1b`, n=214, matching the `Rht-D1b` diagnostic marker; retained raw total explains the old Hap2 artifact | Marker-panel reference only for strict Rht proof, because it combines `Rht-B1` and `Rht-D1` marker states |
| Rht-D1-Zanke2014 | plant height | robust_discovery | tested | Single-marker Rht-D1 target from Zanke2014 Table S2. Directional top is exact `Hap1=D1b`, n=214; raw total top is `Hap3=D1b/Rht-D1a`, n=2 | Single-gene marker-level proof that avoids B1/D1 combination; read direction-aware top, not the n=2 mixed raw top |
| Rht-B1b | plant height | any | blocked | Functional T allele has 0 phenotype-overlap carriers | Cannot validate scoring in current samples |
| TaPIF4 | TGW / plant height | any | blocked | Public repository lacks final per-sample promoter haplotype table and matching phenotype table | Not usable until final data are obtained |

## Literature Decisive Variant Checklist

This table separates the offline literature validation target from the
discovery score. A target supports scoring only when the exact literature
variant or functional haplotype is covered, segregates in the current
sample-phenotype overlap, and is recovered by the top-scored haplotype under
the tested score mode.

| Target | Literature/source | Decisive variant or functional haplotype | Current marker/data | Coverage and segregation | Top-score match | Validation verdict |
|---|---|---|---|---|---|---|
| TaGW2-B1 / TaGW2-6B | Qin et al. 2014, BMC Plant Biology, `https://doi.org/10.1186/1471-2229-14-107` | Favorable promoter haplotype `Hap-6B-1`, diagnostic states `-1709/-721/-83 = A/G/C`; TaGW2 encodes a RING-type E3 ubiquitin ligase and the favored haplotype increases grain width/weight. | Current full-region database uses complete WWWG2B chr6B VCF over `chr6B:291759689-291778752`: 756 Watkins TGW-overlap samples, 71 SNPs. The three diagnostic records are `chr6B:291759689 C>A`, `chr6B:291760677 G>A`, `chr6B:291761315 T>C`. | Covered and segregating, but `A/G/C` is split across many exact full-region haplotypes. | Full-region default/robust: no, strict audit is `present_but_not_top`. Historical 3-marker robust run: yes, exact `A\|G\|C` was top. | Strong regional GW2 proof and useful demonstration that exact literature-marker validation should be separated from full-region discovery ranking. |
| TaGW2-A1 / TaGW2-6A | Jaiswal et al. 2015, PLoS ONE, `https://doi.org/10.1371/journal.pone.0129400` | Promoter `SNP-494` is the causal/expression-regulating marker; favorable allele `A`; four-SNP Hap5 pattern tracked as `SNP-988/SNP-739/SNP-593/SNP-494 = G/A/G/A`. | Local WatSeq chr6A SNP VCF plus remote WheatOmics micro-VCF. | Local data miss `SNP-494`; remote data cover `chr6A:237734341 G>A`, but expected `A` has 0 complete Watkins carriers, so Hap5 is not observed. | No exact validation. Some top haplotypes match the first three promoter SNPs, but not the decisive `SNP-494 A` state. | Gene-level TGW signal only; cannot be used as exact functional-haplotype proof in current samples. |
| TaGW2-D1 | WheatOmics GW2 homolog annotation and current WatSeq A/D run. | No precise published decisive D1 functional marker was configured for this validation; current target used two regional SNPs around `TraesCS6D02G176900`. | WatSeq chr6D SNP VCF; 820 Watkins TGW-overlap samples; 2 retained SNPs. | Region covered sparsely; no literature decisive allele/haplotype available in manifest. | Not applicable. | Low-confidence regional signal only; not a star-gene proof. |
| Rht-B1b | Rht-1 DELLA literature, including Peng et al. 1999 and Ellis et al. 2002 perfect-marker work; WheatOmics annotation used for coordinate. | Canonical semi-dwarf stop-gained DELLA allele, `chr4B:30861571 C>T`, `c.190C>T`, `p.Gln64*`; mechanism is GA-insensitive reduced height through truncated/stabilized DELLA repression of growth. | WheatOmics merged SNP VCF plus Watkins CFLN06 plant-height phenotype. | Marker is present, but expected `T` has 0 phenotype-overlap carriers. | No, blocked before scoring. | Data-blocked in current Watkins panel; not a method-negative result. |
| Rht-D1b | Rht-1 DELLA literature, including Peng et al. 1999 and Ellis et al. 2002 perfect-marker work; WheatOmics annotation used for coordinate. | Canonical semi-dwarf stop-gained DELLA allele, `chr4D:18781242 G>T`, `c.181G>T`, `p.Glu61*`; same GA-insensitive DELLA mechanism as Rht-B1b. | WheatOmics merged SNP VCF plus Watkins CFLN06 plant-height phenotype. | Covered and segregating, but only 5 functional-allele carriers overall; robust top haplotype carrier count is 3. | Yes for raw strict top score: robust top haplotype is exact `T` (`Hap2`, n=3). Direction-aware support shrinkage chooses wild-type `Hap1`, so this remains weak. | Weak exact single-variant positive control: the decisive variant is recovered, but carrier count and R2 are too small for strong proof. |
| Rht-Zanke2014 | Zanke et al. 2014 PLOS ONE Table S2, `https://doi.org/10.1371/journal.pone.0113287`; candidate-gene genotypes use Ellis et al. Rht markers. | Diagnostic `Rht-B1b` and `Rht-D1b` marker states with multi-environment plant-height phenotypes. | Downloaded Table S2 workbook; 368 complete varieties; combined `Rht-B1/Rht-D1` marker panel. | Covered and segregating; `Rht-D1b` occurs in 216 samples and exact `Hap1=B1a\|D1b` has n=214. | Active score axis `directional_total` ranks exact `Rht-D1b` marker combination first; audit status is `matched_top_haplotype` and `matched_directional_top_haplotype`. | Useful marker-panel reference after fixing direction handling, but not strict proof because it combines two Rht genes. |
| Rht-D1-Zanke2014 | Zanke et al. 2014 PLOS ONE Table S2, `https://doi.org/10.1371/journal.pone.0113287`; single-marker extraction from the same candidate-gene table. | Diagnostic `Rht-D1b` marker state only. It is a marker-level proxy for the known Rht-D1b DELLA semi-dwarf allele, not a base-level `chr4D:18781242 G>T` VCF record. | Downloaded Table S2 workbook; 368 complete varieties; only the `Rht-D1` marker column is used. | Covered and segregating; `D1b` occurs in 216 samples, with exact `Hap1=D1b` n=214 and mixed `Hap3=D1b/Rht-D1a` n=2. | Direction-aware top is exact `Hap1=D1b` and audit `directional_validation_status=matched_directional_top_haplotype`. Raw total top is the tiny mixed `Hap3`, so it should not be interpreted as the biological top. | Single-gene marker-level proof that answers the teacher's no-combined-gene concern. It complements, but does not replace, the weak exact-SNP `Rht-D1b` result. |
| VRN-Kiss2014 | Kiss et al. 2014 supplementary marker table; VRN1 biology supported by VRN literature on promoter/intron-1 structural variants. | Diagnostic spring/dominant marker states `VRN-A1=1`, `VRN-B1=1`, `VRN-D1=1`; exact causal events are often promoter or first-intron structural variants, so this is marker-level validation. | Extracted Kiss2014 workbook; 676 samples with VRN marker states and DEV49/DEV59 heading dates. | Covered and segregating; all-three `1\|1\|1` exists but n=2; `VRN-D1=1` single-state haplotype `0\|0\|1` has n=31. | Default: all-three `1\|1\|1` matched top; robust_discovery: top became stable `0\|0\|1`, so all-three is present but not top. | Useful marker-level check; not a strong causal-SV proof. Shows rare-candidate and stable-common ranks should be reported separately. |
| TaPIF4 / TaSG-D1-TaPIF4 | Cao et al. 2024, Nature Communications, `https://doi.org/10.1038/s41467-024-46419-0`; code at `https://github.com/QinZhen1995/CAU-TaSG`. | Primary functional variant in that paper is `TaSG-D1 E286K`, which enhances TaPIF4 phosphorylation/stability under heat stress. TaPIF4 promoter haplotypes include `Del` with 275-bp and 12-bp deletions, and `InDel` with 405-bp deletion plus 1909-bp insertion; these reduce heat-induced TaPIF4 expression. | Current public GitHub exposes scripts and `samlist`; local project does not have final `PIF_hap.txt`, `coverage.martix`, or matching sample-level TGW/height/heat phenotype table. | Not covered in the current database. | No, not run. | Blocked. This target cannot prove scoring until final per-sample promoter haplotypes or raw BAM-derived coverage matrix plus phenotypes are obtained. |

## Detailed Results

### TaGW2-B1 / TaGW2-6B

Literature positive control:
Qin et al. 2014 reported favorable promoter haplotype `Hap-6B-1` with
diagnostic states `-1709=A`, `-721=G`, and `-83=C`. In the current VCF
coordinate system these are `chr6B:291759689 C>A`,
`chr6B:291760677 G>A`, and `chr6B:291761315 T>C`.

Data:
`star_gene_database/wheat_nature_2024/TaGW2-B1-remoteSNP`,
2026-06-23 full-region rerun:
The complete downloaded WWWG2B chr6B VCF was used beyond the three
literature SNPs. The current database at
`star_gene_database/wheat_nature_2024/TaGW2-B1-remoteSNP` now contains
756 Watkins TGW-overlap samples, 71 regional SNPs from
`chr6B:291759689-291778752`, and 112 haplotypes. Variant annotation is
8 promoter SNPs plus 63 intronic SNPs. The integrated HTML displays 23
positions because the main plot only shows positions variable among the
top displayed haplotypes; the database itself has all 71 retained SNPs.

Commands:
`python prepare_wheat2024_tagw2_b1_remote_snp.py --vcf D:\Desktop\data\GW2\chr6B.HARD.SNP.Missing-unphasing.ID.ann.finalSID.1047.allele2_retain.hard_retain.InbreedingCoeff_retain.vcf.gz --min-haplotype-count 1 --max-missing-rate 0.2`

`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP`

`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --score-mode robust_discovery`

Outputs:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP/TaGW2-B1-remoteSNP.html`
and
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery/TaGW2-B1-remoteSNP.html`.

Full-region result:
The literature `A/G/C` promoter pattern is still covered and segregating, but
adding 68 background regional SNPs splits that promoter class into many exact
full-region haplotypes. Strict audit therefore reports the Qin2014 haplotype
as `present_but_not_top` in both default and robust full-region runs. Default
top is rare `Hap43` (n=2, TGW mean `52.544`), and robust top is `Hap2`.
The full-region association remains strong (`P=1.98e-06`, PVE `27.17%`),
but it should be interpreted as discovery-style regional signal, not as an
exact literature-promoter-haplotype recovery.

Earlier 3-marker validation dataset:
Before the full-region rerun, the same target directory was built from only
the three Qin2014 promoter SNPs: 816 Watkins TGW-overlap samples, 3 promoter
SNPs, and 6 haplotypes. The historical records below refer to that 3-marker
positive-control run.

Default command:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP`

Default result:
Top score was `Hap5=C/A|G|C`, n=7, score `4.2581`; exact literature
haplotype was `Hap1=A|G|C`, n=371, score `2.5746`, TGW mean `40.970`.
Strict audit status for the combined literature haplotype was
`contained_in_top_haplotype_not_exact`.

Robust command:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --score-mode robust_discovery`

Robust output:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery/TaGW2-B1-remoteSNP.html`

Robust result:
Top score became exact literature `Hap1=A|G|C`, n=371, score `1.8132`,
sample reliability `0.9488`, ambiguity factor `1.0000`; rare ambiguous
`Hap5=C/A|G|C` dropped to third, n=7, score `0.9733`, reliability `0.2593`,
ambiguity factor `0.7333`. Literature audit status became
`matched_top_haplotype`.

Interpretation:
This is the strongest current positive control for the improved discovery
scoring mode. The scoring formula did not use the literature haplotype label;
the literature comparison was applied only after ranking.

2026-06-17 report UI rerun:
Both score modes were rerun to refresh the HTML reports after adding the
in-page score-mode toggle.

Commands:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP`

`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --score-mode robust_discovery`

Outputs:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP/TaGW2-B1-remoteSNP.html`
and
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery/TaGW2-B1-remoteSNP.html`.

HTML verification:
Each integrated report and each standalone `haplotype_score.html` embeds both
`default` and `robust_discovery` score JSON when the sibling mode result
directory exists. The default report opens in `default`; the robust report opens
in `robust_discovery`. Switching changes only the score plot/statistics/title
and does not use the literature haplotype as a scoring input.

2026-06-20 Post-GWAS evidence panel rerun:
The integrated HTML was refreshed after adding the local Post-GWAS evidence
panel. The report generator now embeds `post_gwas_evidence_json` into the
JavaScript data block, so the panel can render without stopping the downstream
score-mode JavaScript.

Commands:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP`

`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --score-mode robust_discovery`

Outputs:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP/TaGW2-B1-remoteSNP.html`
and
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery/TaGW2-B1-remoteSNP.html`.

Top-scored haplotypes:
Default mode still ranks rare ambiguous `Hap5=C/A|G|C` first, n=7, total
score `4.2581`, so it does not exactly match Qin2014 `Hap-6B-1`.
Robust mode ranks exact literature `Hap1=A|G|C` first, n=371, total score
`1.8132`, sample reliability `0.9488`, matching the functional haplotype.

HTML verification:
Both regenerated integrated reports contain `Local Post-GWAS Evidence`,
`Strict local support`, the top local marker `6B:291,759,689`, and the
Original/Robust score-mode controls. Static HTML inspection confirmed no
unreplaced `{post_gwas_evidence_json}` placeholder remains.

2026-06-20 Post-GWAS evidence layout rerun:
The same default and robust commands were rerun after moving the evidence
panel above the network/GWAS row. This keeps the panel outside the
GWAS-to-gene guide-line corridor, so marker connector lines no longer cross
the evidence cards. Browser geometry check on the robust report confirmed
the DOM order `post-gwas-evidence`, `top-section`, `main-data-section`, with
the evidence panel ending above the GWAS panel.

2026-06-20 GWAS/LD alignment layout rerun:
The same default and robust commands were rerun after tightening the GWAS
guide-line and LD-triangle alignment in the integrated HTML. This was a
layout-only rerun; the discovery scores and literature validation conclusion
were not changed.

Commands:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP`

`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --score-mode robust_discovery`

Outputs:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP/TaGW2-B1-remoteSNP.html`
and
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery/TaGW2-B1-remoteSNP.html`.

Top-scored haplotypes:
Default mode still ranks rare ambiguous `Hap5=C/A|G|C` first, n=7, total
score `4.2581`. Robust mode still ranks exact Qin2014 `Hap1=A|G|C` first,
n=371, total score `1.8132`, sample reliability `0.9488`, matching the
functional haplotype.

HTML verification:
Browser geometry on the robust report showed no remaining `{gwas_*}`
placeholders and no console errors. The GWAS vertical lines now end at the
panel bottom within ~1 px, the gene-structure upward guide lines start at the
GWAS panel bottom instead of crossing above the border, the LD canvas left
edge matches the first visible sequence column (`0 px` delta), and the LD
right edge is within `1 px` of the last visible sequence column. The score
panel now leaves a fixed `15 px` gap before the sequence columns, preventing
the score chart width from pushing the LD triangle to the right.

2026-06-20 Post-GWAS evidence detail layout/bridge-score rerun:
The same default and robust commands were rerun after adapting the newer
auditable evidence detail panels. This was a report-layout and field-mapping
fix only; the discovery score mode and literature validation conclusion were
not changed.

Commands:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP`

`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --score-mode robust_discovery`

Outputs:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP/TaGW2-B1-remoteSNP.html`
and
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery/TaGW2-B1-remoteSNP.html`.

Top-scored haplotypes:
Default mode still ranks rare ambiguous `Hap5=C/A|G|C` first, n=7, total
score `4.2581`. Robust mode still ranks exact Qin2014 `Hap1=A|G|C` first,
n=371, total score `1.8132`, sample reliability `0.9488`, matching the
functional haplotype.

HTML verification:
The new `Top Variant Evidence`, `Variant-Haplotype Bridge`, and
`Confidence Flags` panels are now inside a collapsed
`Review detailed local evidence` disclosure by default. Browser geometry on
the robust report showed the default evidence panel height was reduced from
the previous ~625 px to ~196 px, keeping the network/GWAS row near the first
screen. When expanded, the detail panel remains above the GWAS/gene guide-line
corridor. The bridge now reads haplotype scores from `total` when
`overall_score`/`score` are absent; the robust report shows
`score=Hap1 (1.813)` for the lead evidence rows instead of `score=NA`.

Reliability limits:
The biological validation remains limited to the available 3 promoter SNPs
and 816 Watkins TGW-overlap samples; the current report still has no external
GWAS file input. Literature variants are used only for post-ranking audit,
not as discovery-scoring inputs.

2026-06-21 discovery-safe candidate panel rerun:
The same default and robust commands were rerun after adding the report
features selected for discovery use: `Discovery Candidate List`,
`Score Component Breakdown`, and `Reliability & Population`. Literature
haplotypes and published variant labels are not used by these panels; they
rank and explain candidates only from local haplotype scores, effect summaries,
sample counts, sample reliability, ambiguity penalties, and population/sample
coverage summaries.

Commands:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP`

`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --score-mode robust_discovery`

Outputs:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP/TaGW2-B1-remoteSNP.html`
and
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery/TaGW2-B1-remoteSNP.html`.

Top-scored haplotypes:
Default mode still ranks rare ambiguous `Hap5=C/A|G|C` first, n=7, total
score `4.2581`; exact Qin2014 `Hap1=A|G|C` is second, n=371, score `2.5746`.
Robust mode ranks exact Qin2014 `Hap1=A|G|C` first, n=371, total score
`1.8132`, sample reliability `0.9488`, ambiguity factor `1.0000`; rare
ambiguous `Hap5` drops to third, n=7, score `0.9733`, sample reliability
`0.2593`, ambiguity factor `0.7333`.

HTML verification:
Both regenerated integrated reports contain `Local Candidate Evidence` plus
the three discovery-safe panels. `Score Component Breakdown` is rendered as a
compact `component-matrix` instead of tall per-haplotype cards. Browser
geometry on the robust report showed the candidate strip height was reduced
from about `1014 px` to about `374 px`, the full evidence panel from about
`1220 px` to about `580 px`, and `Review detailed local evidence` remains
collapsed by default. The existing horizontal scroll comes from the fixed-width
large integrated figure canvas, not from these new candidate panels.

2026-06-21 score-mode-aware candidate panel rerun:
The same default and robust commands were rerun after fixing a report UI
consistency bug: the `Discovery Candidate List`, `Score Component Breakdown`,
and `Reliability & Population` panels are now pre-rendered per
`score_mode + phenotype` and updated by the same Original/Robust and phenotype
switches as the score plot. This is a report-state fix only; the discovery
ranking inputs and literature audit conclusion are unchanged.

Commands:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP`

`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --score-mode robust_discovery`

Outputs:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP/TaGW2-B1-remoteSNP.html`
and
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery/TaGW2-B1-remoteSNP.html`.

Top-scored haplotypes:
Default candidate panel ranks rare ambiguous `Hap5=C/A|G|C` first, n=7,
total score `4.2581`; exact Qin2014 `Hap1=A|G|C` is second, n=371, score
`2.5746`. Robust candidate panel ranks exact Qin2014 `Hap1=A|G|C` first,
n=371, total score `1.8132`, sample reliability `0.9488`; rare ambiguous
`Hap5` is third, n=7, score `0.9733`, sample reliability `0.2593`.

HTML/browser verification:
Static inspection of the robust report confirmed `allDiscoveryCandidatePanelData`
contains both `robust_discovery` and sibling `default` panel data, and
`switchScoreMode()` / `switchPhenotype()` call `updateDiscoveryCandidatePanels()`.
In the in-app browser, after reload the robust report opened with
`Score mode: Robust` and candidate/component/reliability first row `Hap1`.
Clicking `Original` changed all three panels to `Score mode: Original` with
first row `Hap5`; clicking `Robust` changed all three panels back to `Hap1`.
This confirms the new discovery-safe explanatory panels no longer show stale
rank/component/reliability data after mode switches.

2026-06-21 Candidate Key Sites discovery panel rerun:
The same default and robust commands were rerun after adding
`Candidate Key Sites`, a discovery-stage panel that asks which displayed
variant sites distinguish the current top-scored haplotype from the other
haplotypes. The panel uses only local haplotype sequences, phenotype means,
sample counts, local marker P values when available, and variant annotation.
It does not use Qin2014 or any other literature functional labels as scoring
or ranking input.

Commands:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP`

`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --score-mode robust_discovery`

Outputs:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP/TaGW2-B1-remoteSNP.html`
and
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery/TaGW2-B1-remoteSNP.html`.

Robust `Candidate Key Sites` result:
The robust report opens with `Hap1=A|G|C` as the top-scored haplotype.
The key-site panel ranks three promoter SNPs:
`291,759,689` (`A vs C=437; C/A=7`, n `371/444`, site-level contrast `2.457`,
P `6.98e-09`, stable), `291,761,315`
(`C vs T=73; T/C=1`, n `371/74`, site-level contrast `5.165`,
shared-allele warning), and `291,760,677`
(`G vs A=7`, n `371/7`, site-level contrast `-1.728`,
low-other-n/shared-allele/low-MAF warning). These correspond to the three
displayed TaGW2-B1 promoter SNPs and show why the top robust haplotype is
locally distinguishable from the other haplotypes.

Default `Candidate Key Sites` result:
Switching to Original changes the same panel to `Hap5=C/A|G|C`, the default
top-scored haplotype. The panel flags `low top n` because `Hap5` has only
7 carriers (`7/809` contrast comparison), which is consistent with the known
default-mode over-ranking of a rare ambiguous haplotype.

HTML/browser verification:
Static inspection confirmed the robust HTML embeds
`top_haplotype_key_sites_html`, `top_haplotype_key_site_positions`,
`highlightKeySiteColumns()`, and the sequence-column `key-site-highlight`
CSS. In the in-app browser, the robust page showed `Candidate Key Sites`
with first row `Hap1`; three positions were highlighted in the haplotype
sequence table. A coordinate click on `Original` changed the key-site panel
to `Hap5` while retaining the same three highlighted display positions,
confirming that key-site rows update with score mode. Screenshot capture of
the very large integrated page timed out in the browser runtime, but DOM and
interaction checks succeeded.

2026-06-21 Candidate Key Sites code-review fix:
The default and robust commands above were rerun after fixing the review
findings in the discovery key-site panel. `phenotype_contrast` is now computed
per site as top-haplotype/top-allele carriers versus samples carrying a
different allele at that specific site, instead of reusing one global top-hap
versus non-top-hap contrast for every row. Candidate-key-site P values can now
be selected per phenotype through `variant_pvalues_by_phenotype`. Low top n,
low different-allele n, shared alleles, high missingness, low MAF, and missing
phenotype evidence now reduce `priority_score` through a reliability factor
instead of only changing the warning label. String marker IDs are preserved in
`top_haplotype_key_site_positions`, and report flag CSS classes are whitelisted
to `pass`, `warn`, or `muted`.

Verification:
`python -m unittest test_star_gene_data.py -v` passed 71 tests;
`python -m py_compile haplotype_phenotype_analysis.py star_gene_validation.py`
passed; `python run_rice_test.py` completed with the known
`LOC_Os01g02660 no_phenotype_match` plus two successful rice targets. Browser
DOM verification on the robust TaGW2-B1 report confirmed `Score mode: Robust`,
phenotype `TGW_mean`, key-site rows for positions `291759689`, `291761315`,
and `291760677`, and highlighted sequence columns for the same positions.

2026-06-21 Discovery Candidate List flag layout rerun:
The default and robust commands above were rerun after fixing truncated
reliability flags in the `Discovery Candidate List`. Long notes such as
`low n, low reliability` are now rendered as separate wrapped chips instead
of one clipped badge, and the candidate/table flag columns allow normal
wrapping. This is a report-layout fix only; the top-scored haplotype and
literature validation conclusion did not change.

Outputs:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP/TaGW2-B1-remoteSNP.html`
and
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery/TaGW2-B1-remoteSNP.html`.

Verification:
`python -m unittest test_star_gene_data.py -v` passed 72 tests;
`python -m py_compile haplotype_phenotype_analysis.py star_gene_validation.py`
passed. Browser DOM geometry on the robust report confirmed every candidate
flag cell had `white-space: normal`, no horizontal overflow, no chip overflow,
and no overlap into the `Seq` column.

### VRN-Kiss2014

Literature positive control:
Kiss et al. 2014 diagnostic marker states for `VRN-A1`, `VRN-B1`, and
`VRN-D1`; the all-spring/dominant marker combination is `1|1|1`.

Data:
`star_gene_database/wheat_nature_2024/VRN-Kiss2014`, 676 complete samples,
3 VRN markers, DEV49/DEV59 heading-date phenotypes.

Default result:
Top-scored haplotype was `Hap7=1|1|1`; literature audit matched all three
spring marker states plus the combined diagnostic haplotype. However, `Hap7`
has only 2 samples, so this is a useful audit check but not strong
population-statistical proof.

Robust result:
Command:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-Kiss2014 --score-mode robust_discovery`

Output:
`star_gene_results/wheat_nature_2024/VRN-Kiss2014__robust_discovery/VRN-Kiss2014.html`

For both `DEV49_mean` and `DEV59_mean`, robust mode ranked `Hap2=0|0|1`
first rather than the all-spring `Hap7=1|1|1`. `Hap2` has 31 samples and
carries the `VRN-D1=1` spring/dominant diagnostic state. `Hap7` has only
2 samples and was demoted by the sample-reliability correction.

`DEV49_mean` robust ranking:
`Hap2` score `1.6641`, n=31, mean `209.677`, reliability `0.6078`;
`Hap7` score `1.0164`, n=2, mean `214.750`, reliability `0.0909`.

`DEV59_mean` robust ranking:
`Hap2` score `1.6537`, n=31, mean `218.371`, reliability `0.6078`;
`Hap7` score `1.0120`, n=2, mean `222.750`, reliability `0.0909`.

Literature audit:
`VRN-D1_spring_diagnostic_marker` is `matched_top_haplotype`, while
`VRN-A1`, `VRN-B1`, and the combined all-three-spring haplotype are
`present_but_not_top`.

Interpretation:
This is not a failure of the discovery mode. It shows a useful tradeoff:
robust discovery ranking prefers a more stable single-locus `VRN-D1` signal
over a very rare all-three-spring literature combination. For future reporting,
VRN should be shown as evidence that stable-common rank and rare-candidate rank
may need to be reported separately.

### Rht-D1b

Literature positive control:
Canonical semi-dwarf `Rht-D1b` stop-gained allele, WheatOmics annotation
`c.181G>T`, `p.Glu61*` in `TraesCS4D01G040400`.

Data:
`star_gene_database/wheat_nature_2024/Rht-D1b`, Watkins CFLN06 plant-height
phenotype overlap.

Default result:
Top-scored haplotype contained the functional T allele and had lower plant
height direction, but carrier count was very small, so confidence is weak.

Robust result:
Command:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Rht-D1b --score-mode robust_discovery`

Output:
`star_gene_results/wheat_nature_2024/Rht-D1b__robust_discovery/Rht-D1b.html`

Robust result:
Top-scored haplotype remained the functional T-allele haplotype `Hap2`,
n=3, mean plant height `84.583`, score `1.0938`, reliability `0.1304`.
`Hap1` has n=797, mean `99.615`, score `0.0000`. Score regression
`R^2=0.0056`, `P=0.0334`. Literature audit status:
`Rht-D1b_stop_gained = matched_top_haplotype`.

2026-06-26 strict rerun:
Commands:
`python prepare_wheat2024_rht1_functional_snps.py --min-haplotype-count 1 --target Rht-D1b`

`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Rht-D1b --score-mode robust_discovery`

Outputs:
`star_gene_results/wheat_nature_2024/Rht-D1b__robust_discovery/Rht-D1b.html`

`star_gene_results/wheat_nature_2024/Rht-D1b__robust_discovery/haplotype_score.html`

The database rebuild used one exact stop-gained SNP and 802 phenotype-overlap
samples. Haplotypes are `Hap1=G` (n=797), `Hap2=T` (n=3), and `Hap3=G/T`
(n=2). The strict top-scored haplotype is `Hap2=T`, exactly matching the
published `Rht-D1b` allele. The audit row records allele counts
`G:797;G/T:2;T:3`, carrier count `5`, `top_haplotype_exact_expected=True`,
and `validation_status=matched_top_haplotype`.

The same run also shows the limitation: the support-shrunk directional top is
`Hap1=G`, not `Hap2=T`, because the true functional class has only three
pure-T top-haplotype samples. Therefore this result is valid evidence that the
scorer can recover the decisive functional variant when it is present, but it
is not strong enough to count as one of the teacher's robust three-proof set.

Interpretation:
Rht-D1b supports that robust mode can still recover a known functional
variant when the effect direction is strong, but it remains weak evidence
because only three top-haplotype carriers are present in the analyzed samples.

### Rht-Zanke2014

Literature positive control:
Zanke et al. 2014 PLOS ONE Table S2 provides variety-level candidate-gene
genotypes for `Rht-B1` and `Rht-D1` plus multi-environment plant-height
phenotypes. The marker states come from the Ellis et al. Rht perfect-marker
scheme: `B1b` and `D1b` are the semi-dwarf functional states.

Data:
`star_gene_database/wheat_nature_2024/Rht-Zanke2014`, 368 complete varieties,
2 diagnostic markers, 4 haplotypes.

Command:
`python prepare_wheat2024_rht_zanke2014.py`

Run command:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Rht-Zanke2014 --score-mode robust_discovery`

Output:
`star_gene_results/wheat_nature_2024/Rht-Zanke2014__robust_discovery/Rht-Zanke2014.html`

Result:
The original robust `total` score ranked the high-plant-height wild-type
haplotype `Hap2=B1a|Rht-D1a` highest. This was the source of the Rht-Zanke2014
concern: the HTML score plot was showing an effect-magnitude axis, not the
expected biological direction for a `decreases_trait` target.

The scorer now preserves raw `total` but uses `directional_total` as the main
score axis when `expected_direction` is known. For `PlantHeight_BLUE`,
`Hap1=B1a|D1b` is the top-scored haplotype (`directional_total=0.9642`, n=214,
mean `82.540 cm`), followed by `Hap3=B1b|Rht-D1a`
(`directional_total=0.5613`, n=26) and `Hap4=B1a|D1b/Rht-D1a`
(`directional_total=0.4813`, n=2). `Hap2=B1a|Rht-D1a` is low after direction
correction (`directional_total=0.3500`, n=126, mean `96.384 cm`). Score
regression is `R^2=0.451`, `P=1.36e-49`, slope `-21.054`.

For `PlantHeight_GAT_2012`, `Hap1=B1a|D1b` remains top
(`directional_total=0.9624`, n=214), while `Hap3` is lower-ranked despite a
more extreme environment-specific mean because support shrinkage penalizes the
smaller group. Score regression is `R^2=0.532`, `P=1.71e-60`, slope `-33.534`.

Literature audit:
`Rht-D1b_diagnostic_marker` is now `matched_top_haplotype` and
`matched_directional_top_haplotype`: top haplotype `Hap1=B1a|D1b`, n=214, all
214 top-haplotype samples carry the expected `D1b` state. `Rht-B1b` remains
present but not top because the B1b class is `Hap3=B1b|Rht-D1a`, n=26.

Interpretation:
Rht-Zanke2014 now supports the scoring method after fixing the direction
handling. The haplotype partition itself matched the literature counts; the
problem was interpreting raw magnitude score as a biological top for a
height-decreasing allele.

### 2026-06-23 Robust full-region rerun and cross-target algorithm check

Purpose:
The full TaGW2-B1 regional SNP panel split the literature promoter haplotype
into many background haplotypes. The algorithm was therefore updated with two
generic safeguards, not target-specific rules: full-sequence position indexing
for regional haplotypes, robust core-position grouping, and direction-aware
stable top selection with low-support extreme filtering. Literature variants
remain post hoc validation labels only and are not used in discovery scoring.

Run command:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --target VRN-D1-Kiss2014 --target Rht-Zanke2014 --score-mode robust_discovery`

Output directories:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery`,
`star_gene_results/wheat_nature_2024/VRN-D1-Kiss2014__robust_discovery`,
and
`star_gene_results/wheat_nature_2024/Rht-Zanke2014__robust_discovery`.

TaGW2-B1 full-region result:
The database has 756 phenotype-overlap samples, 71 regional SNPs, and 112
full-region haplotypes. Score regression remains significant for TGW_mean
(`R^2=0.0124`, `P=0.0022046`), and the output HTML was regenerated at
`TaGW2-B1-remoteSNP.html`. However, the Qin2014 `Hap-6B-1` promoter pattern
`A|G|C` is no longer the raw top or directional top after adding the full
regional background SNPs. The literature haplotype audit remains
`present_but_not_top`; the robust core group is not testable for the full
three-marker pattern because not all literature markers were selected as core
positions. This target is therefore a gene-level/phenotype-signal check under
the full-region database, not a strong functional-haplotype proof.

VRN-D1-Kiss2014 result:
Target gene/marker is `VRN-D1` diagnostic spring marker from Kiss et al. 2014.
Data source is the extracted Springer embedded workbook converted into
`star_gene_database/wheat_nature_2024/VRN-D1-Kiss2014`. Robust mode analyzed
676 samples, 1 diagnostic marker, and 2 haplotypes. Top-scored and
directional top haplotype are both `Hap2` with 38 samples and the expected
spring marker allele `1`. DEV49_mean score regression is `R^2=0.0382`,
`P=2.980468e-07`; DEV59_mean is `R^2=0.0230`, `P=7.488576e-05`.
Literature audit status is `matched_top_haplotype` and
`matched_directional_top_haplotype`. This is a usable independent positive
control, with low confidence caused by modest PVE rather than marker mismatch.

Rht-Zanke2014 result:
Target markers are `Rht-B1` and `Rht-D1` diagnostic states from Zanke et al.
2014 Table S2. Data source is
`star_gene_database/wheat_nature_2024/Rht-Zanke2014` with 368 complete
varieties, 2 markers, and 4 haplotypes. For `PlantHeight_BLUE`, raw top is
the tall wild-type-like core `B1A|RHT-D1A` (`Hap2`, n=126, mean 96.384 cm),
but expected direction is `decreases_trait`, so the stable directional top is
`Hap1=B1a|D1b` with n=214 and mean 82.540 cm. Score regression is
`R^2=0.0842`, `P=1.440617e-08`. Literature audit for
`Rht-D1b_diagnostic_marker` is `matched_directional_top_haplotype`. For
`PlantHeight_GAT_2012`, the stable low-phenotype directional top is
`Hap3=B1b|Rht-D1a`, n=25, mean 72.454 cm, with score regression
`R^2=0.1301`, `P=2.067409e-12`. This remains a usable positive control only
when interpreted with trait direction.

Current interpretation:
The generic changes did not make the method work only for GW2: VRN-D1 and
Rht-Zanke2014 still validate after the same robust rerun. The full-region
TaGW2-B1 result also clarifies a limitation: when many background SNPs are
included, exact literature haplotype matching can fail even when the target
gene still has a significant score-phenotype signal. Therefore future unknown
gene discovery reports should show raw top, direction-aware top, and robust
core groups separately, and should not equate a full-region exact haplotype
rank with recovery of a published promoter haplotype.

### 2026-06-24 Functional sub-haplotype robust rerun

Purpose:
The 71-SNP TaGW2-B1 full-region database showed that background SNPs can split
a known functional promoter haplotype into many exact regional haplotypes. The
algorithm was updated with a generic `functional_haplotype_groups` discovery
layer and support-shrunk direction-aware ranking. Literature variants remain
post hoc audit labels only and are not scoring inputs.

Implemented generic changes:
- `robust_discovery` now scores functional sub-haplotype groups selected from
  local annotation, phenotype signal, EB marker effect, MAF/missingness, LD
  context, and near-boundary gene-body signal.
- `functional_marker` and `diagnostic_marker` annotations are treated as
  functional panel markers.
- Scoring now uses the full regional position set, while HTML can still display
  a smaller top-haplotype variable subset.
- Direction-aware top selection uses support-shrunk phenotype direction, so
  small high/low outlier groups do not dominate stable groups.

Run command:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --target VRN-D1-Kiss2014 --target Rht-Zanke2014 --score-mode robust_discovery`

Output directories:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery`,
`star_gene_results/wheat_nature_2024/VRN-D1-Kiss2014__robust_discovery`,
and
`star_gene_results/wheat_nature_2024/Rht-Zanke2014__robust_discovery`.

TaGW2-B1-remoteSNP result:
The database has 756 phenotype-overlap samples, 71 regional SNPs, and 112
full-region haplotypes. Functional positions selected for TGW_mean are
`291759689;291759713;291759935;291760409;291760469;291760535;291760602;291760677;291761315`,
which include all three Qin2014 audit markers for `Hap-6B-1`
(`291759689=A`, `291760677=G`, `291761315=C`). Score regression is significant
but modest (`R^2=0.0078`, `P=0.0153196`). Raw top functional group is
`C|T|C|T|G|G|G|G|C`, n=329, mean TGW `38.9632`; support-shrunk directional
top functional group is `A|T|C|C|G|G|G|G|C`, n=340, mean TGW `40.9845`.
The literature audit for `Qin2014_Hap-6B-1` reports
`directional_functional_group_validation_status=matched_directional_top_functional_group`.

VRN-D1-Kiss2014 result:
The database has 676 samples, 1 diagnostic marker, and 2 haplotypes.
Functional positions are `[3]`. The top and directional top functional group
is `1`, n=38, matching the spring/dominant `VRN-D1=1` marker. DEV49_mean
score regression is `R^2=0.0382`, `P=2.980468e-07`; DEV59_mean is
`R^2=0.0230`, `P=7.488576e-05`. Literature audit remains
`matched_top_haplotype` and `matched_directional_top_haplotype`.

Rht-Zanke2014 result:
The database has 368 samples, 2 diagnostic markers, and 4 haplotypes.
Functional positions are `[1, 2]`. For `PlantHeight_BLUE`, raw top functional
group is the tall wild-type-like `B1A|RHT-D1A`, n=126, mean `96.3842 cm`, but
the expected direction is `decreases_trait`, so the directional top functional
group is `B1A|D1B`, n=214. Literature audit for
`Rht-D1b_diagnostic_marker` remains `matched_directional_top_haplotype`.
For `PlantHeight_GAT_2012`, the directional top is `B1A|D1B`, n=205.

Current verdict:
This gives three usable positive controls under one generic scoring mode:
TaGW2-B1 full-region validates through direction-aware functional group,
VRN-D1 validates through diagnostic-marker top ranking, and Rht-D1 validates
through direction-aware plant-height ranking. The remaining limitation is that
TaGW2-B1 raw top score still favors a non-literature functional group, so raw
rank alone should not be used as the only discovery criterion for traits where
biological direction matters.

### 2026-06-25 Phenotype-free robust discovery rerun

Purpose:
The previous 2026-06-24 `robust_discovery` implementation used validation
phenotype signal, EB marker effects, effect-size score, and direction-adjusted
totals inside discovery ranking. That is not valid for prediction/discovery.
The algorithm was changed to a local PostGWAS-like site weighting model:
annotation severity, promoter/gene-boundary context, MAF/missingness, LD, and
explicitly external GWAS/eQTL/PostGWAS-like evidence can affect discovery
scores. The current validation phenotype and literature variants remain post
hoc validation only. The code does not call the PostGWAS website at runtime.

Implemented generic changes:
- `robust_discovery` component weights now exclude phenotype-derived
  `eb_effect` and `effect_size`.
- Locally generated per-position phenotype p-values are ignored by scoring
  unless a record is explicitly marked as external evidence.
- Functional/core position selection no longer uses `phenotype_logp` or
  `phenotype_effect`.
- Active `score_axis` is always `total`; `directional_total` remains in JSON
  only as a validation/audit helper.

Run command:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --target VRN-D1-Kiss2014 --target Rht-Zanke2014 --score-mode robust_discovery`

Output directories:
`star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery`,
`star_gene_results/wheat_nature_2024/VRN-D1-Kiss2014__robust_discovery`,
and
`star_gene_results/wheat_nature_2024/Rht-Zanke2014__robust_discovery`.

TaGW2-B1-remoteSNP result:
The database has 756 phenotype-overlap samples, 71 regional SNPs, and 112
full-region haplotypes. The phenotype-free functional positions are
`291759689;291760409;291760469;291760535;291760677;291761315`, which still
include all three Qin2014 diagnostic positions. Raw top full haplotypes are
rare high-weight promoter-background combinations, led by `Hap18` with total
`1.0882`. Raw top functional group is `C|T|G|G|G|C`, n=333, representative
`Hap2`. The Qin2014 `Hap-6B-1` combined promoter haplotype remains
`present_but_not_top` in raw phenotype-free ranking, although two of the three
single diagnostic SNPs are `matched_top_haplotype`. Score regression is low
(`R^2=0.0009`, `P=0.4020`), as expected after removing phenotype leakage.
Current verdict: useful site-level recovery, but not a full positive proof
without external GWAS/eQTL/PostGWAS weights or an allele-direction model.

VRN-D1-Kiss2014 result:
The database has 676 samples, 1 diagnostic marker, and 2 haplotypes. Functional
positions are `[3]`. Raw phenotype-free top is `Hap2` / functional sequence
`1`, n=38, matching the published spring/dominant `VRN-D1=1` marker for both
DEV49_mean and DEV59_mean. Score regression remains DEV49_mean `R^2=0.0382`,
`P=2.980468e-07`; DEV59_mean `R^2=0.0230`, `P=7.488576e-05`. Current verdict:
positive proof under phenotype-free discovery.

Rht-Zanke2014 result:
The database has 368 samples, 2 diagnostic markers, and 4 haplotypes.
Functional positions are `[1, 2]`. Raw phenotype-free full-haplotype top is
`Hap3`, not the stable literature plant-height group. However the functional
group top by support-shrunk rank is `B1A|D1B`, n=214 for PlantHeight_BLUE and
n=205 for PlantHeight_GAT_2012, matching the Rht-D1b diagnostic marker group.
The audit still reports `Rht-B1b_diagnostic_marker=matched_top_haplotype` and
`Rht-D1b_diagnostic_marker=present_but_not_top` at single-marker raw top level.
Current verdict: partial positive proof at functional-group level; raw
full-haplotype top alone is not sufficient.

Current interpretation:
This rerun fixes the scientific leakage problem. It gives one strong
phenotype-free proof (VRN-D1), one partial functional-group proof (Rht-Zanke),
and one site-level but not haplotype-level recovery (TaGW2-B1). To reach the
teacher's "at least 3 positive controls" standard without phenotype leakage,
the next data improvement should add external GWAS/eQTL/PostGWAS-like weights
or richer functional annotations, not reintroduce validation phenotype into
the discovery score.

### 2026-06-26 VRN WheatOmics SNP VCF download

Purpose:
The previous VRN positive control used Kiss2014 diagnostic marker states, not
nucleotide-level VCF genotypes. To avoid using a multi-gene marker panel as
single-gene evidence, VRN-A1, VRN-B1, and VRN-D1 were downloaded as separate
single-gene SNP micro-VCFs from the WheatOmics 1051-accession merged SNP VCF.

Data source:
`https://wheatomics.sdau.edu.cn/jbrowse-1.12.3-release/Chinese_Spring1.0/tracks/vcf/merge.SNP.Missing-unphasing.ID.ann.finalSID.allele2_retain.hard_retain.InbreedingCoeff_retain.clean.anno.vcf.gz`

Reproducible download command:
`python download_wheat2024_vrn_remote_snps.py`

Output directory:
`external_data/wheat_nature_2024/vrn_remote`

Downloaded targets:
- `VRN-A1`, gene `TraesCS5A01G391700`, region
  `chr5A:587409454-587425416`, output
  `VRN-A1.wheatomics_snp.vcf.gz`, 97 SNP records, 1051 samples.
- `VRN-B1`, gene `TraesCS5B01G396600`, region
  `chr5B:573800883-573818070`, output
  `VRN-B1.wheatomics_snp.vcf.gz`, 182 SNP records, 1051 samples.
- `VRN-D1`, gene `TraesCS5D01G401500`, region
  `chr5D:467174609-467186508`, output
  `VRN-D1.wheatomics_snp.vcf.gz`, 62 SNP records, 1051 samples.

Score mode:
Not run yet. This entry records data acquisition only.

Top-scored haplotype:
Not available yet.

Post hoc literature match:
Not available yet. Literature VRN decisive alleles often include promoter,
intron-1 deletion, and copy-number/structural variation. These downloaded
files are SNP-only micro-VCFs, so they can support nucleotide-level SNP
haplotype scoring but cannot by themselves prove recovery of a causal VRN
deletion/SV/CNV allele unless that allele is separately added from an INDEL/SV
source.

Reliability limitation / blocked reason:
Single-gene VCF data are now available for VRN-A1/B1/D1, but the current data
do not yet include VRN INDEL/SV/CNV calls. The next step is to build separate
VCF-derived databases for each VRN gene and run robust discovery scoring
without combining VRN-A1, VRN-B1, and VRN-D1 into one haplotype.

### 2026-06-26 VRN WheatOmics SNP VCF robust-discovery run

Purpose:
Test the newly downloaded VRN SNP-only micro-VCFs as separate single-gene
targets. This avoids the earlier Kiss2014 multi-marker panel issue and does
not combine VRN-A1, VRN-B1, and VRN-D1 into one haplotype.

Preparation command:
`python prepare_wheat2024_vrn_remote_snps.py`

Analysis command:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-A1-remoteSNP --target VRN-B1-remoteSNP --target VRN-D1-remoteSNP --score-mode robust_discovery`

Phenotype:
Watkins CFLN06 growth habit was converted from `GrwHabit_E_sw-CFLN06` into
`GrowthHabitSpringScore_CFLN06`, defined as the fraction of `s` calls in the
spring/winter string. Non-`s/w` values were excluded. `PlantHeight_CFLN06` was
retained only as a covariate/auxiliary column, not as the VRN validation trait.

Output directories:
- `star_gene_results/wheat_nature_2024/VRN-A1-remoteSNP__robust_discovery`
- `star_gene_results/wheat_nature_2024/VRN-B1-remoteSNP__robust_discovery`
- `star_gene_results/wheat_nature_2024/VRN-D1-remoteSNP__robust_discovery`

VRN-A1-remoteSNP result:
The database has 96 SNPs, 18 retained haplotypes, and 663 phenotype-overlap
samples. Haplotype association with spring score is significant but modest:
corrected P `2.1114e-02`, PVE `10.94%`. The top robust-scored haplotype is
`Hap4`, n=39, mean spring score `0.9679`, score `1.1390`; direction is reported
as consistent. However score-vs-phenotype regression is not significant
(`R^2=0.0000`, P `0.8834`). Current verdict: weak/partial positive signal for
ranking a spring-like haplotype, but not a strong proof that the generic score
tracks phenotype across haplotypes.

VRN-B1-remoteSNP result:
The database has 127 SNPs, 26 retained haplotypes, and 511 phenotype-overlap
samples. Haplotype association is strong: corrected P `4.7736e-06`, PVE
`17.18%`, and score-vs-phenotype regression is significant (`R^2=0.0415`, P
`3.4551e-06`). But the slope is negative and the audit reports direction as
inconsistent. The top robust-scored haplotype is `Hap3`, n=15, mean spring score
`0.4333`, score `1.3437`; the top score therefore favors a more winter-like
haplotype rather than the expected spring-like direction. Current verdict:
clear regional SNP signal, but not a positive proof of correct discovery
direction under the current robust score.

VRN-D1-remoteSNP result:
The database has 58 SNPs, 21 retained haplotypes, and 712 phenotype-overlap
samples. Haplotype association is significant: corrected P `7.0250e-05`, PVE
`8.26%`. The top robust-scored haplotype is `Hap12`, n=4, mean spring score
`1.0000`, score `1.1886`; direction is reported as consistent. But the top
haplotype is very rare and score-vs-phenotype regression is not significant
(`R^2=0.0036`, P `0.1094`). Current verdict: partial association evidence, not
a strong proof.

Overall interpretation:
The SNP-only VRN VCF run does not yet give three strong proof cases. It shows
that VRN regional SNP haplotypes carry growth-habit signal, especially at
VRN-B1 and VRN-D1, but the discovery score has two limitations: rare
high-scored haplotypes can dominate, and the score direction can disagree with
the expected biological direction. Because the downloaded VRN data are SNP-only
and do not include the known promoter/intron-1 deletion, CNV, or SV causal
alleles, these results should be treated as regional linked-SNP validation, not
as exact literature functional-variant validation.

### 2026-06-26 WheatOmics public INDEL micro-VCF supplement

Purpose:
Download the other publicly reachable INDEL data needed to check whether the
current star-gene targets can be augmented beyond SNP-only genotypes. The goal
was to find target-region INDEL records for VRN, Rht, and GW2 without
downloading multi-GB whole-genome VCFs.

Downloaded/queried source list:
- WheatOmics JBrowse track list:
  `https://wheatomics.sdau.edu.cn/jbrowse-1.12.3-release/Chinese_Spring1.0/trackList.json`
- WheatOmics VCF directory:
  `https://wheatomics.sdau.edu.cn/jbrowse-1.12.3-release/Chinese_Spring1.0/tracks/vcf/`
- `WEC_INDEL_IWGSCv1.0.eff.vcf.gz`
- `GBS_filtered_Indels_IWGSCv1.0.eff.vcf.gz`
- `WildEmmer10WGS_INDEL_eff.vcf.gz`

Reproducible extraction command:
`python download_wheat2024_indel_microvcfs.py`

Output directory:
`external_data/wheat_nature_2024/indel_microvcfs`

Status table:
`external_data/wheat_nature_2024/indel_microvcfs/microvcf_status.tsv`

Main findings:
- VRN-A1, VRN-B1, and VRN-D1 target regions had 0 INDEL records in all three
  queried WheatOmics public INDEL tracks. This does not solve the VRN
  causal-SV/CNV blocker.
- Rht-B1 had INDEL records in WEC and WildEmmer tracks, including one WEC
  frameshift insertion at `chr4B:30862060` in `TraesCS4B01G043100`, but this
  source has 62 samples and 0 Watkins phenotype-sample overlap.
- TaGW2-A1/B1/D1 had a few INDEL records in GBS or WildEmmer tracks, mostly
  intronic/upstream records, but all queried INDEL sources had 0 Watkins sample
  overlap.
- The queried INDEL tracks therefore cannot be merged into the current Watkins
  phenotype validations as population-level evidence. They are useful only as
  external regional-variant annotations unless a matching phenotype table is
  obtained for those INDEL source panels.

WWWG2B/Earlham status:
WWWG2B API calls that previously supplied OneDrive URLs now return HTTP 526 in
this environment, including `availableTable`, `fileTable`, and
`get_download_url_form_onedrive`. Earlham OpenData requests for WatSeq raw VCF
directories still return `Request Rejected` HTML. Therefore the still-needed
1051-sample WatSeq INDEL/SV/CNV files could not be fetched automatically in
this pass.

Current interpretation:
The additional reachable INDEL data have been downloaded as target-region
micro-VCFs, but they do not provide a stronger proof set because they are not
matched to the Watkins phenotype samples. For VRN, the decisive promoter or
intron-1 deletion/CNV/SV data remain missing from the current accessible
sources.

### 2026-06-26 Rht-D1-Zanke2014 single-marker rerun

Purpose:
Respond to the teacher's concern that `Rht-Zanke2014` combines `Rht-B1` and
`Rht-D1`. This rerun keeps only the Zanke2014 `Rht-D1` diagnostic marker
column, so the haplotypes no longer encode multiple Rht genes.

Literature validation object:
`Rht-D1b` diagnostic marker state from Zanke et al. 2014 Table S2, derived
from the Ellis et al. marker system. This is marker-level evidence for
Rht-D1b, not the base-level `chr4D:18781242 G>T` stop-gained VCF proof.

Data source:
`external_data/literature/rht_zanke2014/Table_S2_candidate_genes_phenotypes.xlsx`

Preparation command:
`python prepare_wheat2024_rht_zanke2014.py --single-marker-targets`

Database output:
`star_gene_database/wheat_nature_2024/Rht-D1-Zanke2014`

Analysis command:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Rht-D1-Zanke2014 --score-mode robust_discovery`

HTML output:
`star_gene_results/wheat_nature_2024/Rht-D1-Zanke2014__robust_discovery/Rht-D1-Zanke2014.html`

Haplotype split:
- `Hap1=D1b`, n=214.
- `Hap2=Rht-D1a`, n=152.
- `Hap3=D1b/Rht-D1a`, n=2.

Score/audit result:
- For `PlantHeight_BLUE`, score regression `R^2=0.3932`, `P=1.3256e-41`.
- For `PlantHeight_GAT_2012`, score regression `R^2=0.4406`, `P=1.0679e-46`.
- Raw robust total top is `Hap3=D1b/Rht-D1a`, score `0.2409`, n=2.
- Direction-aware top is exact `Hap1=D1b`, score `0.2`, n=214.
- `literature_variant_audit.csv` reports allele counts
  `D1b:214;D1b/Rht-D1a:2;Rht-D1a:152`, carrier count `216`, and
  `directional_validation_status=matched_directional_top_haplotype`.

Interpretation:
This target is the cleanest current Rht marker-level validation because it
uses one gene/marker only and recovers the known shorter-height `D1b` class
as the stable direction-aware top. It should be presented alongside the strict
base-level `Rht-D1b` SNP result: the SNP result is exact but weak because
carrier count is tiny, while this Zanke2014 single-marker result is strong at
marker level but not nucleotide-level VCF evidence.

### 2026-06-26 WheatOmics Zhai indel marker annotation supplement

Purpose:
Download the remaining small, reachable WheatOmics indel marker annotation
track for the same teacher star-gene targets, including flanking sequence,
without treating marker annotations as sample-level genotypes.

Reproducible extraction command:
`python download_wheat2024_indel_marker_annotations.py`

Output directory:
`external_data/wheat_nature_2024/zhai_indel_markers`

Status table:
`external_data/wheat_nature_2024/zhai_indel_markers/zhai_indel_marker_status.tsv`

Main findings:
- The script queries the WheatOmics `Indel_marker_from_zhai` JBrowse track and
  saves raw `trackData.json` plus relevant lazy chunks for each chromosome.
- With the default 1 Mb flank, `TaGW2-B1` has 8 marker rows near the target
  interval and `VRN-B1` has 7 marker rows near the target interval.
- `Rht-B1`, `Rht-D1`, `TaGW2-A1`, `TaGW2-D1`, `VRN-A1`, and `VRN-D1` have 0
  marker rows in the same 1 Mb flanking query.
- These rows contain marker coordinates, indel size, polymorphism ratio, and
  primers, but no per-accession genotype matrix. Therefore they are not direct
  inputs to discovery scoring and cannot prove haplotype-score effectiveness
  in the current Watkins phenotype panel.

Current interpretation:
This completes the reachable small-data supplement from WheatOmics for now:
public INDEL micro-VCFs plus Zhai marker annotations are downloaded, but the
still-needed validation-strengthening data are the 1051-sample WatSeq
INDEL/SV/CNV genotype files or another sample-matched genotype/phenotype table.

### 2026-06-26 Phenotype-free site-weighted robust discovery update

Purpose:
Revise `robust_discovery` toward a PostGWAS-like discovery model without using
the current validation phenotype in scoring. Site weights are computed from
variant annotation, explicit external evidence, MAF stability, missingness, LD
pruning, and gene-structure context. Current haplotype means, local trait
effect sizes, local p-values, and validation phenotype direction are excluded
from discovery scoring and remain post hoc validation/audit fields only.

Code changes:
- Added the `site_weighted` component to robust discovery scoring.
- Added auditable `site_weights` and `site_weighting_policy` fields to
  `haplotype_scores.json`; every site-weight row records
  `current_phenotype_used=False`.
- Kept rare high-impact functional sites when MAF is below the ordinary
  background threshold.
- Fixed the snpEff annotation alias `stop_gained -> stop_gain`; this was why
  the strict `Rht-D1b` stop-gained SNP had no site-weight record in the first
  implementation pass.

Regression tests:
`python -m unittest test_star_gene_data.StarGeneDataTests.test_site_weight_keeps_rare_high_impact_functional_sites test_star_gene_data.StarGeneDataTests.test_site_weight_keeps_snp_eff_stop_gained_annotation test_star_gene_data.StarGeneDataTests.test_robust_discovery_site_weighted_score_is_phenotype_free test_star_gene_data.StarGeneDataTests.test_robust_discovery_score_is_unchanged_when_phenotype_is_inverted test_star_gene_data.StarGeneDataTests.test_robust_discovery_uses_only_explicit_external_site_evidence test_star_gene_data.StarGeneDataTests.test_robust_discovery_adds_boundary_gene_body_signal_to_functional_groups test_star_gene_data.StarGeneDataTests.test_robust_discovery_penalizes_tiny_ambiguous_haplotypes test_star_gene_data.StarGeneDataTests.test_expected_decrease_direction_reorders_raw_high_phenotype_score -v`

Analysis command:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --target VRN-D1-Kiss2014 --target Rht-D1b --target Rht-Zanke2014 --score-mode robust_discovery`

Updated validation results:
- `TaGW2-B1-remoteSNP`: 756 samples, 71 regional SNPs, 112 haplotypes.
  `site_weights=39`, all phenotype-free. Raw top remains `Hap2`
  (`total=1.5448`, `site_weighted=0.6662`), while the audit still supports the
  Qin2014 promoter haplotype through directional/function-group evidence rather
  than raw full-region top rank. Score regression: `R^2=0.0031`,
  `P=0.12597`. This is a useful full-region stress test, not a raw-top proof.
- `VRN-D1-Kiss2014`: one diagnostic marker. `site_weights=1`,
  `Hap2` is both raw and directional top (`total=2.0155`,
  `site_weighted=0.6552`). Audit: `matched_top_haplotype` and
  `matched_directional_top_haplotype`. Score regression: `R^2=0.0382`,
  `P=2.98e-07`. This remains a strong marker-level positive proof.
- `Rht-D1b`: exact base-level stop-gained SNP `chr4D:18781242 G>T` now has one
  site-weight row (`annotation=stop_gain`, `MAF=0.003741`,
  `current_phenotype_used=False`). Raw and directional top in the regenerated
  score JSON are `Hap2=T` (`total=1.2372`, `site_weighted=0.1304`). Audit:
  `matched_top_haplotype`, but carrier support is still tiny
  (`Hap2` n=3; total T carriers 5). Score regression: `R^2=0.0056`,
  `P=0.0338`. This is exact but weak evidence.
- `Rht-Zanke2014`: marker-panel target with two functional markers.
  `site_weights=2`. Raw top is `Hap3` for the first trait, but the Rht-D1b
  literature marker remains recovered by direction-aware validation as
  `Hap1` (`matched_directional_top_haplotype`). Because this target combines
  B1 and D1 marker states, use it only as marker-panel/directionality evidence,
  not strict single-gene proof.

Current interpretation:
The site-weighted robust update makes the discovery score phenotype-free while
retaining a PostGWAS-like weighted-site mechanism. It fixes the real Rht-D1b
annotation-loss bug and provides clearer audit fields in the HTML/JSON. The
current positive set is still strongest for `VRN-D1-Kiss2014`, usable but weak
for exact `Rht-D1b`, and conditional for `TaGW2-B1`/`Rht-Zanke2014` because
their strongest support depends on direction-aware or functional-group
post hoc validation rather than raw full-region top rank.

Code-review follow-up:
A review found that `gwas_data` required explicit external evidence flags, but
`variant_info` fields such as `minus_log10_p`, `gwas_pvalue`, and `site_score`
were still accepted as site-weight evidence without the same external/source
gate. This has been fixed: `variant_info` score/logp/p-value fields are used
only when the variant record is explicitly external (`external=True` or an
external `source_type`/`evidence_type`). Regression coverage was added with
`test_site_weight_ignores_unmarked_variant_info_pvalues`.

Follow-up verification:
- `python -m unittest test_star_gene_data.py -v` passed 90 tests.
- `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --target VRN-D1-Kiss2014 --target Rht-D1b --target Rht-Zanke2014 --score-mode robust_discovery` completed.
- In the rerun, all four target `site_weights` have
  `current_phenotype_used=False`; their `external_weight` values are 0 unless
  explicitly external evidence is supplied. The previously recorded validation
  conclusions are unchanged: `VRN-D1-Kiss2014` remains the strongest positive,
  exact `Rht-D1b` remains weak but recovered, and TaGW2/Rht-Zanke2014 remain
  direction-aware or functional-group evidence rather than raw-top-only proof.
### 2026-06-27 Base-level genotype scoring rerun

Purpose:
Tighten the validation standard to true base-level genotype calls. From this
point, strict proof requires per-accession REF/ALT SNP/INDEL/SV genotypes.
Marker labels such as `D1b`, `VRN-D1=1`, `P/W/PAR`, or paper haplotype names
are recorded as marker-level or paper-haplotype evidence only.

Preparation commands:
`python prepare_wheat2024_vrn_remote_snps.py`

`python prepare_wheat2024_rht1_functional_snps.py --min-haplotype-count 1 --target Rht-D1b`

`python prepare_wheat2024_tagw2_b1_remote_snp.py --vcf D:\Desktop\data\GW2\chr6B.HARD.SNP.Missing-unphasing.ID.ann.finalSID.1047.allele2_retain.hard_retain.InbreedingCoeff_retain.vcf.gz --min-haplotype-count 1 --max-missing-rate 0.2`

Analysis commands:
`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-A1-remoteSNP --target VRN-B1-remoteSNP --target VRN-D1-remoteSNP --score-mode robust_discovery`

`python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Rht-D1b --target TaGW2-B1-remoteSNP --score-mode robust_discovery`

Output directories:
- `star_gene_results/wheat_nature_2024/VRN-A1-remoteSNP__robust_discovery`
- `star_gene_results/wheat_nature_2024/VRN-B1-remoteSNP__robust_discovery`
- `star_gene_results/wheat_nature_2024/VRN-D1-remoteSNP__robust_discovery`
- `star_gene_results/wheat_nature_2024/Rht-D1b__robust_discovery`
- `star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery`

Target results:

- `Rht-D1b`: exact causal base SNP `chr4D:18781242 G>T`, `c.181G>T`,
  `p.Glu61*`. Database has 802 samples and 1 SNP. Robust top is `Hap2=T`,
  n=3, with total T carriers 5. `site_weights=1`,
  `current_phenotype_used=False`, score regression `R^2=0.0056`,
  `P=0.0338`. Literature audit: `Rht-D1b_stop_gained=matched_top_haplotype`.
  Verdict: strict base-level proof, but weak because carrier support is tiny.
- `TaGW2-B1-remoteSNP`: full-region base SNP panel from complete chr6B WWWG2B
  VCF. Database has 756 samples, 71 SNPs, and 112 haplotypes. Qin2014
  promoter states are real bases at `chr6B:291759689 C>A`,
  `chr6B:291760677 G>A`, and `chr6B:291761315 T>C`. Robust raw top is
  `Hap2`, n=133; score regression `R^2=0.0031`, `P=0.12597`. Literature audit:
  `Qin2014_Hap-6B-1` is `present_but_not_top` at exact full-region haplotype
  level but `matched_directional_top_haplotype` through directional top
  `Hap1`. Verdict: base-level regional evidence and useful stress test, not
  raw-top-only proof.
- `VRN-A1-remoteSNP`: SNP-only VCF region, 663 samples, 96 SNPs, 18 haplotypes.
  Top is `Hap4`, n=39; score regression `R^2=0.0081`, `P=0.02045`.
  Verdict: real-base regional SNP signal, but not causal VRN SV/CNV proof.
- `VRN-B1-remoteSNP`: SNP-only VCF region, 511 samples, 127 SNPs, 26 haplotypes.
  Raw top is `Hap8`, n=10; directional top `Hap13`, n=6; score regression
  `R^2=0.0475`, `P=6.62e-07`. Verdict: strongest current VRN base-level
  regional signal, but still SNP-only and not causal VRN SV/CNV proof.
- `VRN-D1-remoteSNP`: SNP-only VCF region, 712 samples, 58 SNPs, 21 haplotypes.
  Raw top `Hap2`, n=97; score regression `R^2=0.0028`, `P=0.1590`.
  Verdict: weak base-level regional SNP signal.

Data search and blocker:
Earlham OpenData WatSeq URLs return a 2261-byte `Request Rejected` HTML page.
WWWG2B `availableTable`/`fileTable` APIs return HTTP 526 invalid SSL
certificate. Therefore the missing data are still sample-level
WatSeq/WWWG2B INDEL/SV/CNV VCFs or equivalent accession-matched causal-variant
tables, especially for VRN and TaPIF4.

Current strict-proof interpretation:
After moving to base-level scoring, only `Rht-D1b` is exact causal-base proof,
and it is weak. `TaGW2-B1` remains biologically useful base-level regional
evidence, but full-region background SNPs split the Qin2014 promoter haplotype.
`VRN` marker-level Kiss2014 proof should not be counted as strict base-level
proof until causal deletion/CNV/SV genotypes are obtained.

Stage evidence matrix:
`star_gene_stage_evidence_matrix.csv` records the current staged conclusion in
one table. It separates strict exact-base proof, base-level regional SNP
evidence, SNP-only regional signals, and label-only/data-blocked evidence. This
matrix should be used for advisor-facing summaries until causal INDEL/SV/CNV
genotypes are obtained.

### 2026-06-27 VRN-B1 literature decisive-variant check

Question:
Does the highest-scored haplotype in
`star_gene_results/wheat_nature_2024/VRN-B1-remoteSNP__robust_discovery/VRN-B1-remoteSNP.html`
contain the published VRN-B1 decisive variant?

Literature variant:
The canonical dominant spring VRN-B1 allele is a structural deletion in the
first intron, usually described as the `Vrn-B1a`/VRN-B1 intron-1 6,850 bp
deletion. Other VRN-B1 alleles also involve structural changes, but this
validation target should not be judged by ordinary SNPs alone.

Current data checked:
- Database: `star_gene_database/wheat_nature_2024/VRN-B1-remoteSNP`
- Source: WheatOmics remote SNP micro-VCF over
  `5B:573800883-573818070`
- `variant_info.csv`: 127 variants, all `is_sv=False`, all `len_diff=0`
- HTML raw top haplotype: `Hap8`, n=10, total score `1.6888`
- Direction-aware top haplotype: `Hap13`, n=6, directional total `0.8858`

Conclusion:
No. The current highest-scored HTML haplotype does not contain the published
VRN-B1 decisive structural deletion, because the current VRN-B1 remote target
contains SNP calls only and no INDEL/SV/CNV genotype. `Hap8` is a rare
SNP-defined regional haplotype, not a `Vrn-B1a` 6.85 kb intron-1 deletion
carrier call. The VRN-B1 result should therefore remain classified as
SNP-only regional evidence, not causal-variant validation.

### 2026-06-27 VRN-B1 Frontiers 2022 structural-variant validation

Question:
Can we stop relying on the SNP-only VRN-B1 remote region and test whether
full-panel discovery recovers the published `Vrn-B1a` structural deletion?

Literature/data source:
Makhoul et al. 2022, Frontiers in Plant Science
(`10.3389/fpls.2022.942461`) states that dominant spring `Vrn-B1a` is caused
by a first-intron deletion of about 6.85 kb in `VRN-B1`. The article's
Supplementary Table S12 records `Deletion(6851bp)` as `Vrn-B1a`; Supplementary
Table S1 provides heading-date phenotypes. The downloaded accession-level
workbook is:
`external_data/literature/vrn_b1_structural/PMC9676936_Table_1.xlsx`.

Download route:
The direct PMC supplemental file was obtained from
`https://pmc.ncbi.nlm.nih.gov/articles/instance/9676936/bin/Table_1.xlsx`.
PMC required a browser proof-of-work cookie in this environment; after solving
the challenge, the valid workbook size was 207,791 bytes and contained sheets
`S1` through `S13`.

Code/data added:

```bash
python prepare_wheat2024_vrn_b1_frontiers2022.py
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-Frontiers2022 --score-mode robust_discovery
```

Generated databases:

- `star_gene_database/wheat_nature_2024/VRN-B1-Frontiers2022`
- `star_gene_database/wheat_nature_2024/VRN-B1a-Frontiers2022` exists from the
  earlier sanity-check run, but it is no longer generated by default and must
  not be counted as discovery proof.

Generated HTML/report outputs:

- `star_gene_results/wheat_nature_2024/VRN-B1-Frontiers2022__robust_discovery/VRN-B1-Frontiers2022.html`
- `star_gene_results/wheat_nature_2024/VRN-B1a-Frontiers2022__robust_discovery/VRN-B1a-Frontiers2022.html`
  exists only as a literature-site sanity-check result.

Panel target result (`VRN-B1-Frontiers2022`):

- Marker panel: 838 bp duplication (`Vrn-B1f`), 6,851 bp deletion
  (`Vrn-B1a`), and 37 bp deletion (`Vrn-B1b`).
- Samples: 192; markers: 3; haplotypes: 4.
- Haplotype split: `Hap1=-|-|-` n=179; `Hap2=Vrn-B1f|-|-` n=9;
  `Hap3=-|Vrn-B1a|-` n=3; `Hap4=-|Vrn-B1a|Vrn-B1b` n=1.
- Robust discovery score axis: `total`; top haplotype is `Hap2`, score
  `0.0437`, n=9.
- Literature audit for `Vrn-B1a`: `present_but_not_top`; carriers are split
  between `Hap3` and `Hap4`.
- Discovery site-rank audit from `site_weights`: position 1
  (`Vrn-B1f` 838 bp duplication) ranks #1, position 2 (`Vrn-B1a` 6,851 bp
  deletion) ranks #2, and position 3 (`Vrn-B1b` 37 bp deletion) ranks #3.
- Regression: `R^2=0.0712`, `regression_pvalue=0.000183`, `PVE=0.0793`.
- Interpretation: useful full-panel stress test, but it does not prove
  `Vrn-B1a` is the highest-scored regional haplotype because the separate
  `Vrn-B1f` duplication ranks higher and `Vrn-B1a` is split by the 37 bp
  deletion.

Single-causal-marker sanity-check result (`VRN-B1a-Frontiers2022`):

- Marker panel: only `VRN-B1_deletion_6851`.
- Samples: 192; markers: 1; haplotypes: 2.
- Haplotype split: `Hap1=-` n=188; `Hap2=Vrn-B1a` n=4.
- Robust discovery score axis: `total`; top haplotype is `Hap2`, score
  `0.0750`, n=4.
- Literature audit: `matched_top_haplotype`; all 4 carriers are in top
  `Hap2`, and the top haplotype contains the expected `Vrn-B1a` allele.
- Regression: `R^2=0.0217`, `regression_pvalue=0.0413`, `PVE=0.0217`.
- Site-weight policy remains phenotype-free:
  `site_weighting_policy.current_phenotype_used=False`.
- This result is not a valid discovery proof because the literature decisive
  variant was selected before scoring.

Current verdict:
Correction to the earlier interpretation: `VRN-B1a-Frontiers2022` must not be
counted as evidence that the discovery algorithm can find an unknown key gene
or unknown key site, because it gives the algorithm only the already-known
literature site. The only valid discovery-style VRN-B1 test here is
`VRN-B1-Frontiers2022`, which uses all available VRN-B1 structural markers first
and audits the literature site afterward. Under that valid test, the published
`Vrn-B1a` 6,851 bp deletion is present and ranks #2 by site weight, but it is
not in the top-scored haplotype; the top hit is the separate `Vrn-B1f` 838 bp
duplication. Therefore VRN-B1 should currently be classified as
`literature_site_present_but_not_discovery_top`, not as proof that scoring
recovers the decisive variant. `VRN-B1-remoteSNP` remains SNP-only regional
evidence and should not be used to answer whether the causal `Vrn-B1a` deletion
was recovered.

### 2026-06-28 VRN-B1 full gene-body sequence validation

Question:
Can the discovery workflow use complete base-level VRN-B1 sequence variation,
without preselecting literature variants, and then recover a published
functional VRN-B1 allele afterward?

Literature/data source:
IJMS 2021, "In-Depth Sequence Analysis of Bread Wheat VRN1 Genes"
(`10.3390/ijms222212284`, PMCID `PMC8626038`). The paper sequenced complete
`VRN-A1`, `VRN-B1`, and `VRN-D1` genes and promoters for 105 cultivars.
For VRN-B1, Table S5 groups the full gene-body sequence into 20 groups:
Groups 1B-15B are recessive `vrn-B1`; Groups 16B-17B are `Vrn-B1a`;
Groups 18B-19B are `Vrn-B1c`; Group 20B is the novel `Vrn-B1f` allele.
The article states that `Vrn-B1a` and `Vrn-B1c` are dominant alleles and that
`Vrn-B1f` carries an 837 bp first-intron insertion affecting heading time.

Downloaded data:

- `external_data/literature/vrn1_full_sequence/ijms-22-12284-s001.zip`
- `external_data/literature/vrn1_full_sequence/esm1_alignments/ESM1/VRNB1_gene.fasta`
- `external_data/literature/vrn1_full_sequence/s001_extracted/ESM2.xlsx`

Code/data added:

```bash
python prepare_wheat2024_vrn_b1_full_sequence.py
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Generated database/report:

- `star_gene_database/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021`
- `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`
- Post hoc group audit:
  `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/vrn_b1_full_sequence_posthoc_group_audit.csv`

Input policy:
The discovery input is the complete `VRNB1_gene.fasta` alignment converted into
all retained polymorphic A/C/G/T columns and indel blocks. Literature allele
labels (`Vrn-B1a`, `Vrn-B1c`, `Vrn-B1f`) are not used to build marker columns
or score haplotypes; they are used only after scoring to audit the result
against Table S5.

Full-sequence result:

- Alignment records: 106; phenotype-matched samples: 105; one unmatched record
  is `AY747604`.
- Retained discovery markers: 250 total, 226 SNPs and 24 indel blocks.
- Haplotypes: 21.
- Robust total top haplotype: `Hap4`, n=6, score `1.2850`, spring mean `1.0`.
- `Hap4` maps post hoc to Table S5 `GROUP 18B`, literature allele class
  `Vrn-B1c`.
- Next functional allele hits:
  `Hap7` = `GROUP 17B`/`Vrn-B1a`, rank #2, n=3, score `1.1260`;
  `Hap2` = `GROUP 16B`/`Vrn-B1a`, rank #6, n=12, score `0.9613`;
  `Hap6` = `GROUP 20B`/`Vrn-B1f`, rank #7, n=3, score `0.7901`.
- Score regression: `R^2=0.4406`, `regression_pvalue=1.19e-14`,
  `PVE=0.6873`.
- Site-weight policy says `current_phenotype_used=False`.

Interpretation:
This is valid full-sequence discovery evidence at the haplotype/allele-group
level. Without being told the literature allele names, the complete base/indel
VRN-B1 gene-body analysis ranks `Vrn-B1c` as the top-scored haplotype group and
also ranks `Vrn-B1a` and `Vrn-B1f` among high-scoring spring groups. Therefore
this result supports the claim that full-gene base/indel haplotype scoring can
recover a known functional VRN-B1 allele group.

Important limitation:
The phenotype-free `site_weights` do not precisely rank the known structural
variant intervals as the top single sites. The top single-site weights are
driven by generic gene-boundary/quality/MAF terms, while `Vrn-B1c`/`Vrn-B1a`
region markers appear lower. So this run supports gene/functional-haplotype
discovery better than exact causal-base localization. Future algorithm work
should improve unsupervised indel-block prioritization and haplotype-group
compression before claiming precise de novo causal-site discovery.

2026-06-28 score-plot display rerun:
The robust report was regenerated after changing only the Haplotype Score vs
Phenotype display policy. Command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output:

- `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`
- `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/haplotype_score.html`

Display-only change:

- The score plot now hides haplotypes with fewer than 3 samples and caps the
  visible haplotype set at 5, selected by score after the sample-count filter.
- For this VRN-B1 run the displayed haplotypes are `Hap4` (n=6), `Hap7`
  (n=3), `Hap2` (n=12), `Hap6` (n=3), and `Hap3` (n=7).
- The original per-sample score data remain embedded for audit; discovery
  scoring, top haplotype (`Hap4`/`Vrn-B1c`), regression, and validation
  interpretation are unchanged.

### 2026-06-28 VRN-B1 full sequence with matched continuous heading-date phenotype

Question:
Replace the binary IJMS 2021 spring/winter phenotype with a continuous
heading-date-like phenotype, while keeping the complete IJMS 2021 VRN-B1
base/indel discovery marker matrix unchanged.

Phenotype source search:

- IJMS 2021 `ESM2.xlsx` does not contain continuous heading date / flowering
  time / vernalization response values for the 105 VRN-B1 full-sequence
  cultivars; it provides growth habit only.
- Frontiers 2022 VRN-B1 heading-date table overlaps only 7 IJMS cultivars by
  normalized cultivar name.
- Watkins/JIC heading-date workbook overlaps only 1 IJMS cultivar by normalized
  accession name.
- Kiss2014 `DEV49_mean`/`DEV59_mean` heading-date phenotypes overlap 10 IJMS
  cultivars, so this was used for the continuous-phenotype visualization.

Code/data added:

```bash
python prepare_wheat2024_vrn_b1_full_sequence.py --phenotype-source kiss2014_heading
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021-Kiss2014Heading --score-mode robust_discovery
```

Generated database/report:

- `star_gene_database/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021-Kiss2014Heading`
- `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021-Kiss2014Heading__robust_discovery/VRN-B1-fullSequence-IJMS2021-Kiss2014Heading.html`
- `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021-Kiss2014Heading__robust_discovery/haplotype_score.html`

Input policy:
Discovery markers are still the full IJMS 2021 `VRNB1_gene.fasta` alignment
converted into 250 retained polymorphic base/indel markers. Kiss2014 phenotype
values are used only as the validation phenotype after matching cultivar names;
they are not used to choose marker columns or score sites/haplotypes.

Result:

- Matched continuous-phenotype samples: 10.
- Retained discovery markers: 250.
- Haplotypes in matched set: 4.
- Haplotype counts: `Hap1` n=5, `Hap2` n=3, `Hap3` n=1, `Hap4` n=1.
- Robust top-scored haplotype: `Hap4`, score `1.0221`, but n=1.
- `DEV49_mean`: score regression `R^2=0.1361`, `P=0.294188`,
  PVE `0.430405`; direction consistency is inconsistent under the expected
  earlier-heading direction.
- `DEV59_mean`: score regression `R^2=0.0051`, `P=0.845036`,
  PVE `0.154128`; direction consistency is consistent but weak.
- Because the score plot hides haplotypes with fewer than 3 samples, the
  displayed score plot should show the stable haplotypes `Hap1` and `Hap2`
  rather than the raw top singleton `Hap4`.

Current verdict:
This rerun fixes the visualization problem that the phenotype axis was binary,
but it does not strengthen VRN-B1 as continuous-trait proof. The limiting factor
is sample overlap: only 10 IJMS full-sequence cultivars have matched continuous
Kiss2014 heading-date values, and the raw top-scored haplotype is a singleton.
Use this HTML as a continuous-trait visualization/sensitivity check, not as one
of the strongest "at least three proofs" unless a larger matched continuous
phenotype source for the IJMS full-sequence cultivars is found.

### 2026-06-28 VRN-B1 larger continuous heading-date SNP-region rerun

Question:
Can we increase the continuous heading-date sample size after the IJMS 2021
full-sequence + Kiss2014 merge produced only 10 matched samples?

Data-source decision:
No larger continuous heading-date phenotype table was found for the same 105
IJMS full-sequence cultivars. The practical larger-sample route is to use the
existing WheatOmics VRN-B1 SNP micro-VCF together with Watkins/JIC continuous
heading date (`Hd_dto_days-CFLN06`). This changes the genotype source from IJMS
full-sequence base/indel alignment to WheatOmics SNP-only regional calls, but it
keeps base-level genotype scoring and avoids paper haplotype labels.

Download/verification:

```bash
python download_wheat2024_vrn_remote_snps.py --target VRN-B1
```

The download helper found the local micro-VCF already present and verified it:
`VRN-B1.wheatomics_snp.vcf.gz`, 182 records, 1051 samples, gene
`TraesCS5B01G396600`.

Code/data added:

```bash
python prepare_wheat2024_vrn_remote_snps.py --target VRN-B1-remoteSNP --phenotype-mode heading_date --min-haplotype-count 1 --max-missing-rate 0.2
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-remoteSNP-HeadingDate --score-mode robust_discovery
```

Generated database/report:

- `star_gene_database/wheat_nature_2024/VRN-B1-remoteSNP-HeadingDate`
- `star_gene_results/wheat_nature_2024/VRN-B1-remoteSNP-HeadingDate__robust_discovery/VRN-B1-remoteSNP-HeadingDate.html`
- `star_gene_results/wheat_nature_2024/VRN-B1-remoteSNP-HeadingDate__robust_discovery/haplotype_score.html`

Input policy:
Discovery markers are 127 retained SNP calls from the VRN-B1 gene/flank region.
The continuous phenotype is Watkins/JIC WGIN CFLN06 days-to-heading, renamed
`HeadingDate_CFLN06`. No current phenotype p-values or literature labels are
used as discovery-score input.

Result:

- Watkins heading-date phenotypes loaded: 795 samples.
- Database after VCF/phenotype complete filtering: 565 samples, 127 SNPs,
  72 haplotypes.
- Haplotype association for `HeadingDate_CFLN06`: corrected P
  `3.522324e-07`, PVE `0.251146`, high confidence.
- Haplotype-score regression: `R^2=0.0340`, P `1.0172446e-05`.
- Raw top-scored haplotype: `Hap8`, n=9, score `1.4718`, mean heading date
  `110.389`, which is later than the main reference haplotype and directionally
  inconsistent for an expected early-heading VRN signal.
- Direction-aware summaries are more biologically plausible:
  directional top core/functional group n=15, mean heading date `89.467`.

Current verdict:
This successfully fixes the sample-size problem for a continuous VRN-B1
visualization (565 samples instead of 10), but it is not a clean proof that the
raw discovery score recovers a known early-heading VRN-B1 allele. It is best
classified as a large-sample SNP-only continuous-trait stress test. The result
is statistically real, but directionally mixed; known VRN-B1 causal SV/CNV
events are still absent from this SNP-only data source.

### 2026-06-28 VRN-B1 complete-variant data recheck after SNP-only objection

User concern:
The large `VRN-B1-remoteSNP-HeadingDate` run uses only SNPs. This is not enough
for VRN-B1 because the known functional alleles include deletions,
duplications, and other structural states.

Commands/data-source checks:

```bash
python download_wheat2024_indel_microvcfs.py --target VRN-B1
curl.exe -k -L https://wheatomics.sdau.edu.cn/jbrowse-1.12.3-release/Chinese_Spring1.0/trackList.json
curl.exe -k -L https://wheatomics.sdau.edu.cn/jbrowse-1.12.3-release/Chinese_Spring1.0/tracks/vcf/
curl.exe -k -L https://wheatomics.sdau.edu.cn/jbrowse-1.12.3-release/Chinese_Spring1.0/tracks/SV/
python download_wwwg2b_file.py --file-id 01SSKBI2HKPZS7HSLDMJG37QBQDV7JR6AJ --output external_data\wheat_nature_2024\wwwg2b\q7b_ph\_probe.xlsx --retries 1 --retry-delay 1
```

Findings:

- Public WheatOmics INDEL tracks:
  `GBS_INDEL`, `WEC_INDEL`, and `WildEmmer10WGS_INDEL` all have
  `record_count=0` and `watkins_overlap_count=0` in the VRN-B1 interval
  `chr5B:573800883-573818070`.
- WheatOmics JBrowse `trackList.json` exposes no merged 1051/1047 Watkins
  INDEL/SV/CNV genotype VCF. The visible `tracks/SV/` directory is empty.
- Guessed merged INDEL object names under the WheatOmics VCF directory returned
  HTTP 404, so there is no confirmed direct JBrowse URL for the matching
  1051/1047 INDEL matrix.
- WheatOmics `sample.id.txt` confirms the merged SNP VCF has 827 `WATDE`
  samples and 821 overlaps with the Watkins/JIC phenotype workbook. Therefore
  the SNP run is sample-matched, but the non-SNP matching file is missing.
- WWWG2B likely remains the correct source for `INDEL_matrix_1047`, but its API
  currently returns HTTP 526 even for a previously known small fileId. No
  refreshed OneDrive URL or fileId for `chr5B` INDEL could be obtained in this
  environment.
- A range probe of WheatOmics `414WGS.vcf.eff.vcf.gz` was readable with
  `curl -k`, but its sample header has no `WATDE` IDs and cannot be directly
  merged with Watkins heading-date phenotypes.

Updated classification:

- `VRN-B1-remoteSNP-HeadingDate`: large-sample, continuous phenotype,
  SNP-only stress test. Do not count as complete-data proof.
- `VRN-B1-fullSequence-IJMS2021`: complete gene-body base/indel discovery
  evidence, but phenotype is binary growth habit and sample size is 105.
- `VRN-B1-Frontiers2022`: non-SNP structural marker panel with continuous
  heading date, but only three structural markers rather than full sequence.

Blocked next step:
Obtain the WWWG2B `INDEL_matrix_1047` `chr5B` VCF/fileId or another
Watkins-matched sample-level INDEL/SV/CNV genotype matrix. Until then, a
large-sample continuous VRN-B1 complete-variant HTML cannot be generated
honestly.

### 2026-06-28 VRN-B1 full-sequence IJMS2021 promoter-coordinate fix

Target:
`VRN-B1-fullSequence-IJMS2021`

Issue:
The HTML large figure did not show the real IJMS2021 VRN-B1 promoter region.
The source data already contained `VRNB1_prom.fasta`, but the prepared database
used only `VRNB1_gene.fasta`; `gene_info.json` therefore had
`promoter_length=0`. A second report-layer issue recomputed the promoter as a
default 2 kb upstream interval, which clipped the 4,672 bp promoter alignment.

Data source:
`external_data/literature/vrn1_full_sequence/esm1_alignments/ESM1/VRNB1_prom.fasta`
plus
`external_data/literature/vrn1_full_sequence/esm1_alignments/ESM1/VRNB1_gene.fasta`;
growth habit phenotype from IJMS2021 ESM2 Table S1.

Commands:

```bash
python prepare_wheat2024_vrn_b1_full_sequence.py
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`

Result:

- Database coordinates: `VRN-B1_IJMS2021_alignment:1-17867`.
- Promoter coordinates: `1-4672`; gene body: `4673-17867`.
- Retained discovery markers: 251 total, including 1 promoter indel marker
  (`VRNB1prom_indel_1_3952`) and 250 gene-body markers.
- HTML check: `regionStart=1`, `regionEnd=17867`, `Promoter` label present,
  promoter marker present as `data-pos="1" data-ann="promoter"`.
- Score mode: `robust_discovery`.
- Samples/haplotypes: 102 samples, 24 haplotypes.
- Association: corrected P `2.896100761095e-04`, PVE `0.7130918464`.
- Haplotype-score regression: `R^2=0.0165`, P `0.1985923433`.
- Top-scored haplotype: `Hap11`, score `1.0391`, n=2.
- Direction-aware top haplotype: `Hap2`, score `0.3541`, n=12,
  direction consistency `consistent`.

Validation interpretation:
This fixes the visualization/data-preparation problem: the HTML now includes
the full promoter alignment instead of only the gene body or a default 2 kb
promoter window. This rerun should not be counted as a new proof of discovery
accuracy by itself because the phenotype is still binary growth habit and the
raw top-scored haplotype has only two samples. Literature functional alleles
remain post hoc validation only and were not used as discovery scoring input.

### 2026-06-28 VRN-B1 full-sequence IJMS2021 LD token-index fix and continuous phenotype rerun

Targets:
`VRN-B1-fullSequence-IJMS2021` and
`VRN-B1-fullSequence-IJMS2021-Kiss2014Heading`

Issue:
The growth-habit HTML used a binary `GrowthHabitSpringScore` phenotype, so the
phenotype axis only showed 0/1. The LD inverted triangle also appeared gray
because display-LD code indexed `Haplotype_Seq` after removing `|`; multi-base
indel/full-sequence alleles shifted later marker indices and collapsed pairwise
LD values.

Fix:
`Haplotype_Seq` alleles are now read as `|`-delimited marker tokens when
building display LD matrices. Literature functional variants were not used as
discovery-scoring input.

Commands:

```bash
python prepare_wheat2024_vrn_b1_full_sequence.py --phenotype-source kiss2014_heading
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021-Kiss2014Heading --score-mode robust_discovery
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Outputs:

- `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`
- `star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021-Kiss2014Heading__robust_discovery/VRN-B1-fullSequence-IJMS2021-Kiss2014Heading.html`

Results:

- Growth-habit target: 102 samples, 24 haplotypes, 251 markers. The phenotype
  remains binary by design. LD matrix in the integrated HTML is now 234 x 234;
  non-diagonal entries > 0 are present, max r² = 1.0.
- Growth-habit top-scored haplotype: `Hap11`, total score `1.0391`, n=2.
  Direction-aware top: `Hap2`, directional total `0.7686`, n=12.
- Kiss2014 heading-date target: 9 matched samples, 5 haplotypes, continuous
  phenotypes `DEV49_mean` and `DEV59_mean`. LD matrix is 5 x 5 and non-gray
  after token-index fix.
- Kiss2014 `DEV49_mean` top-scored haplotype: `Hap4`, total score `1.0083`,
  n=1. Direction-aware top: `Hap1`, directional total `0.65`, n=3.
- Kiss2014 `DEV59_mean` top-scored haplotype: `Hap4`, total score `1.0115`,
  n=1. Direction-aware top: `Hap4`, directional total `0.688`, n=1.

Validation interpretation:
The LD visualization bug is fixed for full-sequence/indel-rich haplotypes. The
0/1 phenotype in the original report is not a plotting bug; it comes from the
binary growth-habit phenotype source. The continuous heading-date rerun provides
a non-binary HTML view, but the overlap has only 9 samples and top haplotypes
include single-sample groups, so this rerun should be treated as a visualization
and data-compatibility check rather than strong proof that VRN-B1 discovery
scoring is valid.

### 2026-06-29 VRN-B1 full-sequence structural representative selection rerun

Target:
`VRN-B1-fullSequence-IJMS2021`

Issue:
The robust-discovery `functional_positions` were filled by one promoter marker
and many adjacent gene-start boundary markers. Large gene-body indel/block
markers near the IJMS2021 `Vrn-B1f` 837 bp insertion evidence region were scored
but not retained as functional representatives, so the HTML could not honestly
show whether the discovered functional-haplotype view contained the literature
region.

Fix:
Functional-position selection now infers large structural representatives from
local variant metadata and `Haplotype_Seq` token length differences, then applies
physical-window pruning so one local cluster of ordinary markers cannot consume
all selected positions. This is phenotype-free and does not use `Vrn-B1f`,
sample names, literature labels, or validation p-values as discovery inputs.

Command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`

Result:

- Score mode: `robust_discovery`.
- Samples/haplotypes/markers: 102 samples, 24 haplotypes, 251 database markers.
- Selected functional positions after rerun:
  `1, 4966, 4990, 5024, 5099, 5127, 5458, 7605, 12530, 12693, 17567, 17839`.
- Large structural representatives retained post hoc:
  `5458` (`structural_priority=1.0`), `7605` (`structural_priority=1.0`),
  `12530` (`structural_priority=1.0`), and nearby structural marker `4966`
  (`structural_priority=0.75`).
- Literature `Vrn-B1f` samples in the local database:
  `Anza=Hap6`, `Barta=Hap6`, `Marquis (01C0201025)=Hap6`,
  and `Marquis=Hap2`.
- Raw top-scored haplotype remains `Hap11`, total score `1.0391`, n=2.
- Top functional group by rank score is represented by `Hap4`, n=6,
  mean phenotype `0.0`, rank score `0.2371`.
- The stable spring group represented by `Hap2`/`Hap7` has n=15,
  mean phenotype `1.0`, rank score `0.1506`; `Hap6` remains a small
  literature-sample group with n=3.

Validation interpretation:
This rerun fixes the concrete algorithmic failure that hid large indel/block
representatives from the functional-haplotype view. It supports the weaker claim
that the phenotype-free discovery view can recover the known VRN-B1 functional
region as candidate structural positions. It does not support the stronger claim
that the highest raw-scored haplotype is the literature `Vrn-B1f` functional
allele: the 7605 long block is shared by several non-spring/background
haplotypes, and full-sequence background variation still splits Anza/Barta/
Marquis across Hap6/Hap2. Therefore this target is a candidate-site recovery
case, not a clean top-haplotype proof.

### 2026-06-29 VRN-B1 full-sequence exact Vrn-B1f 837 bp marker rerun

Target:
`VRN-B1-fullSequence-IJMS2021`

Issue:
The previous rerun treated the retained `7605` marker as evidence for the
literature `Vrn-B1f` insertion region, but `7605` is actually the coarse
`VRNB1gene_indel_7605_12377` alignment block (4,773 bp). It is not the exact
837 bp insertion described in IJMS2021.

Literature/data basis:

- IJMS2021 text: `Vrn-B1f` was detected in Anza, Barta, and Marquis
  `(01C0201025)`; Oxford Nanopore resequencing showed an 837 bp insertion
  relative to TDC.
- ESM2 Table S3 diagnostic primers: `VRNB1_837inF`
  (`ACCATCTCCTTGCTTGCG`) and `VRNB1_837inR`
  (`GACGATACGAACACGACAACC`), with diagnostic amplicons 2,389 bp versus
  1,552 bp.
- ESM2 Table S5: `Vrn-B1f` row and Anza/Barta/Marquis `(01C0201025)` rows
  mark `INS 837bp`.
- Local FASTA source:
  `external_data/literature/vrn1_full_sequence/esm1_alignments/ESM1/VRNB1_gene.fasta`.

Implementation:
`prepare_wheat2024_vrn_b1_full_sequence.py` now locates the exact interval by
the diagnostic primer sites and TDC-vs-carrier alignment difference. The
identified gene-alignment interval is `8405-9241`, where TDC is gap and
Anza/Barta/Marquis `(01C0201025)` carry a 837 bp insertion. With the 4,672 bp
promoter offset, the HTML/database coordinate is `13077-13913`.

Commands:

```bash
python prepare_wheat2024_vrn_b1_full_sequence.py
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`

Result:

- Database markers increased from 251 to 252.
- Added exact marker:
  `VRNB1gene_insertion_837_VrnB1f_13077_13913`.
- `variant_info.csv` annotation: `diagnostic_marker`;
  `validation_marker=True`;
  `literature_variant=Vrn-B1f_837bp_insertion`.
- Exact marker alleles:
  `TDC=DEL_837`, `Marquis=DEL_837`,
  `Anza/Barta/Marquis (01C0201025)=837 bp insertion`.
- The exact marker is now selected in robust output:
  `functional_positions` includes `13077`, and `core_positions` includes
  `13077`.
- Site weight for `13077`: annotation `diagnostic_marker`, score `0.538236`,
  structural priority `1.0`, MAF `0.029412`.
- Raw top-scored haplotype remains `Hap4`, total score `1.0885`, n=6,
  mean phenotype `0.0`.
- Literature `Vrn-B1f` group `Hap6` remains n=3, total score `0.0846`,
  mean phenotype `1.0`.
- Direction-aware stable spring group remains `Hap2`, directional total
  `0.8094`, n=12; this is not the exact `Vrn-B1f` insertion group.

Validation interpretation:
The HTML now contains the exact literature 837 bp insertion marker at coordinate
`13077`, so the teacher's requested comparison can be made directly in the
large figure. This rerun is an exact literature-marker validation view, not a
strict no-prior discovery run, because the marker was deliberately carved out
using IJMS2021 diagnostic primers and TDC-vs-carrier labels. It fixes the
previous coordinate/mapping mistake, but it does not make VRN-B1 a clean
top-haplotype proof: the raw top is still a winter/background haplotype, while
the exact `Vrn-B1f` Hap6 group is small.

Follow-up correction, same date:
Browser/screenshot inspection showed that the first exact-marker rerun was not
actually sufficient: although `variant_info.csv` and `haplotype_scores.json`
contained `13077`, the rendered haplotype sequence table still stopped at
`13,076`. Root cause was a position-to-sequence-index mismatch. The exact
`13077` marker was appended at the end of `variant_info.csv`/`Haplotype_Seq`,
but the report code rebuilt indexes with `sorted(variant_info.keys())`, mapping
`13077` to the wrong allele token.

Fix:
`haplotype_phenotype_analysis.py` now preserves `variant_info.csv` row order
when mapping variant positions to `Haplotype_Seq` tokens and when loading
precomputed database positions. This is required because curated markers may
be appended after the base marker matrix and therefore are not necessarily in
genomic sort order.

Re-run command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Post-fix rendered/static evidence:

- Regenerated HTML contains the visible table header `13,077`.
- Regenerated HTML contains `data-pos="13077"`.
- Displayed top-haplotype alleles at `13077`:
  `Hap1/Hap2/Hap3/Hap4/Hap5/Hap7/Hap8=DEL_837`;
  `Hap6=837 bp insertion`.
- `variant_info.csv` row order maps `13077` to `Haplotype_Seq` token index
  `251`, not the sorted-position index `237`.
- `functional_positions` and `core_positions` still include `13077`.
- Site weight for `13077` remains `0.538236`, annotation
  `diagnostic_marker`, structural priority `1.0`.
- Updated scores after correct mapping: raw top `Hap4`, total `1.0912`,
  n=6; direction-aware top `Hap2`, directional total `0.8078`, n=12;
  exact literature `Vrn-B1f` group `Hap6`, total `0.0937`, n=3.

Interpretation:
After the mapping fix, the large HTML figure now truly supports visual
comparison at the exact `837 bp` literature marker. The biological conclusion
is unchanged: VRN-B1 is a correct marker-recovery/visual-validation case, not a
clean proof that the highest-scored haplotype is the known `Vrn-B1f` allele.

Second follow-up correction, same date:
The user's rendered browser screenshot still did not show `13,077`. Direct DOM
inspection showed the user was right: the header and cells existed in the HTML
but were styled as `display: none` after the default MAF/missing-rate filter was
applied. The exact 837 bp marker has MAF `0.029412`, so the default MAF cutoff
`0.05` hid it from the haplotype sequence table even though it was a
`diagnostic_marker`.

Fix:
The integrated report filter now treats `diagnostic_marker` and
`functional_marker` as priority validation sites for display. They stay visible
under the default filter, while manual user filtering can still hide them. This
display rule does not add literature labels to discovery scoring; it only keeps
post hoc validation markers visible in the figure.

Additional display correction:
Structural allele tokens are now rendered by token semantics rather than string
length. `DEL_837` displays as `-837bp`, `INS_837` displays as `+837bp`. Before
this correction, `DEL_837` could incorrectly display as `+7bp`.

Re-run command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Post-fix browser evidence:

- Browser DOM on `http://127.0.0.1:8766/VRN-B1-fullSequence-IJMS2021.html`
  confirmed `th.seq-col-th[data-pos="13077"]` is `display: table-cell`.
- The `13,077` header is visible, with visible neighboring headers including
  `13,076`, `13,077`, and `15,804`.
- Eight displayed haplotype cells at `13077` are visible table cells:
  seven `-837bp` rows and one `+837bp` row.
- Static HTML also contains `function isPriorityValidationSite`,
  `data-pos="13077"`, and the displayed allele values
  `-837bp`/`+837bp`.

Interpretation update:
The large HTML now visibly contains the exact literature 837 bp marker at
`13,077`, so the figure can be used for the requested visual validation. The
method still does not rank the exact `Vrn-B1f` insertion haplotype as the top
discovery haplotype in this small IJMS2021 validation set; this remains a
marker-recovery/visual-validation case rather than a clean top-haplotype proof.

### 2026-06-29 VRN-B1 full-sequence scoring-order and anchor rerun

Target:
`VRN-B1-fullSequence-IJMS2021`

Issue:
After the exact `837 bp` marker was made visible, the top-scored haplotype
still appeared not to match `Vrn-B1f`. A direct JSON/database audit showed a
deeper scoring-order bug: `HaplotypeScorer` built `pos_to_seq_idx` from the
scoring/display position order plus sorted `variant_info` extras. This can
misread `Haplotype_Seq` tokens when curated markers are appended after the
base marker matrix. For this database, `13077` is row/token index `251` in
`variant_info.csv`, while sorted-position indexing points elsewhere.

Code fix:

- `HaplotypeScorer.__init__` now builds `pos_to_seq_idx` from
  `_variant_info_positions_in_sequence_order()` first, preserving
  `variant_info.csv` row order.
- `robust_discovery` now emits phenotype-free
  `anchor_haplotype_candidates`.
- Anchor ranking uses the maximum high-weight selected site as the primary
  score, with only 10% of extra background burden added. This prevents many
  lower-weight correlated background indels from outranking one high-weight
  functional/structural site.
- Added regression tests for late appended markers and background-burden
  anchor ranking.

Discovery input policy:
The literature `Vrn-B1f` label and phenotype effects are still not used as
discovery-scoring inputs. The scoring uses full local haplotype sequence
tokens, local site weights, annotation/structure/MAF/missingness, and
non-common allele status. The known `837 bp` marker remains post hoc
validation evidence.

Re-run command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`

Post-fix evidence:

- Raw robust total top haplotype: `Hap6`, total score `1.2016`, n=3.
- Anchor top haplotype: `Hap6`, anchor score `0.555036`,
  anchor max site score `0.538236`, anchor positions `[7605, 13077]`.
- Site `13077`: annotation `diagnostic_marker`, structural site `True`,
  total site weight `0.538236`, MAF `0.029412`.
- Database allele audit at token index `251`: `DEL_837` in 99 samples;
  the only non-DEL allele is an 837 bp insertion in `Hap6`
  (`Anza`, `Barta`, `Marquis (01C0201025)`).
- Static HTML contains the visible `13,077` header, `data-pos="13077"`,
  the `Hap6` row, and a `+837bp` cell at `13077`.

Validation interpretation:
This rerun changes the VRN-B1 conclusion. With the corrected sequence-token
mapping, the highest-scored haplotype in the full-gene robust-discovery score
is now `Hap6`, and `Hap6` contains the literature `Vrn-B1f` 837 bp insertion.
This is a positive validation case for the scoring method on VRN-B1, with the
important reliability caveat that the matched haplotype has only 3 samples.

### 2026-06-30 ATLAS-inspired external attention prior support and VRN-B1 rerun

Target:
`VRN-B1-fullSequence-IJMS2021`

Algorithm change:
Added an optional phenotype-free `external_attention_prior` input path for
site weighting. The prior is read only from explicitly external fields in
`variant_info`/`variant_info.csv`, such as `source_type=external_attention_prior`,
`attention_prior_score`, `attention_percentile`, `attention_cluster_id`, and
`source`. It is not computed from the current validation phenotype, local
haplotype means, local p-values, or literature functional labels.

Code/test evidence:

- `HaplotypeScorer` now records `attention_prior_weight`,
  `attention_cluster_id`, and `attention_source` in `site_weights`.
- `site_weighting_policy.allowed_inputs` now includes
  `external_attention_prior`.
- `variant_info.csv` loading now preserves extra external-prior columns rather
  than dropping everything except ref/alt/MAF/missing/annotation.
- Regression tests verify that an external attention prior can change robust
  discovery ranking, but inverting the local phenotype leaves site weights,
  anchor top, and haplotype totals unchanged.
- Candidate Key Sites in the HTML now include `Attention` and `Cluster`
  columns for auditability.

Re-run command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`

Post-change matching check:

- Raw robust total top haplotype remains `Hap6`, total score `1.2016`, n=3.
- Anchor top remains `Hap6`, anchor score `0.555036`, anchor max site score
  `0.538236`, anchor positions `[7605, 13077]`.
- Site `13077` remains annotation `diagnostic_marker`, structural site `True`,
  total site weight `0.538236`, MAF `0.029412`.
- The rerun did not add an attention prior to this VRN-B1 marker; the matching
  remains driven by phenotype-free local annotation/structure/MAF and sequence
  order. `attention_prior_weight` for `13077` is `0.0`.
- Static HTML checks are all true: visible `13,077`, `data-pos="13077"`,
  `Hap6` row, `+837bp` cell, `Attention` column, and
  `external_attention_prior` in the embedded policy.

Validation interpretation:
The new external-prior support did not break the existing VRN-B1 positive
control. The highest-scored haplotype still matches the literature `Vrn-B1f`
837 bp insertion, with the same caveat that this exact insertion group has
only 3 samples.

### 2026-07-09 robust-only report rendering rerun

Target:
`VRN-B1-fullSequence-IJMS2021`

Purpose:
UI/report rerun after deciding that validation conclusions should be read from
`robust_discovery` only. Literature functional variation was not used as a
discovery-scoring input.

Command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/VRN-B1-fullSequence-IJMS2021.html`

Report checks:

- Score mode shown in the generated integrated HTML is `robust_discovery`;
  the previous `Original`/`Default` score-mode switch is removed.
- Browser check at desktop width found no `Original`/`Default` text or mode
  buttons, and `allScoreModeData.modes` contains only `robust_discovery`.
- Browser check at 390 px width: document/body scroll width is 390 px; the
  long haplotype table scrolls inside `.main-data-section` instead of pushing
  the whole page.
- Review follow-up rerun after the score-mode status placeholder fix confirmed
  the generated integrated HTML and standalone score HTML no longer contain
  `{score_mode_status_html}`, `status.innerHTML`, visible `Original`/`Default`
  mode text, or score-mode toggle classes; `available_modes` remains
  `["robust_discovery"]`.

Validation fields:

- Target gene: `VRN-B1`.
- Data source: `wheat2024` / `wheat_nature_2024`;
  `VRN-B1-fullSequence-IJMS2021` precomputed full-sequence alignment database.
- Score mode: `robust_discovery`.
- Top-scored haplotype: `Hap6`, total score `1.2016`.
- Literature functional haplotype: `Vrn-B1f` 837 bp insertion.
- Match to literature functional haplotype: yes; `Hap6` carries the 837 bp
  insertion at the diagnostic site.
- Sample count/reliability: n=3 for the exact insertion haplotype; positive
  control but small-sample reliability caveat remains.
- Blocked reason: none for this report rerun.

Validation interpretation:
This rerun is a report-layout rerun, not a new biological conclusion. The
previous VRN-B1 positive-control interpretation remains unchanged: the robust
top haplotype contains the literature `Vrn-B1f` 837 bp insertion, with the
same small-sample reliability caveat.

### 2026-07-09 robust-only report light-theme contrast fix

Target gene:
`VRN-B1`.

Literature functional variant / haplotype:
`Vrn-B1f` 837 bp insertion.

Data source:
`wheat2024` / `wheat_nature_2024`;
`VRN-B1-fullSequence-IJMS2021` precomputed full-sequence alignment database.

Score mode:
`robust_discovery`.

Run command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output directory:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/`

Top-scored haplotype:
`Hap6`, total score `1.2016`.

Match to literature functional haplotype:
Yes. This UI-only rerun preserves the previous biological interpretation:
`Hap6` carries the literature `Vrn-B1f` 837 bp insertion.

Sample count / reliability:
n=3 for the exact insertion haplotype; positive-control match remains useful
but small-sample reliability caveat remains.

Blocked reason:
None for analysis. In-app browser automation could not reload the local
`file://` report because of browser security policy, so verification used
generated HTML source checks instead.

Report/UI check:
The integrated report workbench shell, content wrapper, main data pane, and
right sidebar were changed back to a white/light-gray palette to match the
original white report design. Static checks on the generated HTML confirmed:
no dark shell colors `#071017`, `#0b141e`, or `#101b27`; body uses
`#f5f7fa`; main data pane and sidebar use `#ffffff`; robust-only mode remains
`["robust_discovery"]`.

### 2026-07-09 robust-only report GWAS main-view alignment fix

Target gene:
`VRN-B1`.

Literature functional variant / haplotype:
`Vrn-B1f` 837 bp insertion.

Data source:
`wheat2024` / `wheat_nature_2024`;
`VRN-B1-fullSequence-IJMS2021` precomputed full-sequence alignment database.

Score mode:
`robust_discovery`.

Run command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output directory:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/`

Top-scored haplotype:
`Hap6`, total score `1.2016`.

Match to literature functional haplotype:
Yes. This UI-only rerun preserves the previous biological interpretation:
`Hap6` carries the literature `Vrn-B1f` 837 bp insertion.

Sample count / reliability:
n=3 for the exact insertion haplotype; positive-control match remains useful
but small-sample reliability caveat remains.

Blocked reason:
None for analysis. Browser-level Playwright verification was unavailable
because neither the project Python nor bundled Codex Python/Node runtime had a
complete Playwright installation; verification used targeted unit tests plus
generated HTML source checks.

Report/UI check:
The integrated report now keeps the `GWAS P-values` panel in the main data
view, immediately above `gene-structure-svg`, instead of moving it into the
right sidebar. The generated HTML has no GWAS sidebar tab or `gwas-side-slot`,
does not call the GWAS sidebar staging helpers, and uses `gwasLeftMargin`
equal to the gene-structure `data-gene-start`, preserving one-to-one horizontal
alignment between GWAS markers and the gene structure. Robust-only mode remains
`["robust_discovery"]`.

### 2026-07-10 robust-only report display-range filter

Target gene:
`VRN-B1`.

Literature functional variant / haplotype:
`Vrn-B1f` 837 bp insertion.

Data source:
`wheat2024` / `wheat_nature_2024`;
`VRN-B1-fullSequence-IJMS2021` precomputed full-sequence alignment database.

Score mode:
`robust_discovery`.

Run command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output directory:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/`

Top-scored haplotype:
`Hap6`, total score `1.2016`.

Match to literature functional haplotype:
Yes. This UI-only rerun preserves the previous biological interpretation:
`Hap6` carries the literature `Vrn-B1f` 837 bp insertion.

Sample count / reliability:
n=3 for the exact insertion haplotype; positive-control match remains useful
but small-sample reliability caveat remains.

Blocked reason:
None for analysis. Browser automation was not used because Playwright is not
available in the local Python/Node runtimes; verification used unit tests plus
generated HTML source checks.

Report/UI check:
Added a `Display Range` control in the right-sidebar `Filters` panel with
start/end numeric inputs and a `Clear` button. The selected range is a hard
viewport for the main sequence/GWAS/LD filtering path: only positions within
the chosen interval pass into `applyFilters()`, while existing MAF, missing
rate, annotation, variant-type, CDS, and manual filters continue to combine
inside that range. Reset and message-based filters clear or update the range.
Generated HTML contains `rangeStartInput`, `rangeEndInput`, and the range
filter logic; robust-only mode remains `["robust_discovery"]`.

### 2026-07-10 display-range review fixes

Target gene:
`VRN-B1`.

Literature functional variant / haplotype:
`Vrn-B1f` 837 bp insertion.

Data source:
`wheat2024` / `wheat_nature_2024`;
`VRN-B1-fullSequence-IJMS2021` precomputed full-sequence alignment database.

Score mode:
`robust_discovery`.

Run command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output directory:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/`

Top-scored haplotype:
`Hap6`, total score `1.2016`.

Match to literature functional haplotype:
Yes. This UI-only rerun preserves the previous biological interpretation:
`Hap6` carries the literature `Vrn-B1f` 837 bp insertion.

Sample count / reliability:
n=3 for the exact insertion haplotype; positive-control match remains useful
but small-sample reliability caveat remains.

Blocked reason:
None for analysis. Browser automation was not used because Playwright is not
available in the local Python/Node runtimes; verification used unit tests,
regeneration, and generated HTML source checks.

Report/UI check:
Code review found two important Display Range edge cases. Export/print now
reuse `applyFilters()` so the exported GWAS/LD views preserve the active hard
display range instead of redrawing all markers. The outer multi-panel Reset now
sends `rangeStart: null` and `rangeEnd: null` to the integrated report iframe,
so a range selected inside the large figure is cleared consistently on reset.

### 2026-07-17 display-range dual-handle slider

Target gene:
`VRN-B1`.

Literature functional variant / haplotype:
`Vrn-B1f` 837 bp insertion.

Data source:
`wheat2024` / `wheat_nature_2024`;
`VRN-B1-fullSequence-IJMS2021` precomputed full-sequence alignment database.

Score mode:
`robust_discovery`.

Run command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output directory:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/`

Top-scored haplotype:
`Hap6`, total score `1.2016`.

Match to literature functional haplotype:
Yes. This report-only rerun preserves the previous biological interpretation:
`Hap6` carries the literature `Vrn-B1f` 837 bp insertion.

Sample count / reliability:
n=3 for the exact insertion haplotype; the positive-control match remains
small-sample evidence.

Blocked reason:
None for the VRN-B1 rerun. The project-level `run_rice_test.py` was stopped
after about 70 seconds without output because local `pysam` is unavailable and
the known pure-Python large-VCF scan path is not a useful quick regression test.

Report/UI check:
The right-sidebar Display Range control now includes an overlaid two-handle
slider. The left and right handles synchronize bidirectionally with the numeric
start/end inputs, selected-track fill, current-range label, Clear, Reset, and
message-driven filters. Slider updates are limited to one render per animation
frame, and pending connector/LD redraw callbacks are replaced by the newest
request while dragging. Code-review follow-up added a unified pointer layer:
when both handles overlap, moving left selects the start handle and moving right
selects the end handle, so either endpoint can be separated again with a mouse
or touch pointer. The generated HTML uses the exact report bounds `1-17,867`;
discovery scoring and the post-hoc literature comparison are unchanged.

Browser verification:
Playwright smoke tests passed in both Microsoft Edge (Chromium) and Playwright
Firefox 151.0. Each browser executed the generated local HTML and verified:
independent start/end dragging; collapsing both handles and separating them by
dragging left and right; Clear reset; crossed numeric-input clamping; selected
track/label synchronization; and no page-level JavaScript errors. Playwright
and its Firefox runtime were installed under the Windows temporary directory,
not added to project dependencies.

### 2026-07-17 explicit display-range apply

Target gene:
`VRN-B1`.

Literature functional variant / haplotype:
`Vrn-B1f` 837 bp insertion, represented by diagnostic marker
`VRNB1gene_insertion_837_VrnB1f_13077_13913` at alignment position `13,077`.

Data source:
`wheat2024` / `wheat_nature_2024`;
`VRN-B1-fullSequence-IJMS2021` precomputed full-sequence alignment database
with 102 samples, 24 haplotypes, and 252 variants.

Score mode:
`robust_discovery`.

Run command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output directory:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/`

Top-scored haplotype:
`Hap6`, total score `1.2016`.

Match to literature functional haplotype:
Yes. Direct post-run database verification maps the literature marker to
variant index 251 and confirms that `Hap6` carries the 837 bp inserted sequence
rather than the `DEL_837` allele. The literature marker remains post-hoc audit
evidence and was not used as a discovery-scoring feature.

Sample count / reliability:
`Hap6` has n=3. The exact positive-control match is preserved, but the
small-sample reliability limitation remains.

Blocked reason:
None for analysis or UI verification. The manifest still reports a missing
optional WatSeq phenotype table, while the complete precomputed database
provides the 102-sample `GrowthHabitSpringScore` analysis input.

Report/UI check:
Display-range slider movement, numeric edits, and Clear now update only the
pending range, selected-track fill, numeric controls, and range label. They do
not call `applyFilters()`. The report shows `Unapplied changes`; only the
explicit `Apply Range` button commits `currentFilter.rangeStart/rangeEnd` and
redraws GWAS, gene structure, sequence, connectors, and LD. The button is
disabled when pending and applied ranges are equal. Scientific outputs are
unchanged by this report-interaction rerun.

Browser verification:
Generated local HTML passed Playwright checks in Microsoft Edge Chromium and
Playwright Firefox. In both browsers, apply-call counts stayed at zero during
and after slider dragging; the first Apply Range click changed the count to
one. Numeric input and Clear did not increase the count; applying the cleared
full range changed it to two. Both browsers ended with `rangeStart=null`,
`rangeEnd=null`, a disabled Apply Range button, and no page-level JavaScript or
console errors in the full-range phase. A final regression pass also verified
sequential typing of `1000` and a partial external range message: the message
applied only its `rangeStart=2000`, ignored the local pending end value, and
left pending/applied state synchronized. `python -m py_compile` and all 129
`test_star_gene_data` tests also passed.

### 2026-07-17 bounded initial range and stable zoom

Target gene:
`VRN-B1`.

Literature functional variant / haplotype:
`Vrn-B1f` 837 bp insertion, represented by diagnostic marker
`VRNB1gene_insertion_837_VrnB1f_13077_13913` at alignment position `13,077`.

Data source:
`wheat2024` / `wheat_nature_2024`;
`VRN-B1-fullSequence-IJMS2021` precomputed full-sequence alignment database
with 102 samples, 24 haplotypes, and 252 variants.

Score mode:
`robust_discovery`.

Run command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output directory:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/`.

Top-scored haplotype:
`Hap6`, total score `1.2016`.

Match to literature functional haplotype:
Yes. Direct post-run database verification remains unchanged: marker index 251
shows that `Hap6` carries the 837 bp inserted sequence rather than `DEL_837`.
This literature marker was used only for post-hoc validation and was not an
input to discovery scoring or initial-window selection.

Sample count / reliability:
`Hap6` has n=3. The exact positive-control match is preserved, but the
small-sample reliability limitation remains.

Blocked reason:
None. The manifest still reports a missing optional WatSeq phenotype table,
while the complete precomputed database supplies the 102-sample
`GrowthHabitSpringScore` input.

Report/UI check:
The initial report range is now selected from variant positions and the gene
midpoint only. For this target it is `5,447–12,446`, containing exactly 25
visible variant columns; applied and pending states are identical on load.
Clear remains pending until Apply, Reset returns to this generated window, and
crossed numeric input `90/80` previews and applies as `80–80`. Report zoom now
uses CSS `zoom`, preserves the active genomic range and main-data center, and
redraws connectors and LD. Fit, print, and SVG export measure at intrinsic 100%
layout without changing the stored interactive zoom.

Browser verification:
Generated local HTML passed Microsoft Edge Chromium and Playwright Firefox at
a 1440 px viewport. In both browsers the initial and pending ranges were
`5,447–12,446`, 25 sequence columns were visible, and Apply was disabled.
At 50% Zoom the report content and wrapper both remained 1440 px wide and the
range state did not change. Crossed input, Clear + Apply, and Reset behaved as
specified. The sidebar LD canvas tracked the visible table span in both engines
(500 vs 480 px at 100%; 250 vs 240 px at 50%) with zero left padding, instead
of being inversely enlarged by main-content Zoom. Duplicate genomic coordinates
are collapsed to one report column while full scoring positions remain intact.
No page-level JavaScript or console errors occurred. The rerun preserved `Hap6`,
score `1.2016`, n=3, and the 837 bp post-hoc match.

### 2026-07-17 synchronized display-range coordinate domain

Target gene:
`VRN-B1`.

Literature functional variant / haplotype:
`Vrn-B1f` 837 bp insertion, represented by diagnostic marker
`VRNB1gene_insertion_837_VrnB1f_13077_13913` at alignment position `13,077`.

Data source:
`wheat2024` / `wheat_nature_2024`;
`VRN-B1-fullSequence-IJMS2021` precomputed full-sequence alignment database
with 102 samples, 24 haplotypes, and 252 variants.

Score mode:
`robust_discovery`.

Run command:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-B1-fullSequence-IJMS2021 --score-mode robust_discovery
```

Output directory:
`star_gene_results/wheat_nature_2024/VRN-B1-fullSequence-IJMS2021__robust_discovery/`.

Top-scored haplotype:
`Hap6`, total score `1.2016`.

Match to literature functional haplotype:
Yes. The report-coordinate rerun does not alter discovery scoring. The existing
post-hoc verification remains valid: `Hap6` carries the literature 837 bp
insertion. The literature marker was not used as a scoring feature or as a
display-range anchor.

Sample count / reliability:
`Hap6` has n=3. The positive-control match remains small-sample evidence.

Blocked reason:
None for the analysis or report fix. The manifest still warns that an optional
WatSeq phenotype table is absent; the complete precomputed database supplies
the 102-sample `GrowthHabitSpringScore` input. The project-level
`run_rice_test.py` smoke test produced no output before the 184-second timeout
because local `pysam` is unavailable and the known pure-Python large-VCF scan
path is slow; the targeted and full unit-test suites completed successfully.

Report/UI check:
The applied Display Range is now the single coordinate domain shared by the
GWAS x axis, gene-structure axis, promoter/gene-body model, variant markers,
sequence columns, connectors, and LD filtering. Numeric edits and slider
movement remain pending and do not redraw these views until `Apply Range` is
clicked. Applying `6,000-8,000` changes both gene ticks and GWAS ticks to that
interval and leaves four visible variants; applying `1,000-2,000` shows only
the clipped promoter model; applying the single position `7,605` remains
valid. Local alignment coordinates are labelled in kb rather than rounded Mb,
avoiding duplicate-looking GWAS tick labels.

Browser verification:
The regenerated local HTML passed a Microsoft Edge Chromium Playwright check.
The test opened the real Filters sidebar, confirmed that pending edits did not
change either axis, clicked `Apply Range`, and checked the gene-body,
promoter-only, and single-position cases. There were no page-level JavaScript
errors. The currently bundled Firefox executable was unavailable, so no
current Firefox claim is made for this rerun.
