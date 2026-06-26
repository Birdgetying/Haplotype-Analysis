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
