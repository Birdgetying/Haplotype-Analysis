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

1. `TaGW2-B1-remoteSNP__robust_discovery`: strong top-rank proof. The
   exact Qin2014 `Hap-6B-1 = A|G|C` is covered, segregating, large
   (`n=371`), and top-ranked.
2. `VRN-D1-Kiss2014__robust_discovery`: stable diagnostic-marker proof.
   The top robust single-marker haplotype is `VRN-D1=1` (`Hap2`, n=38)
   and has significant DEV49/DEV59 heading-date score regressions.
3. `Rht-Zanke2014__robust_discovery`: direction-aware proof. The raw
   top score still picks the high-plant-height wild-type haplotype, but the
   stable low-phenotype directional top is `Hap1=B1a|D1b` (`n=214`), and
   the literature audit marks `Rht-D1b` as `matched_directional_top_haplotype`.

| Target | Trait | Score mode | Status | Key result | Conclusion |
|---|---|---|---|---|---|
| TaGW2-B1 / TaGW2-6B | TGW | default | tested | Top score was rare ambiguous `Hap5=C/A\|G\|C`, not exact Qin2014 `A\|G\|C` | Biological signal present, default ranking over-rewarded rare ambiguous haplotype |
| TaGW2-B1 / TaGW2-6B | TGW | robust_discovery | tested | Top score became exact Qin2014 `Hap1=A\|G\|C`, n=371 | Strongest current proof that robust scoring can recover a known functional haplotype without using literature labels in scoring |
| TaGW2-A1 | TGW | default | tested | Gene-level signal exists, but `SNP-494 A` is absent in current complete Watkins samples | Not an exact Hap5 validation |
| TaGW2-D1 | TGW | default | tested | Weak low-confidence signal with only 2 retained SNPs | Limited marker coverage; not strong proof |
| VRN-Kiss2014 | heading date | default | tested | Literature all-spring `1\|1\|1` was top-scored but n=2 | Useful audit check, weak statistical proof |
| VRN-Kiss2014 | heading date | robust_discovery | tested | Top score became `Hap2=0\|0\|1`, n=31; all-spring `1\|1\|1` was no longer top | Robust mode favors stable VRN-D1 signal over rare all-spring haplotype |
| VRN-D1-Kiss2014 | heading date | robust_discovery | tested | Single-marker `VRN-D1=1` top haplotype `Hap2`, n=38; DEV49 P `2.98e-07`, DEV59 P `7.49e-05` | Counts as stable marker-level proof |
| Rht-D1b | plant height | default | tested | Top-scored haplotype contained functional T allele but only 5 carriers | Weak positive control |
| Rht-D1b | plant height | robust_discovery | tested | Top score remained functional T haplotype `Hap2`, n=3 | Weak positive control; carrier count too small |
| Rht-Zanke2014 | plant height | robust_discovery + directional top | tested | Raw top is high-plant-height `Hap2=B1a\|Rht-D1a`; stable low-phenotype directional top is `Hap1=B1a\|D1b`, n=214 | Direction-aware proof and evidence that raw total score must report trait direction for decreases_trait targets |
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
| TaGW2-B1 / TaGW2-6B | Qin et al. 2014, BMC Plant Biology, `https://doi.org/10.1186/1471-2229-14-107` | Favorable promoter haplotype `Hap-6B-1`, diagnostic states `-1709/-721/-83 = A/G/C`; TaGW2 encodes a RING-type E3 ubiquitin ligase and the favored haplotype increases grain width/weight. | WWWG2B chr6B SNP VCF records `chr6B:291759689 C>A`, `chr6B:291760677 G>A`, `chr6B:291761315 T>C`; 816 Watkins TGW-overlap samples. | Covered and segregating; exact `A\|G\|C` haplotype is `Hap1`, n=371. | Default: no, top was rare ambiguous `Hap5=C/A\|G\|C`; robust_discovery: yes, top is exact `Hap1=A\|G\|C`. | Strongest current proof for robust discovery scoring; also proves default mode over-rewarded a tiny ambiguous haplotype. |
| TaGW2-A1 / TaGW2-6A | Jaiswal et al. 2015, PLoS ONE, `https://doi.org/10.1371/journal.pone.0129400` | Promoter `SNP-494` is the causal/expression-regulating marker; favorable allele `A`; four-SNP Hap5 pattern tracked as `SNP-988/SNP-739/SNP-593/SNP-494 = G/A/G/A`. | Local WatSeq chr6A SNP VCF plus remote WheatOmics micro-VCF. | Local data miss `SNP-494`; remote data cover `chr6A:237734341 G>A`, but expected `A` has 0 complete Watkins carriers, so Hap5 is not observed. | No exact validation. Some top haplotypes match the first three promoter SNPs, but not the decisive `SNP-494 A` state. | Gene-level TGW signal only; cannot be used as exact functional-haplotype proof in current samples. |
| TaGW2-D1 | WheatOmics GW2 homolog annotation and current WatSeq A/D run. | No precise published decisive D1 functional marker was configured for this validation; current target used two regional SNPs around `TraesCS6D02G176900`. | WatSeq chr6D SNP VCF; 820 Watkins TGW-overlap samples; 2 retained SNPs. | Region covered sparsely; no literature decisive allele/haplotype available in manifest. | Not applicable. | Low-confidence regional signal only; not a star-gene proof. |
| Rht-B1b | Rht-1 DELLA literature, including Peng et al. 1999 and Ellis et al. 2002 perfect-marker work; WheatOmics annotation used for coordinate. | Canonical semi-dwarf stop-gained DELLA allele, `chr4B:30861571 C>T`, `c.190C>T`, `p.Gln64*`; mechanism is GA-insensitive reduced height through truncated/stabilized DELLA repression of growth. | WheatOmics merged SNP VCF plus Watkins CFLN06 plant-height phenotype. | Marker is present, but expected `T` has 0 phenotype-overlap carriers. | No, blocked before scoring. | Data-blocked in current Watkins panel; not a method-negative result. |
| Rht-D1b | Rht-1 DELLA literature, including Peng et al. 1999 and Ellis et al. 2002 perfect-marker work; WheatOmics annotation used for coordinate. | Canonical semi-dwarf stop-gained DELLA allele, `chr4D:18781242 G>T`, `c.181G>T`, `p.Glu61*`; same GA-insensitive DELLA mechanism as Rht-B1b. | WheatOmics merged SNP VCF plus Watkins CFLN06 plant-height phenotype. | Covered and segregating, but only 5 functional-allele carriers overall; robust top haplotype carrier count is 3. | Yes. Default and robust top haplotype contains the functional `T` allele and has lower plant height. | Weak positive control: exact variant recovered, but sample size and R2 are too small for strong proof. |
| Rht-Zanke2014 | Zanke et al. 2014 PLOS ONE Table S2, `https://doi.org/10.1371/journal.pone.0113287`; candidate-gene genotypes use Ellis et al. Rht markers. | Diagnostic `Rht-B1b` and `Rht-D1b` marker states with multi-environment plant-height phenotypes. | Downloaded Table S2 workbook; 368 complete varieties; combined `Rht-B1/Rht-D1` marker panel. | Covered and segregating; `Rht-D1b` occurs in 216 samples and exact `Hap1=B1a\|D1b` has n=214. | Raw total-score top is high-plant-height wild type, but stable direction-aware top is exact `Rht-D1b` haplotype and audit status is `matched_directional_top_haplotype`. | Counts as proof only under the direction-aware validation view; also documents an algorithmic requirement for decreases_trait targets. |
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
816 Watkins TGW-overlap samples, 3 promoter SNPs, 6 haplotypes.

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
`291,759,689` (`A vs C=437; C/A=7`, n `371/445`, contrast `2.451`,
P `6.98e-09`, stable), `291,761,315`
(`C vs T=73; T/C=1`, shared allele warning), and `291,760,677`
(`G vs A=7`, shared allele/low-MAF warning). These correspond to the three
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
n=3, mean plant height `84.583`, score `1.1013`, reliability `0.1304`.
`Hap1` has n=797, mean `99.615`, score `0.0772`. Score regression
`R^2=0.005`, `P=0.0453`. Literature audit status:
`Rht-D1b_stop_gained = matched_top_haplotype`.

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
For `PlantHeight_BLUE`, raw top score is `Hap2=B1a|Rht-D1a`, n=126,
mean `96.384 cm`; the stable low-phenotype directional top is
`Hap1=B1a|D1b`, n=214, mean `82.540 cm`. Score regression is strong
(`R^2=0.5267`, `P=2.04e-61`), but raw score direction is opposite to the
expected `decreases_trait` direction.

For `PlantHeight_GAT_2012`, raw top score is again `Hap2`, while the stable
low-phenotype directional top is `Hap3=B1b|Rht-D1a`, n=25, mean `72.454 cm`;
score regression `R^2=0.6131`, `P=3.40e-75`.

Literature audit:
`Rht-D1b_diagnostic_marker` has raw `validation_status=present_but_not_top`,
but direction-aware audit reports `directional_top_haplotype=Hap1`,
`directional_top_haplotype_sample_count=214`,
`directional_validation_status=matched_directional_top_haplotype`.

Interpretation:
This is the third usable positive control only after adding direction-aware
validation. It shows that a score designed to find "large effect" haplotypes
must not be interpreted as "beneficial/top" without trait direction. For
unknown discovery, report high-value and low-value candidate ranks separately
when the trait direction is not known in advance.
