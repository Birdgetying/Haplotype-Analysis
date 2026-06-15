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

| Target | Trait | Score mode | Status | Key result | Conclusion |
|---|---|---|---|---|---|
| TaGW2-B1 / TaGW2-6B | TGW | default | tested | Top score was rare ambiguous `Hap5=C/A\|G\|C`, not exact Qin2014 `A\|G\|C` | Biological signal present, default ranking over-rewarded rare ambiguous haplotype |
| TaGW2-B1 / TaGW2-6B | TGW | robust_discovery | tested | Top score became exact Qin2014 `Hap1=A\|G\|C`, n=371 | Strongest current proof that robust scoring can recover a known functional haplotype without using literature labels in scoring |
| TaGW2-A1 | TGW | default | tested | Gene-level signal exists, but `SNP-494 A` is absent in current complete Watkins samples | Not an exact Hap5 validation |
| TaGW2-D1 | TGW | default | tested | Weak low-confidence signal with only 2 retained SNPs | Limited marker coverage; not strong proof |
| VRN-Kiss2014 | heading date | default | tested | Literature all-spring `1\|1\|1` was top-scored but n=2 | Useful audit check, weak statistical proof |
| VRN-Kiss2014 | heading date | robust_discovery | tested | Top score became `Hap2=0\|0\|1`, n=31; all-spring `1\|1\|1` was no longer top | Robust mode favors stable VRN-D1 signal over rare all-spring haplotype |
| Rht-D1b | plant height | default | tested | Top-scored haplotype contained functional T allele but only 5 carriers | Weak positive control |
| Rht-D1b | plant height | robust_discovery | tested | Top score remained functional T haplotype `Hap2`, n=3 | Weak positive control; carrier count too small |
| Rht-B1b | plant height | any | blocked | Functional T allele has 0 phenotype-overlap carriers | Cannot validate scoring in current samples |
| TaPIF4 | TGW / plant height | any | blocked | Public repository lacks final per-sample promoter haplotype table and matching phenotype table | Not usable until final data are obtained |

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
