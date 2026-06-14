# Wheat Star-Gene Positive-Control Audit

Date: 2026-06-14

## Decision rule

A target is accepted as a strong positive control only when all of these are true:

1. The literature functional variant or diagnostic haplotype is explicit.
2. The current data contain sample-level genotype/haplotype calls for that exact variant.
3. The same or overlapping samples have the relevant phenotype.
4. The scored target haplotype contains the literature functional state and the trait direction is consistent.

If a functional variant is absent, monomorphic, or represented only by a nearby/background haplotype, it is a data-blocked or weak result, not a method-negative result.

## Summary table

| Target | Trait | Functional variant to test | Current data status | Current verdict |
|---|---|---|---|---|
| Rht1 / Rht-D1b | Plant height | `chr4D:18781242 G>T`, `c.181G>T`, `p.Glu61*`, stop-gained in `TraesCS4D01G040400` | Present in WheatOmics merged SNP VCF; overlaps Watkins CFLN06 plant height for 802 samples, but only 5 carriers | Weak positive control: top-scored haplotype contains the functional T allele and height direction is lower, but carrier count is too small for strong proof |
| Rht1 / Rht-B1b | Plant height | `chr4B:30861571 C>T`, `c.190C>T`, `p.Gln64*`, stop-gained in `TraesCS4B01G043100` | Present in WheatOmics merged SNP VCF, but 0 T carriers among Watkins phenotype-overlap samples | Data-blocked in Watkins panel; cannot validate scoring with this phenotype set |
| VRN / spring-winter habit | Growth habit / flowering biology | VRN1 promoter and intron-1 variants, especially large deletions/CNV at `VRN-A1/VRN-B1/VRN-D1` | Watkins workbook has `GrwHabit_E_sw-CFLN06`, but no current diagnostic VRN marker table; ordinary SNP-only region data may miss intron deletion/CNV | Potentially good positive control, but not yet runnable without diagnostic marker/SV/CNV genotype calls |
| GW2 / TaGW2-A1 | Thousand grain weight | Jaiswal Hap5 promoter group: `SNP-988/SNP-739/SNP-593/SNP-494 = G/A/G/A`, with `SNP-494` as key expression-regulating SNP | Local SNP VCF misses `SNP-494`; remote SNP VCF contains it but Watkins complete samples are all `G` at `SNP-494`, expected `A` carriers = 0 | Gene-level signal exists for TaGW2-A1, but exact functional Hap5 is not validated in current samples |
| GW2 / TaGW2-B1 | Thousand grain weight | Qin et al. 2014 favorable promoter haplotype `Hap-6B-1`, diagnostic pattern `-1709/-721/-83 = A/G/C` | Functional promoter positions resolved to `chr6B:291759688/291760676/291761314`, but the current WheatOmics remote merged SNP VCF has no usable records for these three sites; local `D:\Desktop\data\GW2` still has chr6B `.csi` only and no chr6B `.vcf.gz` | Data-blocked for exact B1 validation until a chr6B genotype source covering these promoter SNPs is available |
| TaPIF4 | Heat-related grain size / plant architecture context | Author workflow classifies promoter coverage haplotypes from BAM depth, producing `PIF_hap.txt` | GitHub `QinZhen1995/CAU-TaSG` exposes scripts and `samlist`, not final per-sample `PIF_hap.txt` plus matching TGW/height table | Not suitable as current positive control unless raw BAM-derived coverage matrix or final haplotype table is obtained |

## Rht1 result

Implemented and ran `prepare_wheat2024_rht1_functional_snps.py`.

Data source:
- VCF: WheatOmics merged SNP VCF.
- Phenotype: Watkins JIC workbook, sheet `WGIN_Watkins_JIC_CFLN06`, column `PH_M_cm-CFLN06`.
- Output database: `star_gene_database/wheat_nature_2024/Rht-D1b`.
- Output HTML: `star_gene_results/wheat_nature_2024/Rht-D1b/Rht-D1b.html`.
- Audit CSV: `star_gene_results/wheat_nature_2024/Rht-D1b/literature_variant_audit.csv`.

Observed Rht-D1b groups:

| Haplotype | Functional state | n | Mean height |
|---|---:|---:|---:|
| Hap1 | `G` | 797 | 99.62 cm |
| Hap2 | `T` | 3 | 84.58 cm |
| Hap3 | `G/T` | 2 | 90.50 cm |

Scoring/audit:
- `Rht-D1b_stop_gained` is covered and segregating.
- Carrier count is 5.
- Top-scored haplotype is Hap2, score `4.4556`.
- Hap2 contains the literature T allele.
- Validation status is `matched_top_haplotype`.
- Score regression is direction-consistent for lower height, but low-confidence: `R^2=0.006`, `P=0.0288`, `PVE=0.61%`.

Interpretation:
This proves the pipeline can assign high score to the exact Rht-D1b stop-gained allele when it is present, but the Watkins panel has too few carriers to be a strong positive-control proof.

Rht-B1b status:
- The exact functional SNP is found in the remote VCF.
- In the Watkins phenotype-overlap samples the expected T allele has 0 carriers.
- This is a population/data blocker, not evidence against the scoring method.

## GW2 result

Current TaGW2-A1 outputs:
- Local SNP-only: `star_gene_results/wheat_nature_2024/TaGW2-A1/literature_variant_audit.csv`.
- Remote SNP rerun: `star_gene_results/wheat_nature_2024/TaGW2-A1-remoteSNP/literature_variant_audit.csv`.

Strict functional-haplotype audit:
- `SNP-988`, `SNP-739`, and `SNP-593` match the top-scored TaGW2-A1 haplotype.
- Local SNP-only data miss `SNP-494`.
- Remote SNP data cover `SNP-494`, but all 552 complete Watkins samples are `G`; expected `A` carrier count is 0.
- Full `Jaiswal_Hap5 = G/A/G/A` is therefore not observed.

Interpretation:
TaGW2-A1 supports a gene-level signal, but it does not yet prove recovery of the published functional Hap5.

TaGW2-B1 follow-up:
- Qin et al. 2014 reports the favorable natural promoter haplotype `Hap-6B-1`.
- From Fig. 2B, the three diagnostic promoter sites are `-1709=A`, `-721=G`, and `-83=C`.
- The local WheatOmics sequence header for `TraesCS6B02G215300` gives gene start `chr6B:291,761,229` and `CDS=169`, so ATG maps to `chr6B:291,761,397`.
- Therefore the three diagnostic sites map to `chr6B:291,759,688`, `chr6B:291,760,676`, and `chr6B:291,761,314`.
- Added `prepare_wheat2024_tagw2_b1_remote_snp.py` and manifest target `TaGW2-B1-remoteSNP` so this route is reproducible.
- Running the script against the current WheatOmics remote merged SNP VCF blocked with: literature promoter SNPs are missing or not segregating at all three marker IDs.
- The Earlham chr6B SNP URL currently returns a small HTML response (`Content-Type: text/html`, `Content-Length: 2261`) rather than the VCF body, and the local folder has only the `.csi` index. Thus B1 cannot yet validate scoring.

Interpretation:
TaGW2-B1 is still the best GW2 biological target, but current genotype data do not cover its published diagnostic promoter haplotype. This is a data blocker, not a negative method result.

## VRN result

Local phenotype status:
- The Watkins workbook has growth habit in `GrwHabit_E_sw-CFLN06`.
- Most values are coded strings such as `ssss`, `wwww`, `wwss`, etc., which can be converted to spring/winter or mixed classes after defining a rule.

Data blocker:
- The important VRN functional variants are often structural or copy-number variants, not just simple SNPs.
- Current local data do not yet include a diagnostic VRN marker table or verified SV/CNV calls for `VRN-A1/VRN-B1/VRN-D1`.

Interpretation:
VRN is conceptually a strong positive control because the phenotype is categorical and large-effect. It should be the next target only after obtaining diagnostic marker/SV/CNV genotypes.

Potential external route:
Kiss et al. 2014 used diagnostic molecular markers for `VRN-A1`, `VRN-B1`, `VRN-D1`, `PPD-B1`, and `PPD-D1` in a 683-genotype worldwide wheat collection with heading phenotypes. This could become a strong independent VRN/heading positive control if the supplementary per-genotype marker table can be downloaded and converted. The Springer supplementary DOC download timed out from this machine in the current run, so it is not yet an available dataset.

## TaPIF4 result

Public repository checked:
- `https://github.com/QinZhen1995/CAU-TaSG`
- `03.PIF4_Promoter_Hap` contains `GeneToCoverageHeatMap.sh`, `martix2hap.sh`, `plotcoverage.R`, and `samlist`.

Important detail:
- `martix2hap.sh` creates `PIF_hap.txt` from a `coverage.martix`.
- The repository does not provide the final `coverage.martix`, `PIF_hap.txt`, or a same-sample TGW/height phenotype table.

Interpretation:
TaPIF4 is not currently usable as a clean positive control. It can become usable only if the final per-sample promoter haplotype table and matching phenotype data are obtained.

## Ranking for proving scoring effectiveness

1. Rht-D1b: usable but weak; exact functional variant recovered, carrier count too low.
2. GW2-A1: gene-level signal, but exact published Hap5 not recovered because `SNP-494 A` is absent.
3. TaGW2-B1: biologically strong and exact promoter haplotype now mapped, but current VCF sources do not cover the three diagnostic SNPs.
4. VRN: promising, but blocked by missing diagnostic structural/CNV marker calls.
5. TaPIF4: currently unsuitable, because public data expose scripts rather than the final per-sample haplotype/phenotype table.

Best next data action:
Obtain a chr6B genotype source covering the TaGW2-B1 promoter diagnostic SNPs, or obtain a per-genotype VRN diagnostic marker table with matching heading/growth-habit phenotype.
