# Star Gene Validation Data Adaptation Plan

## Goal
Use public data from three crop genomics papers to validate whether the haplotype scoring workflow can recover known star genes/loci.

## Constraints
- Do not delete files without explicit user confirmation.
- Do not guess unresolved coordinates, assemblies, data formats, or trait directions.
- Do not download very large datasets by default; provide reproducible download instructions first.
- Code changes require tests, git commit, and push when network permits.

## Phases

### Phase 1: Project context
Status: complete

Find existing star-gene validation entry points, manifests, docs, and current repository state.

### Phase 2: Public data source research
Status: complete

For each paper, identify official data repositories, file formats, file sizes if available, access terms, and whether target-region or small marker-level data can support the first validation pass.

### Phase 3: Design options
Status: complete

Present 2-3 data ingestion and code adaptation approaches with tradeoffs. Ask the user to approve one before implementation.

### Phase 4: Implementation
Status: complete

After approval, update manifest/download helpers/adapters/tests as needed.

### Phase 5: Verification and git
Status: complete

Run required tests, commit modified files, and push to GitHub if network allows.

### Phase 6: Maize qHKW1/ZmBAM1d positive-control run
Status: complete

Download the small MaizeGo SV/pSV and BLUP phenotype files, inspect whether they contain sample-level qHKW1/ZmBAM1d 8.9 kb indel information and hundred-kernel-weight phenotypes, then either build a precomputed database and run analysis or record the exact data blocker.

Outcome: the lightweight MaizeGo path now runs end to end from a marker matrix, but the count-compatible public marker did not provide a positive validation signal for qHKW1. The exact paper 8.9-kb indel marker remains unresolved in the small public resources.

### Phase 7: Next positive-control target selection
Status: complete

Choose the next validation route: obtain the exact qHKW1 8.9-kb indel marker from a paper-specific source/contact, or move to a rice/wheat target with confirmed sample-level genotype and phenotype files.

Current decision: keep qHKW1/ZmBAM1d as the first priority, but require the paper genotype marker before interpreting the result. Added a strict preparation script for the MaizeGo full `SV.386014.zip` package. It scans only HapMap/SV matrix files for a chr1 paper-window 8.5-9.5 kb candidate and builds the database only if exactly one marker is found. With currently available small MaizeGo matrices, the script finds zero exact candidates, so qHKW1 remains data-blocked rather than method-negative.

### Phase 8: Paper-genotype rerun after full MaizeGo SV package
Status: complete

`SV.386014.zip` is present at `external_data/maize_natgenet_2019/maizego/SV.386014.zip` and has been extracted. The package contains `svs.final.ms.txt` and `svs.final.bs.txt`, which are no-header SV catalogue files, not sample-level accession-by-marker genotype matrices.

Running:

```bash
python prepare_maize2019_qhkw1_paper_genotype.py --force-extract
```

found 6 candidate matrix files from the small MaizeGo resources, 2 SV catalogue files from the full package, 0 sample-level qHKW1 paper-window 8.5-9.5 kb marker candidates, and 0 catalogue candidates in chr1 `30.44-30.54 Mb`. An expanded chr1 `27-31 Mb` catalogue scan found 8 MS catalogue records of 8.5-9.5 kb, but these have no sample genotype columns and cannot be used to build haplotypes or score associations.

Outcome: qHKW1/ZmBAM1d remains data-blocked. Do not run or interpret `python run_star_gene_validation.py --run-analysis --paper maize2019 --target qHKW1` until a per-accession genotype table for the paper 8.9-kb indel is available.

### Phase 9: Next actionable validation route
Status: complete

Pick a positive-control target with confirmed sample-level genotype plus phenotype. Options:
- obtain the author/MaizeGo region-analysis qHKW1 8.9-kb indel per-accession genotypes;
- pivot to rice 18K or wheat WatSeq target-region data where public files expose sample-level genotypes for known genes;
- create a narrow target-region extraction plan once the relevant VCF/matrix filenames are confirmed.

Current qHKW1 acquisition status: Nature supplementary files and public MaizeGo downloads do not expose the exact Fig. 4h per-accession genotype table. The public MaizeGo matrices include nearby 8.5-9.5 kb sample-level SV rows, including `chr1_28599370_28608270_deletion_8900`, but this row does not match the paper group sizes and is not sufficient as the positive-control table. Best next action is to request the exact qHKW1/ZmBAM1d Fig. 4h indel genotype table from the MaizeGo/Yan lab contact for specific gene or region analysis.

User decision: pause the maize qHKW1 route and proceed with the other two papers first.

### Phase 10: Rice/Wheat positive-control execution
Status: in_progress

Wheat Q7B-PH is the first completed positive-control route from the remaining two papers. The minimal first pass uses WWWG2B Figure 3g source data (`Accession`, `allele`, `PH_M_cm`) and treats the paper-defined `allele` as the Q7B-PH haplotype label. This is a paper-source haplotype positive control, not a chr7B physical-interval VCF discovery run.

Commands:

```bash
python prepare_wheat2024_q7b_ph_figure3g.py
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Q7B-HT
```

Current Q7B-PH result: 139 plot-level observations, 3 haplotypes, 1 paper-defined marker, corrected association p-value `0.02403`, PVE `29.63%`, HaplotypeScorer regression `R^2=0.2627`, regression p-value `1.119e-10`, high confidence, and direction consistency for the expected height-increasing haplotype.

TaGW2 A/D WatSeq VCF run is also complete for the locally downloaded chr6A and chr6D SNP VCFs. `TaGW2-A1` (`TraesCS6A02G189300`) has 709 samples, 55 SNPs, 19 haplotypes, corrected association P value `7.07e-06`, PVE `16.30%`, HaplotypeScorer `R^2=0.0431`, score-regression P value `2.416e-08`, and moderate confidence. `TaGW2-D1` (`TraesCS6D02G176900`) has 820 samples, 2 SNPs, 3 haplotypes, corrected association P value `1.0`, PVE `3.95%`, HaplotypeScorer `R^2=0.0257`, score-regression P value `3.893e-06`, and low confidence. This is A/D-only; `TaGW2-B1` remains excluded until the chr6B VCF download succeeds.

TaGW2-A1 Hap8 functional follow-up is complete for the available SNP data. `Hap8` is not a one-to-one proxy for the published local two-SNP promoter state: the `-739/-593` state maps to `chr6A:237,734,096/237,734,242` and occurs in 66 samples across `Hap7`, `Hap8`, `Hap15`, and `Hap17`, while full `Hap8` has 25 samples. The Hap8 background/full Hap8 effect on `TGW_mean` is stronger (`-7.55 g`, Welch `P=3.36e-07`) than the two-SNP state alone (`-3.18 g`, `P=7.99e-04`). Country origin is a confounder because the two-SNP state is enriched in China, but a within-country centered test still leaves a Hap8 effect (`-4.74 g`, `P=1.33e-04`). INDEL/SV has not been merged: no local chr6A INDEL/SV VCF exists, and the Earlham chr6A INDEL VCF is about 7.94 GiB with extraction requiring `bcftools`/`tabix` or `pysam`.

Rice Science 2024 remains data-blocked locally. Requests from this machine to Figshare article/API/ndownloader for `10.6084/m9.figshare.19166475` returned HTTP 403, and no `external_data/rice_science_2024` files are currently present. Do not claim rice OsMADS22/OsFTL1 validation until the per-accession genotype and phenotype files are available and inspected.

## Open Questions
- Which validation depth should be prioritized: full-paper datasets, target-region VCFs, or marker/phenotype positive-control tables?
- Whether the user wants local large downloads on this machine or only scripts/instructions.
- For maize, the exact paper 8.9-kb qHKW1 sample-level marker was not identified in the small MaizeGo resources or the downloaded full `SV.386014.zip`; the count-compatible marker and catalogue-only records are not sufficient as proof.
- For rice, the next blocker is data acquisition rather than code: need a reachable Figshare mirror, browser/manual download, or user-provided `NAM_variations` files before target adapters can be written.
- For wheat, the next stronger test is chr7B VCF/INDEL extraction around Q7B-PH after confirming the physical interval; the current Figure 3g result validates scoring against paper-defined haplotype labels only.

## Errors Encountered
- Figshare API and Science article pages returned HTTP 403 from local PowerShell requests. Used accessible public pages/search snippets and official article/data URLs instead.
- On 2026-06-05, Figshare API, private-link API, and ndownloader URLs for article `19166475` still returned HTTP 403 from local PowerShell requests.
- Earlham WatSeq root and guessed subdirectory URLs returned HTTP 404. Paper confirms the OpenData path, but exact downloadable VCF object names should be obtained from WWWG2B/Earlham metadata or portal interaction before hard-coding.

### Phase 11: Literature functional-variant audit table
Status: complete

Implemented a manifest-driven audit that writes `literature_variant_audit.csv` and `literature_variant_audit.json` for analyzed star-gene targets with configured literature variants or literature haplotypes.

Current TaGW2-A1 audit configuration tracks the Jaiswal promoter SNP group:

- `SNP-988` = `chr6A:237,733,847`, expected allele `G`
- `SNP-739` = `chr6A:237,734,096`, expected allele `A`
- `SNP-593` = `chr6A:237,734,242`, expected allele `G`
- `SNP-494` = `chr6A:237,734,341`, expected allele `A`
- `Jaiswal_Hap5` = `G/A/G/A` across the four markers

The audit logic reports whether each literature marker is covered by the current database, whether the expected allele or full literature haplotype segregates in the analyzed samples, whether the top-scored haplotype contains it, and a strict `validation_status`. Monomorphic or absent expected alleles are not accepted as validation even if the top haplotype carries the reference/background state.

### Phase 12: Wheat teacher star-gene triage
Status: complete

Teacher-provided wheat star-gene candidates were audited one by one: Rht1 for plant height, VRN for spring/winter growth habit, GW2 for grain weight, and TaPIF4 for grain weight/height context.

Implemented `prepare_wheat2024_rht1_functional_snps.py` to build Rht-B1b/Rht-D1b functional-SNP positive-control databases from the WheatOmics merged SNP VCF and Watkins CFLN06 plant-height phenotypes. The exact Rht-B1b stop-gained SNP is `chr4B:30861571 C>T` (`c.190C>T`, `p.Gln64*`), and the exact Rht-D1b stop-gained SNP is `chr4D:18781242 G>T` (`c.181G>T`, `p.Glu61*`).

Current Rht1 result: `Rht-B1b` is data-blocked in the Watkins phenotype-overlap samples because the expected T allele has 0 carriers. `Rht-D1b` builds and analyzes, but only 5 carriers are present. The top-scored haplotype (`Hap2`) contains the functional T allele and the direction is lower plant height, but the result is low-confidence (`R^2=0.006`, `P=0.0288`, `PVE=0.61%`), so it is a weak positive control rather than strong proof.

Current VRN result: Watkins has `GrwHabit_E_sw-CFLN06`, but no verified diagnostic VRN structural/CNV marker table is available. VRN is promising only after obtaining `VRN-A1/VRN-B1/VRN-D1` diagnostic marker/SV/CNV genotypes.

Current GW2 result remains unchanged but stricter: TaGW2-A1 shows gene-level signal, yet the exact Jaiswal Hap5 is not validated because `SNP-494 A` is absent in current complete Watkins samples. TaGW2-B1 remains high priority after obtaining the chr6B VCF main file and exact functional marker.

Current TaPIF4 result: the public `QinZhen1995/CAU-TaSG` GitHub repository exposes promoter-haplotype scripts and `samlist`, but not the final per-sample `PIF_hap.txt`, `coverage.martix`, or matching phenotype table. It is not currently suitable as a positive control.

Detailed report: `wheat_star_gene_positive_control_audit.md`.

### Phase 13: TaGW2-B1 and VRN stronger positive-control follow-up
Status: complete

TaGW2-B1 is now downloaded, rebuilt, and analyzed. Qin et al. 2014 reports the favorable natural promoter haplotype `Hap-6B-1` with diagnostic states `-1709=A`, `-721=G`, and `-83=C`. Using the local WheatOmics `TraesCS6B02G215300` sequence header (`chr6B:291,761,229-291,778,752`, `CDS=169`), ATG maps to `chr6B:291,761,397`. The first coordinate pass used `chr6B:291,759,688/291,760,676/291,761,314`, but the complete WWWG2B 1047 chr6B SNP VCF contains the corresponding variant records at `chr6B:291,759,689 C>A`, `chr6B:291,760,677 G>A`, and `chr6B:291,761,315 T>C`; the script and manifest now use these VCF coordinates.

Implemented `download_wwwg2b_file.py` for resumable WWWG2B/OneDrive downloads. It refreshes the temporary URL before each `curl -C -` attempt and moves HTML/non-gzip responses to `.invalid`. The complete chr6B `SNP_matrix_1047` VCF is present at `D:\Desktop\data\GW2\chr6B.HARD.SNP.Missing-unphasing.ID.ann.finalSID.1047.allele2_retain.hard_retain.InbreedingCoeff_retain.vcf.gz` (`2,117,273,732` bytes), with `.csi` index.

Implemented and updated `prepare_wheat2024_tagw2_b1_remote_snp.py` and manifest target `TaGW2-B1-remoteSNP`. The script is strict about all three literature promoter SNPs being present and segregating, and now uses indexed `bcftools` extraction for local `.vcf.gz` files to avoid scanning a 2.1 GB gzip linearly. Running it against the downloaded chr6B VCF built 816 samples, 3 variants, and 6 haplotypes.

Real-data TaGW2-B1 result: `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP` generated HTML and audit outputs. The exact Qin2014 favorable pattern `A|G|C` is `Hap1` (n=371) and has the highest observed TGW mean (`40.970 g`), while `C|G|C` is lower (`39.049 g`) and `C|G|T` is much lower (`34.997 g`). This supports the biological direction of the literature haplotype. However, the default HaplotypeScorer top rank is `Hap5 = C/A|G|C` (n=7, score `4.2581`, TGW mean `39.610 g`), not the exact literature haplotype. The strict literature audit marks the combined Qin2014 haplotype as `contained_in_top_haplotype_not_exact`, with exact matching haplotypes `Hap1:371`. Therefore TaGW2-B1 is a strong positive-control dataset for functional direction but also evidence that the top-rank scoring rule needs refinement before being used alone for unknown-gene discovery.

### Phase 14: Robust discovery scoring mode

Status: complete

Implemented `score_mode=robust_discovery` as a new scoring mode rather than changing the historical default. This mode is designed for discovery ranking and does not use literature-positive labels. It reduces over-weighted duplicate rare/special signals, applies sample reliability `n/(n+20)` to unstable components, penalizes ambiguous allele tokens such as `C/A`, and writes reliability diagnostics into `haplotype_scores.json`.

`run_star_gene_validation.py --score-mode robust_discovery` is now wired through the analyzer and writes to a suffix directory, e.g. `star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery`, so default and robust outputs can be compared directly.

Real-data TaGW2-B1 robust result: `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --score-mode robust_discovery` ranks exact Qin2014 `Hap1=A|G|C` first (`n=371`, score `1.8132`, reliability `0.9488`, ambiguity `1.0`) and demotes rare ambiguous `Hap5=C/A|G|C` to third (`n=7`, score `0.9733`, reliability `0.2593`, ambiguity `0.7333`). The robust literature audit marks `Qin2014_Hap-6B-1` as `matched_top_haplotype`. This is the clearest current evidence that the scoring method can recover a published star-gene functional haplotype after adding discovery-oriented stability controls.

VRN-Kiss2014 is now implemented and runnable. The first local Kiss DOC was incomplete (`3,133,440` bytes); the complete Springer static DOC (`11032_2014_34_MOESM1_ESM.doc`, `8,219,648` bytes) was downloaded with resumable `curl -L -C -`. Its embedded Excel Workbook was extracted to `external_data/literature/vrn_kiss2014/Kiss2014_embedded_workbook.xls`.

Implemented `prepare_wheat2024_vrn_kiss2014.py` and manifest target `VRN-Kiss2014`. The default database uses the three VRN diagnostic marker states (`VRN-A1`, `VRN-B1`, `VRN-D1`) from Kiss et al. 2014 and heading-date phenotypes DEV49/DEV59. Running `python prepare_wheat2024_vrn_kiss2014.py` built 676 complete samples, 3 markers, and 7 haplotypes. Running `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-Kiss2014` generated HTML and literature audit outputs.

Current VRN result: DEV49_mean association P `1.01e-05`, PVE `4.83%`, score R² `0.0262`, score P `2.38e-05`; DEV59_mean association P `0.00168`, PVE `3.13%`, score R² `0.0147`, score P `0.00159`. The top-scored haplotype is `Hap7 = 1|1|1`, and the literature audit marks all three VRN spring diagnostic states plus the combined `1|1|1` diagnostic haplotype as `matched_top_haplotype`. However, `Hap7` has only 2 samples, so this is a positive algorithmic check but not strong population-statistical proof. The most stable effect is `Hap2 = 0|0|1` (VRN-D1 spring state only), which heads earlier than `Hap1 = 0|0|0` by `-3.36` days for DEV49 and `-2.98` days for DEV59.

Next action: rerun the robust discovery mode on VRN-Kiss2014 and Rht-D1b, then compare default vs robust across the positive-control panel. If robust keeps GW2-B1 correct but over-penalizes rare true positives such as VRN all-spring haplotype, the next refinement should report both a stable common-haplotype rank and a rare-candidate rank instead of forcing one ranking to serve both purposes.

### Phase 15: Literature decisive-variant checklist

Status: complete

Updated `star_gene_validation_record.md` with a dedicated checklist that links each teacher star-gene target to the exact literature validation object and the current score result. The checklist now separates strong proof, weak proof, marker-level proof, and data-blocked cases:

- TaGW2-B1/TaGW2-6B: strongest positive control. Qin2014 `Hap-6B-1` (`A/G/C` at `-1709/-721/-83`) is covered and segregating; robust discovery ranks exact `Hap1=A|G|C` first.
- Rht-D1b: exact stop-gained `chr4D:18781242 G>T`, `c.181G>T`, `p.Glu61*` is recovered by the top haplotype, but carrier count is too small for strong proof.
- VRN-Kiss2014: diagnostic marker-level validation only; the exact causal VRN structural variants are not encoded, and robust mode favors the more stable `VRN-D1=1` signal.
- TaGW2-A1: Jaiswal Hap5/SNP-494 is not validated because the expected `SNP-494 A` allele is absent in current complete Watkins samples.
- Rht-B1b: exact stop-gained marker is present in the source VCF but has 0 phenotype-overlap carriers.
- TaPIF4/TaSG-D1-TaPIF4: literature mechanism and promoter InDel classes were checked, but public GitHub currently exposes scripts/samlist rather than final per-sample promoter haplotypes plus matching phenotypes.

### Phase 16: Three-positive-control proof set and direction-aware validation

Status: complete

Implemented the additional proof route needed for the teacher's "at least 3" requirement and documented the resulting evidence:

- Proof #1 remains `TaGW2-B1-remoteSNP__robust_discovery`: exact Qin2014 `Hap-6B-1=A|G|C` is robust top-ranked and audit-matched.
- Proof #2 is `VRN-D1-Kiss2014__robust_discovery`: single-marker `VRN-D1=1` is robust top-ranked with significant DEV49/DEV59 heading-date score regressions.
- Proof #3 is `Rht-Zanke2014__robust_discovery` after direction-aware scoring: Zanke et al. 2014 PLOS ONE Table S2 provides Rht-B1/Rht-D1 marker genotypes and plant heights for 368 varieties. The active `directional_total` score now ranks `Hap1=B1a|D1b` first (n=214), and literature audit marks `Rht-D1b` as both `matched_top_haplotype` and `matched_directional_top_haplotype`. The retained raw `total` explains the old high-plant-height Hap2 artifact.

Code changes in this phase:

- Added `prepare_wheat2024_rht_zanke2014.py`.
- Added manifest target `Rht-Zanke2014`.
- Added `directional_top_haplotype*` fields to validation summary.
- Added `directional_validation_status` and related fields to `literature_variant_audit.csv/json`.
- Added tests for Rht Zanke2014 preparation, direction-aware top selection, and direction-aware literature audit.

TaPIF4 was checked but not used as a proof. `Supplementary_Data_5` has TaPIF4 haplotypes for 331 accessions but no matching phenotype table; the official Source Data ZIP is reachable but downloads were incomplete/corrupt during this session.

### Phase 17: Rht strict single-variant clarification

Status: complete

After the teacher objected that `Rht-Zanke2014` combines multiple Rht genes, the strict Rht proof route was narrowed back to one gene and one functional base variant.

Commands:

```bash
python prepare_wheat2024_rht1_functional_snps.py --min-haplotype-count 1 --target Rht-D1b
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Rht-D1b --score-mode robust_discovery
```

Result:

`Rht-D1b` uses the exact literature stop-gained SNP `chr4D:18781242 G>T` (`c.181G>T`, `p.Glu61*`) from the WheatOmics merged SNP VCF and Watkins CFLN06 plant-height phenotype overlap. The 2026-06-26 rerun rebuilt one marker and 802 samples, then regenerated `star_gene_results/wheat_nature_2024/Rht-D1b__robust_discovery/Rht-D1b.html`.

The strict top raw robust haplotype is `Hap2=T` with n=3, mean `84.583 cm`, score `1.0938`, and audit status `Rht-D1b_stop_gained=matched_top_haplotype`. Total functional-allele carriers are 5 (`G:797;G/T:2;T:3`). Because support-shrunk directional top selects the large wild-type `Hap1=G`, this remains weak exact-SNP validation rather than strong proof.

Decision:

`Rht-Zanke2014` remains useful as a marker-panel/directionality reference, but it should not be used as strict Rht proof because it combines `Rht-B1` and `Rht-D1` marker states.

Implementation note:

The first 2026-06-26 `--target Rht-D1b` run also selected `Rht-Zanke2014` because alias matching treated the marker-panel target's `Rht-D1b` alias as equal to the exact target. `StarGeneValidator.iter_targets()` now gives exact `target_id` filters precedence over alias matches, with a regression test in `test_star_gene_data.py`.

### Phase 18: Wheat INDEL/SV supplement acquisition

Status: complete

After the user asked to download the other needed data, the reachable public
WheatOmics INDEL resources were downloaded or extracted:

- `download_wheat2024_indel_microvcfs.py` extracts target-region micro-VCFs
  from `WEC_INDEL_IWGSCv1.0.eff.vcf.gz`,
  `GBS_filtered_Indels_IWGSCv1.0.eff.vcf.gz`, and
  `WildEmmer10WGS_INDEL_eff.vcf.gz`.
- `download_wheat2024_indel_marker_annotations.py` extracts nearby marker
  annotations from the WheatOmics `Indel_marker_from_zhai` JBrowse track.

Outputs:

- `external_data/wheat_nature_2024/indel_microvcfs/microvcf_status.tsv`
- `external_data/wheat_nature_2024/zhai_indel_markers/zhai_indel_marker_status.tsv`

Result:

The public INDEL VCF tracks do not strengthen the Watkins validations because
they have 0 Watkins phenotype-sample overlap. VRN target regions have 0 public
INDEL VCF records. The Zhai marker annotation track has nearby marker rows for
`TaGW2-B1` and `VRN-B1` with a 1 Mb flank, but it contains marker/primer
annotations rather than per-accession genotype calls.

Blocked data:

The still-needed sample-level WatSeq INDEL/SV/CNV genotype files could not be
downloaded automatically in this pass. WWWG2B API endpoints currently return
HTTP 526, and Earlham WatSeq URLs return HTML rejection pages rather than VCF
data.

### Phase 19: Rht-D1 single-marker Zanke2014 clarification

Status: complete

After the teacher objected that `Rht-Zanke2014` combines `Rht-B1` and
`Rht-D1`, the Zanke2014 Table S2 adapter was extended with
`--single-marker-targets`. This builds `Rht-B1-Zanke2014` and
`Rht-D1-Zanke2014` from the same workbook while keeping each target to one
marker column.

Commands:

```bash
python prepare_wheat2024_rht_zanke2014.py --single-marker-targets
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Rht-D1-Zanke2014 --score-mode robust_discovery
```

Output:

`star_gene_results/wheat_nature_2024/Rht-D1-Zanke2014__robust_discovery/Rht-D1-Zanke2014.html`

Result:

`Rht-D1-Zanke2014` has 368 varieties and one marker. The haplotypes are
`Hap1=D1b` n=214, `Hap2=Rht-D1a` n=152, and mixed `Hap3=D1b/Rht-D1a` n=2.
The raw robust total top is the two-sample mixed `Hap3`, but the
direction-aware top is exact `Hap1=D1b` with 214 samples. The literature audit
reports `Rht-D1b_diagnostic_marker` as
`matched_directional_top_haplotype`.

Interpretation:

This answers the "do not combine multiple genes" concern at marker level. It
does not replace the exact SNP `Rht-D1b` route, because Zanke2014 encodes a
diagnostic marker state rather than the nucleotide `chr4D:18781242 G>T`
record. The two Rht routes should be described together: exact SNP evidence is
weak but base-level; Zanke2014 single-marker evidence is stronger but
marker-level.

### Phase 20: Phenotype-free weighted-site robust scoring

Status: complete

Goal:

Move `robust_discovery` closer to the PostGWAS-style model requested by the
user: score weighted sites first, then let weighted sites contribute to
haplotype discovery. The discovery score must not use the current validation
phenotype.

Implemented:

- Added `site_weighted` as a robust-discovery component.
- Added `site_weights` and `site_weighting_policy` to
  `haplotype_scores.json` for auditability.
- Site weights use annotation severity, explicit external evidence,
  MAF stability, missingness, LD pruning, and gene-structure context.
- Current phenotype, haplotype means, current-trait effect sizes, and local
  current-trait p-values are forbidden inputs.
- Rare high-impact functional sites are retained even when MAF is below the
  ordinary background threshold.
- Normalized snpEff `stop_gained` to internal `stop_gain`, fixing the strict
  `Rht-D1b` site-weight omission.

Verification command:

```bash
python -m unittest test_star_gene_data.StarGeneDataTests.test_site_weight_keeps_rare_high_impact_functional_sites test_star_gene_data.StarGeneDataTests.test_site_weight_keeps_snp_eff_stop_gained_annotation test_star_gene_data.StarGeneDataTests.test_robust_discovery_site_weighted_score_is_phenotype_free test_star_gene_data.StarGeneDataTests.test_robust_discovery_score_is_unchanged_when_phenotype_is_inverted test_star_gene_data.StarGeneDataTests.test_robust_discovery_uses_only_explicit_external_site_evidence test_star_gene_data.StarGeneDataTests.test_robust_discovery_adds_boundary_gene_body_signal_to_functional_groups test_star_gene_data.StarGeneDataTests.test_robust_discovery_penalizes_tiny_ambiguous_haplotypes test_star_gene_data.StarGeneDataTests.test_expected_decrease_direction_reorders_raw_high_phenotype_score -v
python -m py_compile haplotype_phenotype_analysis.py star_gene_validation.py star_gene_literature_audit.py test_star_gene_data.py
```

Validation rerun:

```bash
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --target VRN-D1-Kiss2014 --target Rht-D1b --target Rht-Zanke2014 --score-mode robust_discovery
```

Result:

- `VRN-D1-Kiss2014`: positive marker-level proof; top and directional top are
  `Hap2`, with one phenotype-free diagnostic-marker site weight.
- `Rht-D1b`: exact base-level proof restored at the scoring layer; the
  `chr4D:18781242 G>T` stop-gained SNP is now included as a phenotype-free
  `stop_gain` site weight and raw top is `Hap2=T`, but support is still weak.
- `TaGW2-B1-remoteSNP`: still a full-region stress test. Qin2014 evidence is
  best read through directional/function-group audit, not raw top score.
- `Rht-Zanke2014`: useful marker-panel directionality evidence, but not strict
  single-gene proof because B1 and D1 markers are combined in this target.

### Phase 21: Base-level genotype scoring and data search

Status: in_progress

Goal:

Use true per-accession REF/ALT genotype calls for validation whenever possible.
Targets built from marker labels such as `D1b`, `VRN-D1=1`, `P/W/PAR`, or
paper haplotype names must be treated as marker-level evidence only. They can
help explain biology, but they cannot be counted as strict base-level proof.

Current base-level rerun:

```bash
python prepare_wheat2024_vrn_remote_snps.py
python prepare_wheat2024_rht1_functional_snps.py --min-haplotype-count 1 --target Rht-D1b
python prepare_wheat2024_tagw2_b1_remote_snp.py --vcf D:\Desktop\data\GW2\chr6B.HARD.SNP.Missing-unphasing.ID.ann.finalSID.1047.allele2_retain.hard_retain.InbreedingCoeff_retain.vcf.gz --min-haplotype-count 1 --max-missing-rate 0.2
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-A1-remoteSNP --target VRN-B1-remoteSNP --target VRN-D1-remoteSNP --score-mode robust_discovery
python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Rht-D1b --target TaGW2-B1-remoteSNP --score-mode robust_discovery
```

Output:

- `star_gene_results/wheat_nature_2024/VRN-A1-remoteSNP__robust_discovery`
- `star_gene_results/wheat_nature_2024/VRN-B1-remoteSNP__robust_discovery`
- `star_gene_results/wheat_nature_2024/VRN-D1-remoteSNP__robust_discovery`
- `star_gene_results/wheat_nature_2024/Rht-D1b__robust_discovery`
- `star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP__robust_discovery`

Interim result:

- `Rht-D1b` is strict base-level evidence: exact SNP
  `chr4D:18781242 G>T` is recovered as top robust haplotype `Hap2=T`, but only
  5 carriers overlap phenotype data, so it remains weak proof.
- `TaGW2-B1-remoteSNP` is base-level regional SNP evidence: 71 SNPs and 756
  samples. Literature Qin2014 promoter alleles are covered, but full-region
  raw top is not the exact literature promoter haplotype; use this as
  direction-aware/full-region stress evidence, not raw-top-only proof.
- `VRN-A1/B1/D1-remoteSNP` are base-level SNP-only regional runs. They detect
  growth-habit signal, especially `VRN-B1`, but they do not include the known
  VRN promoter/intron-1 deletion or CNV/SV alleles and therefore cannot yet
  replace causal VRN validation.

Data search status:

WWWG2B APIs currently return HTTP 526, and Earlham OpenData WatSeq URLs return
`Request Rejected` HTML. The next data priority remains sample-level
WatSeq/WWWG2B INDEL/SV/CNV VCFs or another accession-matched causal-variant
table for VRN and TaPIF4.
