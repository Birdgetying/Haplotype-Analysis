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

VRN-Kiss2014 is now implemented and runnable. The first local Kiss DOC was incomplete (`3,133,440` bytes); the complete Springer static DOC (`11032_2014_34_MOESM1_ESM.doc`, `8,219,648` bytes) was downloaded with resumable `curl -L -C -`. Its embedded Excel Workbook was extracted to `external_data/literature/vrn_kiss2014/Kiss2014_embedded_workbook.xls`.

Implemented `prepare_wheat2024_vrn_kiss2014.py` and manifest target `VRN-Kiss2014`. The default database uses the three VRN diagnostic marker states (`VRN-A1`, `VRN-B1`, `VRN-D1`) from Kiss et al. 2014 and heading-date phenotypes DEV49/DEV59. Running `python prepare_wheat2024_vrn_kiss2014.py` built 676 complete samples, 3 markers, and 7 haplotypes. Running `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-Kiss2014` generated HTML and literature audit outputs.

Current VRN result: DEV49_mean association P `1.01e-05`, PVE `4.83%`, score R² `0.0262`, score P `2.38e-05`; DEV59_mean association P `0.00168`, PVE `3.13%`, score R² `0.0147`, score P `0.00159`. The top-scored haplotype is `Hap7 = 1|1|1`, and the literature audit marks all three VRN spring diagnostic states plus the combined `1|1|1` diagnostic haplotype as `matched_top_haplotype`. However, `Hap7` has only 2 samples, so this is a positive algorithmic check but not strong population-statistical proof. The most stable effect is `Hap2 = 0|0|1` (VRN-D1 spring state only), which heads earlier than `Hap1 = 0|0|0` by `-3.36` days for DEV49 and `-2.98` days for DEV59.

Next action: refine the scoring interpretation so exact literature haplotype recovery, sample count, and phenotype direction are weighted more strongly than small heterozygous/ambiguous top-ranked haplotypes; then rerun TaGW2-B1, VRN-Kiss2014, and Rht-D1b as the regression panel.
