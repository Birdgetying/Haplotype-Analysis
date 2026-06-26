# Star Gene Validation Findings

## Repository Context
- Existing framework files: `star_gene_validation.py`, `run_star_gene_validation.py`, `star_gene_manifest.json`.
- Paper wrappers already exist for rice2024, maize2019, and wheat2024.
- Default validation is check-only and conservative.
- Manifest currently records data sources and target loci, but many coordinates and concrete local file paths are unresolved.
- Worktree is dirty with many pre-existing modified/generated files. Future commits should stage only files changed for this task.

## Data Source Research
- Rice Science 2024 (`10.1126/science.adm8762`): Data availability points the 18K-rice genotype matrix to Figshare `10.6084/m9.figshare.19166475`, title `NAM_variations`, total download size shown as 12.64 GB, license shown as CC BY 4.0. Raw sequencing projects are listed separately and are not a practical first-pass download.
- Rice Science 2024 positive-control traits need correction before implementation: the abstract reports `OsMADS22` for panicle number and `OsFTL1` for heading date. The current manifest lists heading date and plant height for `OsMADS22`, which is likely wrong for this validation.
- Maize Nature Genetics 2019 (`10.1038/s41588-019-0427-6`): Paper data availability points to NCBI/GSA/CNGB project accessions and MaizeGo Resources. MaizeGo page has direct small SV map downloads: `download/all.sv/eqtl.psv.19707.zip` (1.88 MB), `GWAS.pSV.21081.zip` (2.96 MB), `pSV.80614.zip` (8.63 MB), `SV.382254.zip` (31.37 MB), and a larger Baidu `SV.386014` (102.42 MB). It also has `download/blup_traits_final.csv` (117 KB) for agronomic BLUP traits.
- Current HTTP checks confirmed the MaizeGo Resources HTML lists the expected SV/pSV files and `blup_traits_final.csv`; HEAD checks for `SV.382254.zip` and `blup_traits_final.csv` returned content lengths consistent with the page sizes.
- MaizeGo Resources page states that some data may require contacting Dr. Jianbing Yan for further global analysis or specific gene/region analysis. This matters for qHKW1/ZmBAM1d if the direct SV map files do not include sample-level allele calls for the 8.9 kb indel.
- Wheat Nature 2024 (`10.1038/s41586-024-07682-9`): Data availability states whole-genome sequencing is at NGDC/GSA BioProject `PRJCA019636`, accession `CRA012590`. Variation matrices, annotations, wheat HapMap, phenotyping data, genetic maps, association results, tagSNPs and KASP markers are in WWWG2B. IBSpy tables, haplotypes, long-range tilling paths, VCF files and raw phenotype data are available via WWWG2B and Earlham OpenData path `https://opendata.earlham.ac.uk/wheat/under_license/toronto/WatSeq_2023-09-15_landrace_modern_Variation_Data/`.
- Current `python run_star_gene_validation.py --check-only --all --no-download` selects six targets and stops safely with `missing_required_files` plus `coordinates_unresolved` for all targets.

## Design Notes
- Existing analyzer supports loading precomputed per-target databases. This is the safest adaptation interface for heterogeneous paper data: convert each paper's available genotype/marker matrix into `gene_info.json`, `haplotype_data.csv`, `haplotype_samples.csv`, `phenotype_data.csv`, and `variant_info.csv`, then call the existing analyzer and scorer.
- Implemented design: `star_gene_data.py` stores lightweight data metadata and generates safe small-file download commands. `build_star_gene_database.py` converts a confirmed small sample-by-marker table into the existing database format.

## Maize qHKW1 / ZmBAM1d First-Pass Result
- Downloaded and inspected MaizeGo small files under `external_data/maize_natgenet_2019/maizego/` (ignored by git): `SV.382254`, `pSV.80614`, `GWAS.pSV.21081`, and `blup_traits_final.csv`.
- `blup_traits_final.csv` contains `100grainweight` for 508 rows; three samples (`BY843`, `SY1032`, `SY1035`) use `-999` / `-999.0` as missing values and must be excluded before any association or scoring.
- Paper text confirms the biological positive control: qHKW1 maps to chr1 between M49 and M40 (~177 kb), candidate gene `ZmBAM1d`, and an 8.9-kb insertion is reported as positively associated with HKW (`n=261` for 0 kb, `n=170` for 8.9 kb, `P < 0.05`).
- In the public MaizeGo matrices, a direct query of the paper figure window `chr1:30.44-30.54 Mb` found no sample-level records in `pSV_MS`, `SV_MS`, or `GWAS_pSV_MS`. The public files may use transformed/assembly-specific coordinates or may not expose the exact paper 8.9-kb marker in the small resources.
- The marker `chr1_27976568_27994848_deletion_18280` in the MS pSV/SV matrices has raw genotype counts close to the paper figure (`NN=261`, `TT=169`, plus `AA=91`) and was therefore used as the most count-compatible candidate marker for a reproducible first-pass run.
- After filtering `-999` phenotypes and retaining `NN` as a valid SV state, this count-compatible marker did not reproduce the paper association in MaizeGo BLUP `100grainweight`: three haplotypes (`NN=246`, `TT=154`, `AA=88`), corrected association P value `0.3914`, PVE `0.47%`, HaplotypeScorer regression `R^2=0.0008`, score-regression P value `0.5303`, low confidence.
- Conclusion for current maize pass: the framework can now run the public small-marker data end to end, but this marker does not prove haplotype scoring effectiveness. We should not claim a qHKW1 positive validation until the exact paper 8.9-kb indel sample-level marker is confirmed or another star gene with a confirmed public marker is used.

## Maize qHKW1 Paper-Genotype Blocker
- The Nature supplement files available locally do not include a sample-level genotype table for the Fig. 4h 8.9-kb qHKW1/ZmBAM1d indel. Supplementary Table 16 is primer/fine-mapping information, not per-accession genotype calls.
- The MaizeGo Resources page exposes `SV.386014` as a Baidu share at `https://pan.baidu.com/s/10ieQpWGTEC805K4sI4RHOg`. The Baidu page metadata confirms `SV.386014.zip`, `fs_id=1016053090657798`, `share_uk=2047402568`, `shareid=2759250563`, and size `107,394,066` bytes. However, direct automated download attempts are blocked by Baidu verification/login mechanics; guessed normal MaizeGo HTTP URLs for `SV.386014.zip` return 404.
- Added `prepare_maize2019_qhkw1_paper_genotype.py` to make the next rerun reproducible. It looks for `external_data/maize_natgenet_2019/maizego/SV.386014.zip`, extracts it when present, scans only MaizeGo/HapMap matrix files, writes `qHKW1_paper_8p9kb_candidates.tsv`, and builds `star_gene_database/maize_natgenet_2019/qHKW1` only if exactly one chr1 paper-window 8.5-9.5 kb candidate marker is found.
- Running that strict preparation script against the currently available small MaizeGo matrices found zero exact qHKW1 paper-window 8.5-9.5 kb candidates. Therefore the current qHKW1 status is data-blocked, not evidence that haplotype scoring cannot recover a real star gene.
- Interpretation rule for this validation: a negative result from an unmatched or count-compatible substitute marker must not be used to judge the method. The positive-control test only starts after the sample-level paper marker is present and verified.

## Maize SV.386014 Full-Package Rerun
- User-provided `SV.386014.zip` is present at `external_data/maize_natgenet_2019/maizego/SV.386014.zip` with size `107,394,066` bytes.
- Extracted full-package files are `SV.386014/SV.386014/svs.final.ms.txt` (`221,334,895` bytes) and `SV.386014/SV.386014/svs.final.bs.txt` (`192,796,318` bytes).
- These two files are no-header 14-column SV catalogues. Example fields are chromosome, start, end, variant type, SV length, strand, anchor id, anchor start/end, three numeric fields, score, and sequence. They do not have sample/accession genotype columns.
- Re-running `python prepare_maize2019_qhkw1_paper_genotype.py --force-extract` after the full package was present found 6 candidate matrix files, 2 catalogue files, 0 sample-level chr1 `30.44-30.54 Mb` 8.5-9.5 kb qHKW1 candidates, and 0 catalogue candidates in the same paper window.
- Expanded catalogue scan over chr1 `27-31 Mb` with 8.5-9.5 kb SV length found 8 MS records and 0 BS records; all are catalogue intervals without per-accession genotypes. Chr1-wide 8.5-9.5 kb SV counts are 269 in `svs.final.ms.txt` and 258 in `svs.final.bs.txt`.
- Code now writes an additional diagnostic TSV, `qHKW1_paper_8p9kb_catalogue_candidates.tsv`, for no-sample catalogue candidates. Catalogue hits are intentionally not accepted for database construction.
- Current conclusion is unchanged but stronger: qHKW1/ZmBAM1d remains data-blocked because the downloaded full MaizeGo package is not the needed sample-level paper-marker genotype table.

## Maize qHKW1 Paper-Marker Acquisition Route
- Nature Genetics article data availability lists raw/resequencing data under PRJNA531553, CRA001363, and CNP0000418, and states that the SV map plus Supplementary Fig. 9 step results are available through MaizeGo Resources. The article supplementary downloads are MOESM1-6 only; none is a source-data table for Fig. 4h.
- Downloaded the remaining Nature supplementary files MOESM2, MOESM4, and MOESM5. They are reporting summary, Iso-Seq/RNA-Seq summary, and eQTL enrichment summary; they do not contain qHKW1 per-accession genotypes.
- The MaizeGo Resources page explicitly provides contact `yjianbing@mail.hzau.edu.cn` and says the group can share data for specific gene or region analysis. This is the best route for the exact Fig. 4h qHKW1 8.9-kb indel genotype table.
- Public MaizeGo sample-level matrices do contain qHKW1-nearby 8.5-9.5 kb SV rows. The strongest candidate is `chr1_28599370_28608270_deletion_8900` in `SV.382254/SV.382254/MS.step1.org.txt`, with 521 sample states `TT=272` and `NN=249`; after overlap with local BLUP HKW phenotypes excluding `-999`, the groups are `TT=253` and `NN=235`, Welch t-test `P=0.15997`.
- This public row may represent a nearby 8.9-kb SV but does not match the paper Fig. 4h group sizes (`0 kb n=261`, `8.9 kb n=170`) and should not be treated as the exact positive-control genotype table without author confirmation.
- Wrote ignored local diagnostic files `external_data/maize_natgenet_2019/qHKW1_nearby_8p5_9p5kb_public_marker_summary.tsv` and `external_data/maize_natgenet_2019/qHKW1_nearby_8p5_9p5kb_public_genotypes.tsv` for manual review of public nearby candidates.

## Wheat Q7B-PH Figure 3g Positive Control
- WWWG2B provides small OneDrive-backed source files through `get_download_url_form_onedrive`. Useful file IDs are Figure 3d Q7B-PH LOD scores `01SSKBI2HAVQT54U7YRJCKGAE46A3GT5ML`, Figure 3g NIL Q7B-PH field data `01SSKBI2GUDON5GUU65VH2NJUWF3IRSKM7`, and Watkins JIC phenotype workbook `01SSKBI2HKPZS7HSLDMJG37QBQDV7JR6AJ`.
- Added `prepare_wheat2024_q7b_ph_figure3g.py`, which converts the Figure 3g workbook into a precomputed Q7B-HT database using plot-level rows by default. It writes 139 plot-level `SampleID`s such as `WL0019__plot0001`, marker `Q7B-PH_allele`, phenotype `PH_M_cm`, and `gene_info.json` source `wwwg2b_q7b_ph_figure3g`.
- The Figure 3g source has 139 complete plot-level observations and three paper-defined alleles. Group means from the analysis are Hap1/Par `n=62`, mean `101.161 cm`; Hap2/P `n=46`, mean `86.439 cm`; Hap3/W `n=31`, mean `92.994 cm`.
- Running `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Q7B-HT` produced a usable positive-control result: corrected association p-value `0.0240306`, PVE `29.63%`, HaplotypeScorer regression `R^2=0.2627`, regression p-value `1.119e-10`, high confidence, and expected-direction consistency for the height-increasing haplotype.
- Effect table confirms the top/reference haplotype is height-increasing: Hap2 effect `-14.722 cm` vs Hap1, and Hap3 effect `-8.168 cm` vs Hap1, both significant in the generated effect table.
- Interpretation limit: this validates scoring on a paper-defined haplotype label from source data. It does not yet prove the method can discover Q7B-PH from an unlabelled chr7B VCF scan. A stronger next wheat step requires confirming the physical Q7B-PH interval and extracting chr7B SNP/INDEL VCF calls.

## Rice Science 2024 Current Blocker
- Local `external_data/rice_science_2024` has no files available for analysis.
- On 2026-06-05, local PowerShell requests to `https://api.figshare.com/v2/articles/19166475`, `https://api.figshare.com/v2/articles/19166475?private_link=12978737918eecb74903`, and `https://figshare.com/ndownloader/articles/19166475?private_link=12978737918eecb74903` all returned HTTP 403.
- Therefore rice OsMADS22/OsFTL1 validation remains data-blocked. Do not interpret rice targets until `NAM_variations` or another official per-accession genotype plus phenotype source is downloaded and inspected.
- User-provided local rice data is present at `D:\Desktop\data\水稻`. It includes chromosome-wise `rice4k_chr01-12_geno.txt.gz` genotype tables, raw genotype tables, VCFs, SnpEff/Coovar/PolyPhen/SIFT H5 annotations, and phenotype tables with `Heading_date`, `Plant_height`, `Grain_weight`, `Grain_length`, `Grain_width`, and related columns. The genotype header uses sample IDs such as `B001...`, and the phenotype table uses IDs such as `C001...`, which are compatible with the existing local rice scripts for overlapping samples.
- This local rice dataset can support local positive-control runs such as DEP1/GW7 through `run_rice_paper_genes.py`, but it should not yet be described as the exact Science 2024 `OsMADS22`/`OsFTL1` validation route until those target gene IDs/coordinates and their paper-specific trait definitions are resolved.

## Wheat TaGW2 Grain-Weight Target
- WheatOmics gene search for `GW2` returns grain-weight/size/width targets: `TraesCS6A02G189300` = `TaGW2-A1` at `chr6A:237.734-237.759 Mb`, `TraesCS6B02G215300` = `TaGW2-B1` at `chr6B:291.761-291.764 Mb`, and `TraesCS6D02G176900` = `TaGW2-D1` at `chr6D:175.712-175.721 Mb`.
- WheatOmics `geneDetail.py` confirms the Chinese Spring 1.0 coordinates: A1 `237,734,651-237,760,058`, B1 `291,761,229-291,778,752`, D1 `175,712,228-175,721,507`; all are annotated as Protein SIP5 with RING-type zinc-finger domains.
- WheatOmics known-gene result summarizes the cited mechanism: TaGW2-B1 and TaGW2-D1 affect TGW through grain width/length, B1 has the stronger effect, and B1/D1 modulate cell number and length in the developing grain outer pericarp.
- Local wheat phenotypes are currently not suitable for GW2 grain-weight validation: `D:\Desktop\wheat_data\Phe.raw2` has columns `PH`, `SL`, `SN`, `SGN`, `FL`, `FW`, `NBS`, `AL` but no confirmed TGW/HKW/grain-weight column; `Phe.txt` has `TFW_DSI`, not grain weight.
- Coordinate caution: local VCF/GFF are CS-IAAS-style `Chr*` coordinates. A quick CS-IAAS GFF check near the WheatOmics Chinese Spring 1.0 GW2 coordinates found `CSIAAS6AG0515400HC` near A1, two small `CSIAAS6BG06472/47400HC` genes near B1, and no gene near D1. Do not claim a direct CS-IAAS GW2 run until the TraesCS-to-CSIAAS coordinate mapping is resolved or a WheatOmics/IWGSC-coordinate genotype file is downloaded.
- User-provided files under `D:\Desktop\data\GW2` were inspected on 2026-06-06. The three target directories contain `Sequence_download.txt` files for `TraesCS6A02G189300`, `TraesCS6B02G215300`, and `TraesCS6D02G176900`; these are FASTA-style Chinese Spring 1.0 genomic/gene/CDS/protein sequence downloads, not sample-level genotype matrices or VCFs. `TraesCS6A02G189300\chart.png` is an expression plot, not a genotype or phenotype table.
- Therefore the current `D:\Desktop\data\GW2` download is useful as gene metadata but is not sufficient for TaGW2 haplotype scoring or grain-weight association. Required missing inputs are a per-accession genotype matrix/VCF for the TaGW2 regions in a matching coordinate system and a same-sample TGW/HKW/grain-weight phenotype table.

## Wheat TaGW2 A/D WatSeq VCF Run
- User-provided WatSeq chr6A and chr6D SNP VCFs are present under `D:\Desktop\data\GW2`, with matching `.csi` indexes. chr6B VCF is still not locally available because download speed/failures block the file, so the current run is explicitly A/D-only and must not be presented as complete TaGW2 validation.
- Added `prepare_wheat2024_tagw2_ad.py`, which builds precomputed databases from the local chr6A/chr6D VCFs plus the Watkins JIC phenotype workbook `Watkins_Collection_WGIN_WISP_DFW_watseq_phenotype_data_JIC.xlsx`. The trait used for validation is `TGW_mean`, computed from available CFLN10 and CFLN14 thousand-grain-weight columns.
- `TaGW2-A1` database source is `watseq_tagw2_ad_vcf`, region `6A:237,732,651-237,760,058` including 2 kb promoter, 709 phenotype-overlap samples, 55 SNPs, and 19 retained haplotypes. Variant annotations from the precomputed metadata are 9 promoter SNPs and 46 intronic SNPs.
- `TaGW2-D1` database source is `watseq_tagw2_ad_vcf`, region `6D:175,710,228-175,721,507` including 2 kb promoter, 820 phenotype-overlap samples, 2 SNPs, and 3 retained haplotypes. Variant annotations are 1 promoter SNP and 1 intronic SNP.
- Running `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-A1 --target TaGW2-D1` completed both A/D targets and wrote combined summaries to `star_gene_results/validation_summary.json` and `.csv`.
- `TaGW2-A1` is the stronger positive-control result: corrected haplotype association P value `7.07e-06`, PVE `16.30%`, HaplotypeScorer regression `R^2=0.0431`, score-regression P value `2.416e-08`, confidence `moderate`. The output HTML files are `star_gene_results/wheat_nature_2024/TaGW2-A1/TaGW2-A1.html` and `haplotype_score.html`.
- `TaGW2-D1` has a detectable but weak/low-confidence scoring signal: corrected haplotype association P value `1.0`, PVE `3.95%`, HaplotypeScorer regression `R^2=0.0257`, score-regression P value `3.893e-06`, confidence `low`. The top scored haplotype has only 16 samples, and the region has only 2 retained SNPs; this should be interpreted as limited marker coverage rather than strong evidence.
- The A/D-only conclusion is: the workflow can recover a known grain-weight gene signal for TaGW2-A1 from target-region WatSeq SNPs and TGW phenotypes, while TaGW2-D1 remains weak with current SNP coverage. The biologically important `TaGW2-B1` is still missing and should be added when chr6B VCF download succeeds.

## Wheat TaGW2-A1 Hap8 Functional-Marker Diagnosis
- Added `analyze_tagw2_hap8_functional_groups.py` to separate the current TaGW2-A1 Hap8 signal into a published promoter-SNP component and a Hap8 background-SNP component. Outputs are under `star_gene_results/wheat_nature_2024/TaGW2-A1/tagw2_hap8_*`.
- Coordinate conversion used the local WheatOmics sequence header for `TraesCS6A02G189300`: gene start `237,734,651`, CDS offset `185`, so ATG is `237,734,835`. This maps published promoter offsets `-739` and `-593` to `chr6A:237,734,096` and `chr6A:237,734,242`; both are present in the current SNP database and are Hap8-defining relative to Hap1.
- Hap8 is not one-to-one with the published two-SNP promoter state. The local `-739/-593` Hap8-like pair is present in 66 Watkins samples across `Hap7=29`, `Hap8=25`, `Hap15=8`, and `Hap17=4`. Full Hap8 is only the 25-sample background subset.
- Trait split shows the stronger negative TGW signal comes from the Hap8 background SNPs, not just the two published promoter SNPs: `PublishedPair` vs others has `TGW_mean` effect `-3.18 g`, Welch `P=7.99e-04`; `Hap8Background/FullHap8` vs others has effect `-7.55 g`, Welch `P=3.36e-07`; `Published+Background` vs `PublishedOnly` has effect `-7.09 g`, Welch `P=1.87e-05`.
- Country composition is non-random: the two-SNP published pair is enriched in China (`32/66`; Fisher `P=1.01e-15`, BH FDR `3.23e-14`), and Hap8 itself is mostly China (`18/25`). After within-country centering, the Hap8 background signal remains negative for `TGW_mean` with effect `-4.74 g`, Welch `P=1.33e-04`, so country origin explains part but not all of the observed Hap8 effect.
- The current SNP-only database cannot confirm the later Jaiswal four-SNP high-TGW promoter haplotype because `SNP-494` maps to `chr6A:237,734,341` and is absent from `variant_info.csv`. It also cannot test the reported `G2373A` splice mutation or coding indel/SV variants without an INDEL/SV VCF and coordinate reconciliation.
- Local `D:\Desktop\data\GW2` has no chr6A INDEL/SV VCF. The Earlham chr6A INDEL VCF URL is reachable by HEAD with `Content-Length: 8,519,486,608` bytes and its `.csi` is `462,684` bytes, but a range HEAD test returned `HTTP 200` rather than `206`; this machine also lacks `pysam`, `bcftools`, and `tabix`. Download and extraction instructions were written to `tagw2_hap8_indel_sv_download_instructions.md`.

## Wheat TaGW2-A1 Remote SNP Check
- On 2026-06-14, WSL `bcftools 1.19` and `tabix 1.19` are available locally.
- Manual `bcftools` region queries must use the true contig name `chr6A`. The earlier `6A` form gives empty output against the local chr6A VCF even for known present variants.
- Local filtered WatSeq chr6A SNP VCF has `SNP-988` (`chr6A:237,733,847 G>A`), `SNP-739` (`chr6A:237,734,096 G>A`), and `SNP-593` (`chr6A:237,734,242 A>G`), but not Jaiswal `SNP-494` (`chr6A:237,734,341 G>A`).
- Remote WheatOmics merged SNP VCF does contain `SNP-494`: `chr6A 237734341 chr6A_237734341 G A`, annotated as `c.-494G>A` upstream of the TaGW2-A1 transcript. It also contains nearby promoter SNPs at `237734277` (`c.-558C>T`) and `237734327` (`c.-508T>C`).
- The remote merged SNP VCF has 1051 samples, including 827 `WATDE*` samples, so it can potentially rebuild a TaGW2-A1 SNP-only positive-control database with the Watkins TGW phenotype workbook without downloading the full 24 GB VCF.
- A micro-VCF for the TaGW2-A1 region was extracted to `external_data/wheat_nature_2024/tagw2_remote/TaGW2-A1.remote_wheatomics_237732651_237760058.vcf.gz`; it has 192 records and includes `SNP-494`.
- Current interpretation: the existing TaGW2-A1 HTML cannot fully validate the published Jaiswal functional haplotype because the local filtered VCF omitted `SNP-494`; a remote-SNP rerun is the next stronger GW2 test.
- Earlham chr6A INDEL endpoint currently returns `Request Rejected` HTML to range/body requests from this machine, so INDEL/SV augmentation is still blocked unless the endpoint becomes healthy or the file is downloaded manually/HPC-side.

## Literature Functional-Variant Audit
- Added generic manifest-driven literature audit outputs. For targets with `literature_variants` or `literature_haplotypes`, analysis now writes `literature_variant_audit.csv` and `literature_variant_audit.json` in the target result directory.
- Audit statuses are intentionally strict. A known functional marker must be present in the database, must segregate in the current analysis samples, and the top-scored haplotype must contain the expected allele/pattern before it is called `matched_top_haplotype`.
- The local WatSeq `TaGW2-A1` audit has five records. The first three Jaiswal promoter SNPs are matched by the top-scored haplotype, but `SNP-494` is `missing_from_database`; therefore the full `Jaiswal_Hap5` row is also `missing_from_database`.
- The remote-SNP `TaGW2-A1-remoteSNP` audit also has five records. It covers `SNP-494`, but the current 552 complete Watkins samples are all `G` at `chr6A:237,734,341`, so the expected `A` allele has `carrier_count=0`, `segregating_in_current_samples=False`, and `validation_status=present_not_segregating`. The full `Jaiswal_Hap5` pattern is `haplotype_not_observed`.
- Interpretation: the top-scored TaGW2-A1 haplotypes agree with the first three promoter SNPs, but current data still cannot prove recovery of the published Jaiswal causal `SNP-494` or full Hap5 because the expected `A` allele is absent among analyzed samples.

## Wheat Teacher Star-Gene Audit
- Teacher-provided targets were triaged in `wheat_star_gene_positive_control_audit.md`: Rht1, VRN, GW2, and TaPIF4.
- Rht1 exact functional SNPs were confirmed from WheatOmics remote merged SNP VCF annotations: `Rht-B1b` is `chr4B:30861571 C>T`, `c.190C>T`, `p.Gln64*`; `Rht-D1b` is `chr4D:18781242 G>T`, `c.181G>T`, `p.Glu61*`. These are stop-gained DELLA variants and are better positive-control markers than the earlier broad Rht-A1 region result.
- Added `prepare_wheat2024_rht1_functional_snps.py`. It builds precomputed databases from the exact Rht-B1b/Rht-D1b SNPs and Watkins CFLN06 plant-height phenotypes, and writes `external_data/wheat_nature_2024/rht1_functional/prepare_status.json`. It continues when one target is blocked instead of aborting the whole run.
- Real data status: `Rht-B1b` is blocked because the expected T allele has 0 carriers in Watkins phenotype-overlap samples. `Rht-D1b` builds with 802 samples and 1 functional SNP, but only 5 carriers.
- `Rht-D1b` analysis result: top-scored haplotype `Hap2` contains the literature T allele and is shorter (`Hap2` n=3 mean `84.58 cm`; `Hap1` G n=797 mean `99.62 cm`; `Hap3` G/T n=2 mean `90.50 cm`). The audit status is `matched_top_haplotype`, direction is consistent, but confidence is low (`R^2=0.006`, score regression `P=0.0288`, `PVE=0.61%`), so it is a weak positive control.
- VRN status: the Watkins workbook contains growth habit column `GrwHabit_E_sw-CFLN06`, but VRN validation needs diagnostic structural/CNV marker genotypes for `VRN-A1/VRN-B1/VRN-D1`. SNP-only region data are not enough because key VRN alleles are often promoter/intron-1 deletions or copy-number variants.
- TaPIF4 status: the public `QinZhen1995/CAU-TaSG` repository contains scripts under `03.PIF4_Promoter_Hap` and `samlist`, but not the final `coverage.martix`, `PIF_hap.txt`, or matching TGW/height phenotype table. It is currently not suitable as a clean positive control.
- Current ranking for proving scoring effectiveness among teacher targets: `Rht-D1b` is usable but weak; `GW2-A1` is gene-level only and not exact-Hap5 validated; `VRN` is promising but data-blocked; `TaPIF4` is unsuitable without final haplotype/phenotype files.

## Wheat TaGW2-B1 Functional-Haplotype Follow-up
- Qin et al. 2014 BMC Plant Biology Fig. 2B reports `TaGW2-6B` natural promoter haplotypes; the favorable `Hap-6B-1` diagnostic states are `-1709=A`, `-721=G`, and `-83=C`.
- The local WheatOmics `TraesCS6B02G215300` sequence header has `chr6B:291,761,229-291,778,752` and `CDS=169`. Therefore ATG is `chr6B:291,761,397`. The first coordinate pass used `291,759,688/291,760,676/291,761,314`, but the complete WWWG2B 1047 chr6B SNP VCF contains the corresponding records at `chr6B:291,759,689 C>A`, `chr6B:291,760,677 G>A`, and `chr6B:291,761,315 T>C`. The manifest and preparation script now use these VCF coordinates.
- Added `download_wwwg2b_file.py`, a resumable WWWG2B downloader that refreshes the OneDrive temporary URL, uses `curl -C -`, and moves HTML/non-gzip responses to `.invalid`. It downloaded the complete `SNP_matrix_1047` chr6B VCF to `D:\Desktop\data\GW2\chr6B.HARD.SNP.Missing-unphasing.ID.ann.finalSID.1047.allele2_retain.hard_retain.InbreedingCoeff_retain.vcf.gz` (`2,117,273,732` bytes).
- Updated `prepare_wheat2024_tagw2_b1_remote_snp.py` to use indexed `bcftools` region extraction for local `.vcf.gz` files when `.csi`/`.tbi` is present, avoiding a slow linear scan through the 2.1 GB VCF.
- Running `python prepare_wheat2024_tagw2_b1_remote_snp.py --vcf D:\Desktop\data\GW2\chr6B.HARD.SNP.Missing-unphasing.ID.ann.finalSID.1047.allele2_retain.hard_retain.InbreedingCoeff_retain.vcf.gz --min-haplotype-count 1` built `star_gene_database/wheat_nature_2024/TaGW2-B1-remoteSNP` with 816 phenotype-overlap samples, 3 variants, and 6 haplotypes.
- Running `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP` generated `star_gene_results/wheat_nature_2024/TaGW2-B1-remoteSNP/TaGW2-B1-remoteSNP.html`, `haplotype_score.html`, and `literature_variant_audit.csv`.
- TaGW2-B1 result: the exact Qin2014 `Hap-6B-1` pattern `A|G|C` is `Hap1` with 371 samples and the highest TGW mean (`40.970 g`). `C|G|C` is lower (`Hap2`, n=364, `39.049 g`), and `C|G|T` is much lower (`Hap3`, n=66, `34.997 g`), so the known functional haplotype has the expected favorable phenotypic direction.
- The default HaplotypeScorer top haplotype is not the exact Qin2014 haplotype: top score is `Hap5 = C/A|G|C` (n=7, score `4.2581`, TGW mean `39.610 g`). It contains the expected alleles under a loose allele-membership check, but strict audit marks the full literature haplotype as `contained_in_top_haplotype_not_exact`; exact matching haplotypes are `Hap1:371`.
- Interpretation: TaGW2-B1 is the strongest positive-control dataset so far for confirming biological direction, but it also exposes an algorithmic issue. The current scoring over-rewards a tiny heterozygous/ambiguous haplotype relative to the large exact literature favorable haplotype, so evidence supports improving the scoring/audit criteria before using top rank alone to nominate unknown genes.
- `score_mode=robust_discovery` fixes this positive-control failure without using literature labels in the score. The mode adds sample-size reliability and ambiguous-token penalties, reduces duplicate rare/special-signal weights, and keeps `default` unchanged. On TaGW2-B1, robust mode ranks exact Qin2014 `Hap1=A|G|C` first (`n=371`, score `1.8132`) and drops rare ambiguous `Hap5=C/A|G|C` to third (`n=7`, score `0.9733`, reliability `0.2593`, ambiguity `0.7333`). The robust audit status for `Qin2014_Hap-6B-1` is `matched_top_haplotype`.
- 2026-06-23 full-region rerun: `prepare_wheat2024_tagw2_b1_remote_snp.py` now extracts all polymorphic SNPs from the complete chr6B VCF over `chr6B:291759689-291778752`, while still requiring the three Qin2014 promoter SNPs to be present. The rebuilt database has 756 phenotype-overlap samples, 71 regional SNPs, and 112 haplotypes; annotation counts are 8 promoter SNPs and 63 intronic SNPs. The regenerated integrated HTML shows 23 variants because the main plot only displays positions variable among the top displayed haplotypes.
- Full-region interpretation: the Qin2014 `A/G/C` promoter pattern remains covered and segregating, but 68 additional background SNPs split that promoter class into many exact regional haplotypes. Therefore the full-region strict literature audit reports `present_but_not_top` in both default and robust runs. This does not mean the complete VCF failed; it means exact functional promoter validation and full-region discovery ranking answer different questions.

## Wheat VRN External Positive-Control Route
- Watkins has `GrwHabit_E_sw-CFLN06`, but no local diagnostic `VRN-A1/VRN-B1/VRN-D1` genotype table.
- Kiss et al. 2014 is a promising external route because it used diagnostic molecular markers for `VRN-A1`, `VRN-B1`, `VRN-D1`, `PPD-B1`, and `PPD-D1` in 683 wheat genotypes with heading phenotypes.
- The first local `Kiss2014_MOESM1.doc` was incomplete (`3,133,440` bytes) and `olefile` reported an incomplete OLE sector. The correct Springer static file is `11032_2014_34_MOESM1_ESM.doc`, `8,219,648` bytes, downloaded from `https://static-content.springer.com/esm/art%3A10.1007%2Fs11032-014-0034-2/MediaObjects/11032_2014_34_MOESM1_ESM.doc`; the server supports byte ranges, so `curl -L -C -` works for resume.
- The complete DOC contains an embedded Excel Workbook stream. Extracted workbook path: `external_data/literature/vrn_kiss2014/Kiss2014_embedded_workbook.xls`.
- The embedded workbook has per-genotype tables `stepwise_reg2011_2012a` and `stepwise_reg2011_2012b`. After filtering, there are 676 complete genotypes with diagnostic marker states for `VRN-A1`, `VRN-B1`, `VRN-D1`, `PPD-D1`, and `PPDB1`, plus 2011/2012 DEV49 and DEV59 heading dates.
- Added `prepare_wheat2024_vrn_kiss2014.py`, manifest target `VRN-Kiss2014`, and tests. The default database uses only the three VRN markers (`VRN-A1/VRN-B1/VRN-D1`) as the haplotype panel and keeps PPD markers in intermediate summaries to avoid mixing photoperiod effects into the main VRN positive control.
- Running `python prepare_wheat2024_vrn_kiss2014.py` built `star_gene_database/wheat_nature_2024/VRN-Kiss2014` with 676 samples, 3 diagnostic markers, and 7 VRN haplotypes: `Hap1=0|0|0` n=575, `Hap2=0|0|1` n=31, `Hap3=0|1|0` n=27, `Hap4=1|0|0` n=24, `Hap5=1|1|0` n=12, `Hap6=0|1|1` n=5, `Hap7=1|1|1` n=2.
- Running `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target VRN-Kiss2014` produced `star_gene_results/wheat_nature_2024/VRN-Kiss2014/VRN-Kiss2014.html` and `haplotype_score.html`. DEV49_mean association P value is `1.01e-05`, PVE `4.83%`, score regression `R^2=0.0262`, P `2.38e-05`, low confidence. DEV59_mean association P value is `0.00168`, PVE `3.13%`, score regression `R^2=0.0147`, P `0.00159`, low confidence.
- Literature audit now has four matched rows for VRN: all three spring diagnostic markers (`VRN-A1=1`, `VRN-B1=1`, `VRN-D1=1`) and the combined `1|1|1` diagnostic haplotype are covered, segregating, and present in the top-scored haplotype `Hap7`. However, `Hap7` has only 2 samples, so it is a positive algorithmic check but not a strong population-statistical proof.
- The most stable VRN phenotypic effect is `Hap2=0|0|1` (VRN-D1 spring state only), not the rare all-three spring combination. `Hap2` is significantly earlier than `Hap1=0|0|0`: DEV49_mean effect `-3.36` days, P `0.0009`; DEV59_mean effect `-2.98` days, P `0.0022`.

## Star-Gene Decisive-Variant Checklist
- Added the checklist to `star_gene_validation_record.md` so each teacher-requested target now has an explicit literature variant/haplotype, current marker/data status, top-score match status, and validation verdict.
- Literature/source checks used Crossref metadata for Qin et al. 2014 (`TaGW2-6B`, `doi:10.1186/1471-2229-14-107`), Jaiswal et al. 2015 (`TaGW2-A1`, `doi:10.1371/journal.pone.0129400`), Peng et al. 1999 (`Rht-1`, `doi:10.1038/22307`), Ellis et al. 2002 (`Rht-B1b/Rht-D1b` perfect markers, `doi:10.1007/s00122-002-1048-4`), Kiss et al. 2014 (`VRN/PPD`, `doi:10.1007/s11032-014-0034-2`), and Cao et al. 2024 (`TaSG-D1-TaPIF4`, `doi:10.1038/s41467-024-46419-0`).
- Current strongest proof is TaGW2-B1 under `robust_discovery`: exact Qin2014 `Hap-6B-1` (`A|G|C`) is covered, segregating, large (`n=371`), and top-ranked without using literature labels in scoring.
- Current weak proof is Rht-D1b: the exact functional stop-gained T allele is recovered by the top haplotype and has the expected shorter-height direction, but only a few carriers overlap the phenotype.
- VRN-Kiss2014 is marker-level validation, not causal-SV validation. It supports the audit workflow and shows the need to report rare-candidate rank separately from stable-common rank.
- TaGW2-A1, Rht-B1b, TaGW2-D1, and TaPIF4 do not currently prove scoring effectiveness because the decisive allele is absent, no decisive marker is configured, or final sample-level haplotype/phenotype data are missing.

## Three-Proof Positive-Control Set
- `TaGW2-B1-remoteSNP__robust_discovery` is proof #1: exact Qin2014 favorable `Hap-6B-1 = A|G|C` is the raw top-ranked robust haplotype with n=371 and audit `matched_top_haplotype`.
- `VRN-D1-Kiss2014__robust_discovery` is proof #2 at diagnostic-marker level: single-marker `VRN-D1=1` is top-ranked (`Hap2`, n=38), audit `matched_top_haplotype`, with DEV49 score regression P `2.98e-07` and DEV59 P `7.49e-05`.
- `Rht-Zanke2014__robust_discovery` is proof #3 under direction-aware validation: Zanke et al. 2014 Table S2 provides 368 varieties with Rht-B1/Rht-D1 marker states and plant-height phenotypes. Raw robust top is the high-plant-height wild-type haplotype, but the stable low-phenotype directional top is `Hap1=B1a|D1b` (n=214, BLUE mean `82.540 cm`), and the literature audit now reports `Rht-D1b_diagnostic_marker` as `matched_directional_top_haplotype`.
- The Rht result is also an algorithmic finding: for decreases-trait genes, raw "highest total score" is not necessarily the desired biological allele. The validation summary now records `directional_top_haplotype` fields, and the literature audit records `directional_validation_status`.
- TaPIF4 remains blocked. `Supplementary_Data_5_accessions_haplotypes.xlsx` contains 331 accessions and TaPIF4 haplotype labels, but no matching phenotype table. The official `MOESM11` Source Data ZIP is 82,112,605 bytes; repeated curl downloads were incomplete/corrupt in this session, so it was not used as evidence.

## 2026-06-23 Robust Scoring Cross-Target Check
- Fixed a full-region scoring bug: haplotype sequences store all `variant_info` positions, while displayed/scored positions can be filtered. The scorer now maps genomic position to full-sequence index before reading alleles, so regional haplotypes are not shifted when only a subset of positions is displayed.
- Added discovery-only robust core grouping for `robust_discovery`. Core positions are selected from phenotype signal, EB marker effect, functional annotation weight, MAF stability, missingness, and LD pruning; literature markers are not used as scoring inputs.
- Added stable direction-aware top selection for both full haplotypes and robust core groups. Low-support extreme groups are filtered when higher-support alternatives exist, preventing full-region singleton/near-singleton haplotypes from dominating the direction-aware validation summary.
- Re-ran, in one sequential command to avoid `validation_summary.csv` overwrite races:
  `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --target VRN-D1-Kiss2014 --target Rht-Zanke2014 --score-mode robust_discovery`.
- `TaGW2-B1-remoteSNP__robust_discovery`: 756 samples, 71 regional SNPs, 112 haplotypes. TGW_mean score regression remains significant (`R^2=0.0124`, `P=0.0022046`), but Qin2014 `Hap-6B-1=A|G|C` is not recovered as full-region top after background SNPs split the promoter class. This is no longer a strong functional-haplotype proof under the full-region database.
- `VRN-D1-Kiss2014__robust_discovery`: 676 samples, 1 marker, 2 haplotypes. Top and directional top are both `Hap2` carrying `VRN-D1=1` with 38 samples; DEV49_mean `R^2=0.0382`, `P=2.98e-07`; DEV59_mean `R^2=0.0230`, `P=7.49e-05`; audit is `matched_top_haplotype` and `matched_directional_top_haplotype`.
- `Rht-Zanke2014__robust_discovery`: 368 samples, 2 markers, 4 haplotypes. Raw top is the tall wild-type-like haplotype, but direction-aware validation recovers the shorter functional class: BLUE directional top `Hap1=B1a|D1b`, n=214, mean 82.540 cm, audit `Rht-D1b_diagnostic_marker=matched_directional_top_haplotype`; score regression `R^2=0.0842`, `P=1.44e-08`.
- Current proof status: the algorithmic changes are generic and validated on VRN-D1-Kiss2014 and Rht-Zanke2014, not only GW2. GW2 full-region remains a useful stress test showing the need to report raw top, direction-aware top, and core groups separately.

## 2026-06-24 Functional Sub-Haplotype Robust Scoring
- Added a generic `functional_haplotype_groups` layer for `robust_discovery`. It is selected from local annotations, phenotype marker signal, EB marker effect, MAF/missingness, LD context, and near-boundary gene-body signal; literature variants are still only used by the post hoc audit.
- Fixed the full-region scoring input path so `HaplotypeScorer` sees all regional positions instead of only the HTML display subset. This matters for TaGW2-B1, where the HTML may display 23 variants while the database has 71 retained SNPs.
- Added `functional_marker` and `diagnostic_marker` annotation support so marker-panel targets such as VRN-D1 and Rht-Zanke2014 produce functional groups instead of empty functional outputs.
- Added support-shrunk direction-aware ranking. This prevents small outlier groups from beating stable groups when their phenotype mean is only slightly more extreme.
- Re-ran the three proof targets with:
  `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target TaGW2-B1-remoteSNP --target VRN-D1-Kiss2014 --target Rht-Zanke2014 --score-mode robust_discovery`.
- Final proof status:
  TaGW2-B1 full-region is now positive through `matched_directional_top_functional_group` for Qin2014 `Hap-6B-1`; VRN-D1 remains `matched_top_haplotype`/`matched_directional_top_haplotype`; Rht-D1b remains `matched_directional_top_haplotype`.
- TaGW2-B1 limitation remains important: raw top functional group is not the Qin2014 group, and score regression is modest (`R^2=0.0078`, `P=0.0153`). This supports using direction-aware functional groups as validation evidence, not claiming that raw top score alone has solved full-region functional-haplotype discovery.

## 2026-06-25 Rht-Zanke2014 Direction Fix
- Rechecked the teacher concern that Rht-Zanke2014 Hap1/Hap3/Hap4 clustered away from Hap2 in the score plot. The haplotype partition itself matches the Zanke2014 marker panel: `Hap1=B1a|D1b` n=214, `Hap2=B1a|Rht-D1a` n=126, `Hap3=B1b|Rht-D1a` n=26, and `Hap4=B1a|D1b/Rht-D1a` n=2.
- The issue was score interpretation, not marker parsing. The retained raw robust `total` ranked the tall wild-type `Hap2` highest because effect magnitude was not direction-aware for a `decreases_trait` target.
- Added an active `directional_total` score axis when `expected_direction` is known. It preserves raw `total` for diagnosis but uses support-shrunk phenotype direction in the main HTML/summary/audit rank.
- Re-ran `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Rht-Zanke2014 --score-mode robust_discovery`. For `PlantHeight_BLUE`, top is now `Hap1=B1a|D1b`, `directional_total=0.9642`, n=214, `R^2=0.451`, `P=1.36e-49`; for `PlantHeight_GAT_2012`, top remains `Hap1`, `directional_total=0.9624`, `R^2=0.532`, `P=1.71e-60`.
- `literature_variant_audit.csv` now reports `Rht-D1b_diagnostic_marker` as both `matched_top_haplotype` and `matched_directional_top_haplotype`. `Rht-B1b` remains present but not top because that class is `Hap3=B1b|Rht-D1a`, n=26.

## 2026-06-26 Rht single-variant clarification
- Rebuilt and reran the strict single-variant `Rht-D1b` positive control with `python prepare_wheat2024_rht1_functional_snps.py --min-haplotype-count 1 --target Rht-D1b` and `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Rht-D1b --score-mode robust_discovery`.
- The exact decisive SNP remains `chr4D:18781242 G>T`, annotated `c.181G>T`, `p.Glu61*`. The rebuilt database has 802 phenotype-overlap samples and one stop-gained marker.
- New `Rht-D1b__robust_discovery` output: `Hap1=G` n=797, mean `99.615 cm`; `Hap2=T` n=3, mean `84.583 cm`, score `1.0938`; `Hap3=G/T` n=2, mean `90.500 cm`. Score regression is `R^2=0.0056`, `P=0.0334`, confidence low.
- `literature_variant_audit.csv` reports `Rht-D1b_stop_gained=matched_top_haplotype`, allele counts `G:797;G/T:2;T:3`, carrier count `5`, and `top_haplotype_exact_expected=True`.
- Support-shrunk direction-aware ranking selects the large wild-type `Hap1`, so `Rht-D1b` remains a weak exact single-SNP positive control rather than a strong proof. `Rht-Zanke2014` should be described as marker-panel reference only because it combines `Rht-B1` and `Rht-D1`.
- Fixed target selection so an exact `--target Rht-D1b` request does not also select `Rht-Zanke2014` merely because that target lists `Rht-D1b` as an alias. Exact `target_id` matches now take precedence over alias matches.

## 2026-06-26 Rht-D1-Zanke2014 single-marker clarification
- Extended `prepare_wheat2024_rht_zanke2014.py` with `--single-marker-targets`, producing one-marker databases for `Rht-B1-Zanke2014` and `Rht-D1-Zanke2014` from Zanke2014 Table S2.
- Added manifest target `Rht-D1-Zanke2014` so the analysis and literature audit can run without selecting the combined `Rht-B1/Rht-D1` panel.
- Ran `python prepare_wheat2024_rht_zanke2014.py --single-marker-targets` and `python run_star_gene_validation.py --run-analysis --paper wheat2024 --target Rht-D1-Zanke2014 --score-mode robust_discovery`.
- New output: `star_gene_results/wheat_nature_2024/Rht-D1-Zanke2014__robust_discovery/Rht-D1-Zanke2014.html`.
- Haplotype split is one marker only: `Hap1=D1b` n=214, `Hap2=Rht-D1a` n=152, `Hap3=D1b/Rht-D1a` n=2. This removes the old multi-gene `B1|D1` interpretation problem.
- The raw robust total top is the two-sample mixed `Hap3`, but the stable direction-aware top is exact `Hap1=D1b`, n=214. The audit marks `Rht-D1b_diagnostic_marker` as `matched_directional_top_haplotype`.
- Interpretation: this is a useful single-gene marker-level proof for Rht-D1, while the exact nucleotide SNP route remains weak but base-level evidence.

## 2026-06-26 WheatOmics INDEL supplement download
- Downloaded WheatOmics `trackList.json` and VCF directory listing from `https://wheatomics.sdau.edu.cn/jbrowse-1.12.3-release/Chinese_Spring1.0/`.
- Added and ran `download_wheat2024_indel_microvcfs.py`, which extracts target-region micro-VCFs from `WEC_INDEL_IWGSCv1.0.eff.vcf.gz`, `GBS_filtered_Indels_IWGSCv1.0.eff.vcf.gz`, and `WildEmmer10WGS_INDEL_eff.vcf.gz`.
- Outputs are under `external_data/wheat_nature_2024/indel_microvcfs`; the summary is `microvcf_status.tsv`.
- VRN-A1/VRN-B1/VRN-D1 regions had 0 INDEL records in all three public WheatOmics INDEL tracks, so the VRN causal deletion/CNV/SV blocker remains.
- Rht-B1, TaGW2-A1, TaGW2-B1, and TaGW2-D1 regions have a few INDEL records in WEC/GBS/WildEmmer tracks, but all queried sources have 0 overlap with Watkins phenotype samples. These data cannot directly strengthen the current Watkins-based validations.
- WWWG2B API currently returns HTTP 526 for available table/file table/download URL endpoints, and Earlham OpenData WatSeq raw VCF URLs still return `Request Rejected` HTML. The 1051-sample WatSeq INDEL/SV files remain inaccessible from this environment in this pass.
- Added and ran `download_wheat2024_indel_marker_annotations.py` for WheatOmics `Indel_marker_from_zhai`. It saves raw JBrowse chunks and target-window marker annotations under `external_data/wheat_nature_2024/zhai_indel_markers`.
- With the default 1 Mb flank, exact target windows plus nearby sequence contain 8 Zhai marker rows near `TaGW2-B1` and 7 near `VRN-B1`; `Rht-B1`, `Rht-D1`, `TaGW2-A1`, `TaGW2-D1`, `VRN-A1`, and `VRN-D1` have 0 marker rows. These are primer/marker annotations, not per-accession genotype calls, so they are useful for manual follow-up but cannot be fed into current haplotype scoring.

## 2026-06-26 Phenotype-Free Site-Weighted Robust Discovery
- Implemented a PostGWAS-like `site_weighted` component for `robust_discovery`.
  The component is phenotype-free: it uses annotation severity, explicit
  external evidence, MAF stability, missingness, LD pruning, and gene-structure
  context, but not the current validation phenotype, local haplotype means, or
  local p-values.
- Added `site_weights` and `site_weighting_policy` to score JSON and report
  labels so reviewers can see that `current_phenotype_used=False` for each
  weighted site.
- Root cause for the first strict `Rht-D1b` failure was annotation vocabulary,
  not absence of the variant. The real VCF/database annotation was snpEff
  `stop_gained`, while the scorer only retained rare high-impact `stop_gain`.
  The scorer now normalizes `stop_gained`, `stopgain`, and `nonsense` to
  internal `stop_gain`.
- Re-ran four robust targets:
  `TaGW2-B1-remoteSNP`, `VRN-D1-Kiss2014`, `Rht-D1b`, and `Rht-Zanke2014`.
  `Rht-D1b` now has one site-weight row at `chr4D:18781242`, normalized to
  `stop_gain`, and `Hap2=T` is recovered as the raw top score.
- Validation interpretation after this update: `VRN-D1-Kiss2014` remains the
  strongest marker-level proof; strict `Rht-D1b` is exact but weak because only
  5 T carriers overlap phenotype data; `TaGW2-B1` and `Rht-Zanke2014` are
  still best described through direction-aware/functional-group validation,
  not raw full-region top rank alone.
- Code-review follow-up: `variant_info` score/logp/p-value fields now require
  explicit external evidence marking before contributing to `site_weighted`.
  This closes the leak where an unmarked local `minus_log10_p`, `gwas_pvalue`,
  or `site_score` could previously raise a site weight despite `gwas_data`
  being correctly gated.
