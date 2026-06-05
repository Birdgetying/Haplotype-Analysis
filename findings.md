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
