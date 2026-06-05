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
