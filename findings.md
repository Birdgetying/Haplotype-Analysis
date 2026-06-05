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
