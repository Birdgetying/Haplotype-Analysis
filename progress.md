# Star Gene Validation Progress

## 2026-06-05
- Read project instructions from AGENTS.md provided by user.
- Checked git status and found a dirty worktree with many pre-existing changes.
- Inspected core validation files and docs.
- Created task tracking files for this multi-step research and implementation task.
- Checked public rice and maize data source availability.
- Identified MaizeGo direct download links suitable for lightweight first-pass validation.
- Checked wheat article data availability and confirmed WWWG2B/Earlham/NGDC data routes.
- Ran current check-only validator for all manifest targets; all targets are blocked by missing local inputs and unresolved coordinates, as expected.
- Added standard-library tests for lightweight data metadata, CLI download printing, and marker-matrix database conversion.
- Added `star_gene_data.py` and `build_star_gene_database.py`.
- Updated manifest rice positive-control traits and data-source notes.
- Ran `python -m unittest test_star_gene_data.py -v`: passed.
- Ran `python run_star_gene_validation.py --check-only --all --no-download`: passed, all targets remain safely blocked by missing local inputs/unresolved coordinates.
- Ran `python run_star_gene_validation.py --print-downloads` for maize/rice/wheat: passed.
- Ran `python run_rice_test.py`: exit code 0; one existing rice test gene had no phenotype sample overlap and two succeeded.
- Re-ran `python -m unittest test_star_gene_data.py -v`: passed 8 tests.
- Re-ran `python run_star_gene_validation.py --check-only --all --no-download`: passed; all six manifest targets were safely blocked by missing local inputs/unresolved coordinates.
- Re-ran `python run_star_gene_validation.py --print-downloads --paper maize2019`: printed the four MaizeGo small-file download commands.
- Re-ran `python run_rice_test.py`: exit code 0; status distribution was two `success` genes and one `no_phenotype_match` gene due to no overlap between haplotype sample IDs and phenotype sample IDs.
- Re-checked MaizeGo Resources and HEAD responses for two representative small files; direct MaizeGo links remain reachable from this machine.
