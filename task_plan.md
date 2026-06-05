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
Status: in_progress

Choose the next validation route: obtain the exact qHKW1 8.9-kb indel marker from a paper-specific source/contact, or move to a rice/wheat target with confirmed sample-level genotype and phenotype files.

Current decision: keep qHKW1/ZmBAM1d as the first priority, but require the paper genotype marker before interpreting the result. Added a strict preparation script for the MaizeGo full `SV.386014.zip` package. It scans only HapMap/SV matrix files for a chr1 paper-window 8.5-9.5 kb candidate and builds the database only if exactly one marker is found. With currently available small MaizeGo matrices, the script finds zero exact candidates, so qHKW1 remains data-blocked rather than method-negative.

### Phase 8: Paper-genotype rerun after full MaizeGo SV package
Status: pending

Place `SV.386014.zip` at `external_data/maize_natgenet_2019/maizego/SV.386014.zip`, then run:

```bash
python prepare_maize2019_qhkw1_paper_genotype.py
python run_star_gene_validation.py --run-analysis --paper maize2019 --target qHKW1
```

Only if the preparation script finds the paper-compatible 8.9 kb indel should the qHKW1 run be interpreted as a positive-control validation.

## Open Questions
- Which validation depth should be prioritized: full-paper datasets, target-region VCFs, or marker/phenotype positive-control tables?
- Whether the user wants local large downloads on this machine or only scripts/instructions.
- For maize, the exact paper 8.9-kb qHKW1 marker was not identified in the small MaizeGo resources; the count-compatible marker is not sufficient as proof.
- `SV.386014.zip` exists on a Baidu share, but direct automated download is blocked by Baidu verification/login mechanics from this environment.

## Errors Encountered
- Figshare API and Science article pages returned HTTP 403 from local PowerShell requests. Used accessible public pages/search snippets and official article/data URLs instead.
- Earlham WatSeq root and guessed subdirectory URLs returned HTTP 404. Paper confirms the OpenData path, but exact downloadable VCF object names should be obtained from WWWG2B/Earlham metadata or portal interaction before hard-coding.
