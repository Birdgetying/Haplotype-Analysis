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

## Open Questions
- Which validation depth should be prioritized: full-paper datasets, target-region VCFs, or marker/phenotype positive-control tables?
- Whether the user wants local large downloads on this machine or only scripts/instructions.

## Errors Encountered
- Figshare API and Science article pages returned HTTP 403 from local PowerShell requests. Used accessible public pages/search snippets and official article/data URLs instead.
- Earlham WatSeq root and guessed subdirectory URLs returned HTTP 404. Paper confirms the OpenData path, but exact downloadable VCF object names should be obtained from WWWG2B/Earlham metadata or portal interaction before hard-coding.
