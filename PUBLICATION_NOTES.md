# Public deposit notes

Prepared for https://github.com/Jasonlee335/myproject on 2026-10-30.

## Changes from the supplementary archive

1. Copied existing R scripts and aggregate outputs into a separate directory.
   Original manuscripts, responses, source data, and supplementary archive files
   were not edited.
2. Removed computer-specific R-library overrides and replaced computer-specific
   default input/output paths with working-directory defaults. The runner still
   requires explicit `DXA_DATA_ROOT` and `DXA_PROJECT_DIR` and sets `DXA_CODE_DIR`.
3. Corrected Step 09's Step 09b source-file path to resolve from the code directory
   rather than the separately configured output directory.
4. Added a package check, checksum verifier, data requirements, citation metadata,
   and Git ignore rules for private inputs and generated person-level outputs.
5. Documented the consolidated CHARLS input dependency and multiple-imputation
   adjustment set. No model formulas, seeds, or numerical results were changed.

## Verification scope

The public package is checked for R syntax, declared dependencies, checksums, and
common credential patterns. CSV schemas were inspected: committed outputs are
model/cohort summaries, bootstrap summaries, or variable metadata, not participant
records. Aggregate CSVs are preserved byte for byte from the supplied archive.

The complete raw-data workflow was not rerun solely for this deposit. Data access
conditions still apply. `CHARLS.csv` is not an official release, and its upstream
consolidation code is not included. Dependency versions are documented by archived
sessions, not a lockfile. No repository DOI or manuscript acceptance is claimed.

No permissive license was selected on the authors' behalf. No raw/cleaned
participant records, secrets, manuscript drafts, or reviewer letters are included.
