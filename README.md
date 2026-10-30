# OSTA capacity-constrained DXA referral

R code and aggregate outputs for US/Korean DXA validation and a separate analysis
of self-reported hip fractures in CHARLS.

Repository: https://github.com/Jasonlee335/myproject

## Contents

- Steps 01-07: cleaning, survey validation, calibration-drift decomposition,
  capacity-constrained referral, CHARLS outcomes, and sensitivity analyses.
- Step 08: exploratory age checks, age standardization, and age-weight comparisons.
- Step 09: paired contrasts, age-specific AUC summaries, FRAX metadata audit,
  Figure 2, and supplementary Figure S5.
- Step 09b: supplementary Figure S4 layout, invoked by Step 09.
- `aggregate_results/`: cohort/model summaries, bootstrap summaries, model
  coefficients, archived session information, and input checksums.

No participant records, raw survey files, cleaned person-level datasets, fitted
person-level model objects, manuscript drafts, or reviewer correspondence are included.

## Reproducibility boundary

This is a public copy of the revision archive. Portability changes are documented
in `PUBLICATION_NOTES.md`. Statistical formulas, definitions, seeds, and archived
numerical estimates were not changed for this deposit.

Full re-execution requires locally obtained NHANES, KNHANES, CHARLS, and Harmonized
CHARLS inputs. Step 01 also requires the author's existing consolidated longitudinal
`CHARLS.csv`. This is not an official CHARLS download filename and is not distributed.
The repository does not currently reconstruct that consolidated input from official
raw modules alone. This upstream dependency limits end-to-end reproducibility; it
does not permit redistribution of participant data. See `DATA_REQUIREMENTS.md`.

## Requirements

The archived analyses used R 4.5.0. Required packages are listed in
`check_requirements.R`; archived session files are in `aggregate_results/`.
Exact dependency versions are not locked, so different R/package versions may
produce small numerical differences. Packages are not installed automatically.

Run from the repository root:

```r
source("check_requirements.R", encoding = "UTF-8")
```

Install missing packages through an approved CRAN mirror, then repeat the check.

## Run the workflow

Obtain source data under their providers' terms. Keep data and generated outputs
outside the Git checkout. Set paths to your local directories:

```r
Sys.setenv(
  DXA_DATA_ROOT = "/path/to/authorized_source_data",
  DXA_PROJECT_DIR = "/path/to/generated_analysis_outputs"
)
source("run_reproduction.R", encoding = "UTF-8")
```

Run from the repository root. The runner uses 1,000 bootstrap replicates and 20
imputations, and stops when required inputs or packages are unavailable. Step 01
may download missing public NHANES components; it does not obtain restricted-access
KNHANES or CHARLS data. Full execution is computationally intensive and writes
cleaned person-level data only to the local output directory.

Outputs beneath `DXA_PROJECT_DIR` include `data_clean/`, `quality_control/`,
`step2_results/` through `step7_results/`, `reviewer1_revision/analysis_results/`,
and `reviewer1_revision/figures/`. The committed `aggregate_results/` is an archived
reference, not an automatically updated copy of a new run's results.

To re-run only revision analyses after prior-stage outputs exist:

```r
Sys.setenv(
  DXA_DATA_ROOT = "/path/to/authorized_source_data",
  DXA_PROJECT_DIR = "/path/to/generated_analysis_outputs",
  DXA_CODE_DIR = normalizePath(getwd(), winslash = "/"),
  REVISION_BOOTSTRAP_REPS = "1000"
)
source("08_reviewer1_revision.R", encoding = "UTF-8")
source("09_revision_synthesis.R", encoding = "UTF-8")
```

The optional `audit_frax_inputs.R` reads variable metadata from local modules.
Review generated data and input manifests before sharing them: they can contain
participant information or local filesystem paths.

## Interpretation

- Primary imaging comparisons use common-US-reference femoral-neck osteoporosis.
  Some archived scenario labels call these Korean rows sensitivity analyses; use
  actual outcome/reference columns rather than historical labels alone.
- Local Korean-reference results are secondary local-service analyses.
- CHARLS has no DXA. Its self-reported fracture outcome is a separate prognostic
  analysis, not validation of osteoporosis status in Chinese women.
- OSTA is deterministic. C-indices are apparent in-sample estimates and do not
  establish information independent of age and weight.
- Age-allocation analyses are post hoc and exploratory. Fractional tie allocation
  attains empirical weighted capacity in expectation, not future clinical benefit.
- Model-performance intervals condition on fitted coefficients; the additional
  source-target bootstrap propagates cutoff uncertainty, not full model development.
- The multiple-imputation Cox estimate additionally adjusts for BMI, current
  smoking, and ever drinking. It is not the imputed counterpart of the primary
  four-covariate adjustment model.
- FRAX without BMD was not computed. The audit documents incomplete harmonization;
  age-weight models are not substitutes for a FRAX comparison.
- `Table15_final_go_no_go.csv` is an exploratory feasibility summary, not a
  preregistration, independent validation, or clinical decision rule.

## Verification and citation

`FILE_MANIFEST.json` records SHA-256 checksums. `SOURCE_PROVENANCE.json` maps copied
files to the archive checksums. See `PUBLICATION_NOTES.md` for verification scope
and `CITATION.cff` for authorship metadata. No Zenodo DOI has been assigned, and
this repository does not claim acceptance of the associated manuscript.
Use a commit permalink or version tag to identify the reviewed version.

## Data rights and code licensing

Data remain subject to NHANES, KNHANES, CHARLS, and Harmonized CHARLS access and
redistribution terms. No dataset license is granted here. The authors have not
selected an open-source software license for this deposit; public visibility is
not a blanket reuse license.
