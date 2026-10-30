# Set DXA_DATA_ROOT to the directory containing the licensed source datasets.
# Set DXA_PROJECT_DIR to a writable directory for clean data and all outputs.
code_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (Sys.getenv("DXA_DATA_ROOT") == "") stop("Set DXA_DATA_ROOT before running")
if (Sys.getenv("DXA_PROJECT_DIR") == "") stop("Set DXA_PROJECT_DIR before running")
if (!dir.exists(Sys.getenv("DXA_DATA_ROOT"))) stop("DXA_DATA_ROOT does not exist")
source(file.path(code_dir, "check_requirements.R"), encoding = "UTF-8")
Sys.setenv(DXA_CODE_DIR = code_dir)
dir.create(Sys.getenv("DXA_PROJECT_DIR"), recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(Sys.getenv("DXA_PROJECT_DIR"), winslash = "/", mustWork = TRUE)
if (identical(tolower(output_dir), tolower(code_dir))) {
  stop("Keep generated analysis outputs outside the code directory")
}
Sys.setenv(STEP3_BOOTSTRAP_REPS="1000", STEP5_BOOTSTRAP_REPS="1000",
           STEP6_BOOTSTRAP_REPS="1000", STEP7_BOOTSTRAP_REPS="1000",
           STEP7_IMPUTATIONS="20", REVISION_BOOTSTRAP_REPS="1000")
scripts <- c(
  "01_data_cleaning.R", "02_descriptive_shift_analysis.R",
  "03_bidirectional_validation.R", "04_calibration_drift_decomposition.R",
  "05_capacity_constrained_referral.R", "06_charls_prospective_validation.R",
  "07_robustness_and_synthesis.R", "08_reviewer1_revision.R",
  "09_revision_synthesis.R"
)
stopifnot(all(file.exists(file.path(code_dir, scripts))))
for (file in file.path(code_dir, scripts)) {
  message("Running ", basename(file))
  source(file, encoding="UTF-8", local=new.env(parent=globalenv()))
}
