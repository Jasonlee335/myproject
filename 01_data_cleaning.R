options(stringsAsFactors = FALSE, scipen = 999)

# STEP 1: Reproducible cleaning for NHANES, KNHANES, and CHARLS.
# This script does not fit any outcome model.

required_packages <- c("haven", "dplyr", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    ". Run install.packages(c(",
    paste(sprintf("'%s'", missing_packages), collapse = ", "), "))."
  )
}

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tibble)
})

# -----------------------------------------------------------------------------
# 1. Paths and frozen analysis definitions
# -----------------------------------------------------------------------------

root <- Sys.getenv("DXA_DATA_ROOT", unset = file.path(getwd(), "data_raw"))
project_dir <- Sys.getenv("DXA_PROJECT_DIR", unset = file.path(root, "DXA_transportability_R"))
clean_dir <- file.path(project_dir, "data_clean")
qc_dir <- file.path(project_dir, "quality_control")
dir.create(clean_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

nh_root <- file.path(root, "nhanes_osteoporosis_screening", "raw")
kn_root <- file.path(root, "KNHANES", "raw")
charls_long_file <- file.path(root, "CHARLS.csv")
charls_raw_root <- file.path(root, "charls_official", "raw_extracted")
charls_original_root <- file.path(root, "charls原版数据")

age_min <- 50
age_max <- 79

# NHANES III young non-Hispanic White female reference values.
# Femoral neck is the primary harmonized diagnostic site.
reference_values <- tibble(
  site = c("femoral_neck", "total_hip", "lumbar_spine"),
  mean_bmd = c(0.858, 0.942, 1.040),
  sd_bmd = c(0.120, 0.122, 0.110)
)
write.csv(reference_values, file.path(qc_dir, "00_tscore_reference_values.csv"),
          row.names = FALSE, na = "")

nh_periods <- tribble(
  ~cycle, ~start_year, ~analysis_role, ~demo_file, ~bmx_file,
  ~femur_file, ~spine_file, ~rhq_file, ~weight_variable, ~weight_divisor,
  "2005-2006", 2005L, "development", "DEMO_D.XPT", "BMX_D.XPT",
  "DXXFEM_D.XPT", "DXXSPN_D.XPT", "RHQ_D.XPT", "WTMEC2YR", 3,
  "2007-2008", 2007L, "development", "DEMO_E.XPT", "BMX_E.XPT",
  "DXXFEM_E.XPT", "DXXSPN_E.XPT", "RHQ_E.XPT", "WTMEC2YR", 3,
  "2009-2010", 2009L, "development", "DEMO_F.XPT", "BMX_F.XPT",
  "DXXFEM_F.XPT", "DXXSPN_F.XPT", "RHQ_F.XPT", "WTMEC2YR", 3,
  "2013-2014", 2013L, "temporal_validation_1", "DEMO_H.XPT", "BMX_H.XPT",
  "DXXFEM_H.XPT", "DXXSPN_H.XPT", "RHQ_H.XPT", "WTMEC2YR", 1,
  "2017-March 2020", 2017L, "temporal_validation_2", "P_DEMO.XPT", "P_BMX.XPT",
  "P_DXXFEM.XPT", "P_DXXSPN.XPT", "P_RHQ.XPT", "WTMECPRP", 1
)

kn_years <- c(2008:2011, 2024L)

wave_year <- c(`1` = 2011L, `2` = 2013L, `3` = 2015L,
               `4` = 2018L, `5` = 2020L)

# -----------------------------------------------------------------------------
# 2. General cleaning and audit helpers
# -----------------------------------------------------------------------------

assert_file <- function(path) {
  if (!file.exists(path)) stop("Required input file not found: ", path)
  invisible(path)
}

as_num <- function(x) suppressWarnings(as.numeric(x))

clean_range <- function(x, lower, upper) {
  x <- as_num(x)
  x[!is.finite(x) | x < lower | x > upper] <- NA_real_
  x
}

clean_binary <- function(x, zero_values = 0, one_values = 1) {
  x <- as_num(x)
  out <- rep(NA_integer_, length(x))
  out[x %in% zero_values] <- 0L
  out[x %in% one_values] <- 1L
  out
}

clean_positive_weight <- function(x) {
  x <- as_num(x)
  x[!is.finite(x) | x <= 0] <- NA_real_
  x
}

clean_id <- function(x) {
  out <- trimws(as.character(x))
  out[out %in% c("", ".", "NA", "<NA>")] <- NA_character_
  out
}

get_column <- function(data, candidates, default = NA_real_) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) return(rep(default, nrow(data)))
  data[[hit[1]]]
}

min_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

max_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

any_or_na <- function(x) {
  if (all(is.na(x))) return(NA_integer_)
  as.integer(any(x == 1, na.rm = TRUE))
}

any_tscore_osteoporosis <- function(...) {
  scores <- cbind(...)
  apply(scores, 1, function(z) {
    if (all(is.na(z))) return(NA_integer_)
    as.integer(any(z <= -2.5, na.rm = TRUE))
  })
}

complete_tscore_osteoporosis <- function(...) {
  scores <- cbind(...)
  apply(scores, 1, function(z) {
    if (any(is.na(z))) return(NA_integer_)
    as.integer(any(z <= -2.5))
  })
}

duplicate_key_report <- function(data, keys, dataset) {
  key_data <- data %>% select(all_of(keys))
  duplicated_any <- duplicated(key_data) | duplicated(key_data, fromLast = TRUE)
  tibble(
    dataset = dataset,
    key = paste(keys, collapse = "+"),
    rows = nrow(data),
    missing_key_rows = sum(!complete.cases(key_data)),
    duplicate_key_rows = sum(duplicated_any),
    duplicate_key_rate = mean(duplicated_any)
  )
}

missingness_report <- function(data, dataset) {
  tibble(
    dataset = dataset,
    variable = names(data),
    n = nrow(data),
    missing_n = vapply(data, function(x) sum(is.na(x)), numeric(1)),
    missing_pct = vapply(data, function(x) mean(is.na(x)) * 100, numeric(1)),
    distinct_n = vapply(data, function(x) length(unique(x[!is.na(x)])), numeric(1))
  )
}

flow_row <- function(cohort, step, data, condition = rep(TRUE, nrow(data))) {
  condition[is.na(condition)] <- FALSE
  tibble(cohort = cohort, step = step, n = sum(condition))
}

file_manifest <- function(paths) {
  paths <- unique(paths[file.exists(paths)])
  info <- file.info(paths)
  tibble(
    file = normalizePath(paths, winslash = "/", mustWork = TRUE),
    bytes = info$size,
    modified = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    md5 = unname(tools::md5sum(paths))
  )
}

write_both <- function(data, stem) {
  saveRDS(data, file.path(clean_dir, paste0(stem, ".rds")), compress = "xz")
  write.csv(data, file.path(clean_dir, paste0(stem, ".csv")),
            row.names = FALSE, na = "")
}

# -----------------------------------------------------------------------------
# 3. NHANES ingestion and cleaning
# -----------------------------------------------------------------------------

download_nhanes_component <- function(start_year, filename) {
  destination <- file.path(nh_root, filename)
  if (file.exists(destination) && file.info(destination)$size > 1000) {
    return(destination)
  }
  url_filename <- sub("\\.XPT$", ".xpt", filename)
  url <- sprintf(
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%s/DataFiles/%s",
    start_year, url_filename
  )
  message("Downloading official NHANES file: ", url)
  tryCatch(
    download.file(url, destination, mode = "wb", quiet = TRUE),
    error = function(e) warning("Could not download ", filename, ": ", e$message)
  )
  if (!file.exists(destination) || file.info(destination)$size <= 1000) {
    warning("NHANES component unavailable: ", destination)
    return(NA_character_)
  }
  destination
}

nhanes_input_files <- character()
nhanes_join_qc <- list()

clean_nhanes_period <- function(spec) {
  files <- c(
    demo = download_nhanes_component(spec$start_year, spec$demo_file),
    bmx = download_nhanes_component(spec$start_year, spec$bmx_file),
    femur = download_nhanes_component(spec$start_year, spec$femur_file),
    spine = download_nhanes_component(spec$start_year, spec$spine_file),
    rhq = download_nhanes_component(spec$start_year, spec$rhq_file)
  )
  if (any(is.na(files))) stop("At least one required NHANES component is unavailable for ", spec$cycle)
  invisible(lapply(files, assert_file))
  nhanes_input_files <<- c(nhanes_input_files, unname(files))

  demo <- read_xpt(files[["demo"]])
  bmx <- read_xpt(files[["bmx"]]) %>%
    select(any_of(c("SEQN", "BMXWT", "BMXBMI")))
  femur <- read_xpt(files[["femur"]]) %>%
    select(any_of(c("SEQN", "DXXOFBMD", "DXXNKBMD", "DXXFMBCC")))
  spine <- read_xpt(files[["spine"]]) %>%
    select(any_of(c("SEQN", "DXXOSBMD", "DXXOSBCC")))
  if (!is.na(files[["rhq"]]) && file.exists(files[["rhq"]])) {
    rhq <- read_xpt(files[["rhq"]]) %>%
      select(any_of(c("SEQN", "RHQ031", "RHD042", "RHQ060", "RHD280",
                      "RHQ291", "RHQ305", "RHQ332")))
  } else {
    rhq <- tibble(SEQN = demo$SEQN)
  }

  joined <- demo %>%
    left_join(bmx, by = "SEQN") %>%
    left_join(femur, by = "SEQN") %>%
    left_join(spine, by = "SEQN") %>%
    left_join(rhq, by = "SEQN")

  nhanes_join_qc[[spec$cycle]] <<- tibble(
    cohort = "NHANES",
    cycle = spec$cycle,
    analysis_role = spec$analysis_role,
    demo_rows = nrow(demo),
    demo_unique_ids = n_distinct(demo$SEQN),
    bmx_match_pct = mean(demo$SEQN %in% bmx$SEQN) * 100,
    femur_match_pct = mean(demo$SEQN %in% femur$SEQN) * 100,
    spine_match_pct = mean(demo$SEQN %in% spine$SEQN) * 100,
    reproductive_match_pct = mean(demo$SEQN %in% rhq$SEQN) * 100,
    rows_after_join = nrow(joined)
  )

  joined %>%
    transmute(
      cohort = "NHANES",
      cycle = spec$cycle,
      analysis_role = spec$analysis_role,
      id = clean_id(SEQN),
      age = clean_range(RIDAGEYR, 20, 85),
      female = clean_binary(RIAGENDR, zero_values = 1, one_values = 2),
      weight_kg = clean_range(BMXWT, 25, 250),
      bmi = clean_range(BMXBMI, 10, 80),
      survey_weight_raw = clean_positive_weight(
        get_column(joined, spec$weight_variable)
      ),
      survey_weight_pooled = clean_positive_weight(
        get_column(joined, spec$weight_variable)
      ) / spec$weight_divisor,
      survey_weight_scope = case_when(
        spec$analysis_role == "development" ~ "NHANES 2005-2010 pooled",
        TRUE ~ spec$cycle
      ),
      strata = if_else(!is.na(SDMVSTRA), paste0(spec$cycle, "_", SDMVSTRA),
                       NA_character_),
      psu = if_else(!is.na(SDMVPSU), paste0(spec$cycle, "_", SDMVPSU),
                    NA_character_),
      fn_bmd = clean_range(get_column(joined, "DXXNKBMD"), 0.30, 1.80),
      total_hip_bmd = clean_range(get_column(joined, "DXXOFBMD"), 0.30, 1.80),
      lumbar_spine_bmd = clean_range(get_column(joined, "DXXOSBMD"), 0.30, 2.00),
      femur_invalidity_code = as_num(get_column(joined, "DXXFMBCC")),
      spine_invalidity_code = as_num(get_column(joined, "DXXOSBCC")),
      period_past_12_months = as_num(get_column(joined, "RHQ031")),
      no_period_reason = as_num(get_column(joined, "RHD042")),
      age_last_period = clean_range(get_column(joined, "RHQ060"), 8, 75),
      hysterectomy = clean_binary(get_column(joined, "RHD280"),
                                  zero_values = 2, one_values = 1),
      bilateral_oophorectomy = clean_binary(get_column(joined, "RHQ305"),
                                             zero_values = 2, one_values = 1),
      postmenopausal_or_hysterectomy = case_when(
        female == 1L & no_period_reason == 7 ~ 1L,
        female == 1L & period_past_12_months == 1 ~ 0L,
        TRUE ~ NA_integer_
      ),
      natural_postmenopause_strict = case_when(
        female == 1L & no_period_reason == 7 & hysterectomy == 0L ~ 1L,
        female == 1L & period_past_12_months == 1 ~ 0L,
        TRUE ~ NA_integer_
      )
    ) %>%
    mutate(
      fn_t_us = (fn_bmd - 0.858) / 0.120,
      total_hip_t_us = (total_hip_bmd - 0.942) / 0.122,
      lumbar_spine_t_us = (lumbar_spine_bmd - 1.040) / 0.110,
      osteoporosis_fn_us = if_else(is.na(fn_t_us), NA_integer_,
                                   as.integer(fn_t_us <= -2.5)),
      osteoporosis_total_hip_us = if_else(is.na(total_hip_t_us), NA_integer_,
                                          as.integer(total_hip_t_us <= -2.5)),
      osteoporosis_lumbar_us = if_else(is.na(lumbar_spine_t_us), NA_integer_,
                                       as.integer(lumbar_spine_t_us <= -2.5)),
      osteoporosis_any_us = any_tscore_osteoporosis(
        fn_t_us, total_hip_t_us, lumbar_spine_t_us
      ),
      osteoporosis_any_complete_us = complete_tscore_osteoporosis(
        fn_t_us, total_hip_t_us, lumbar_spine_t_us
      ),
      complete_three_site_bmd = complete.cases(
        fn_bmd, total_hip_bmd, lumbar_spine_bmd
      ),
      osta_raw = if_else(!is.na(age) & !is.na(weight_kg),
                         0.2 * (weight_kg - age), NA_real_),
      osta = trunc(osta_raw),
      eligible_age = age >= age_min & age <= age_max,
      analysis_eligible = female == 1L & eligible_age & !is.na(weight_kg) &
        !is.na(fn_bmd) & !is.na(survey_weight_pooled) &
        !is.na(strata) & !is.na(psu)
    )
}

nhanes_clean <- bind_rows(lapply(seq_len(nrow(nh_periods)), function(i) {
  clean_nhanes_period(nh_periods[i, ])
}))

nhanes_analysis <- nhanes_clean %>% filter(analysis_eligible)

make_nhanes_flow <- function(data, label) {
  bind_rows(
    flow_row(label, "01 merged demographic records", data),
    flow_row(label, "02 women", data, data$female == 1L),
    flow_row(label, "03 women age 50-79", data,
             data$female == 1L & data$eligible_age),
    flow_row(label, "04 valid measured weight", data,
             data$female == 1L & data$eligible_age & !is.na(data$weight_kg)),
    flow_row(label, "05 valid femoral-neck BMD", data,
             data$female == 1L & data$eligible_age & !is.na(data$weight_kg) &
               !is.na(data$fn_bmd)),
    flow_row(label, "06 positive survey weight and design", data,
             data$analysis_eligible)
  )
}

nhanes_flow <- bind_rows(
  make_nhanes_flow(
    nhanes_clean %>% filter(analysis_role == "development"),
    "NHANES development (2005-2010)"
  ),
  make_nhanes_flow(
    nhanes_clean %>% filter(cycle == "2013-2014"),
    "NHANES temporal validation (2013-2014)"
  ),
  make_nhanes_flow(
    nhanes_clean %>% filter(cycle == "2017-March 2020"),
    "NHANES temporal validation (2017-March 2020)"
  )
)

# -----------------------------------------------------------------------------
# 4. KNHANES ingestion and cleaning
# -----------------------------------------------------------------------------

knhanes_input_files <- character()
knhanes_join_qc <- list()

clean_knhanes_year <- function(year) {
  yy <- substr(as.character(year), 3, 4)
  all_file <- file.path(kn_root, as.character(year), paste0("hn", yy, "_all.sas7bdat"))
  dxa_file <- file.path(kn_root, as.character(year), paste0("hn", yy, "_dxa.sas7bdat"))
  assert_file(all_file)
  all_data <- read_sas(all_file)
  all_id <- if ("ID" %in% names(all_data)) "ID" else "id"

  if (year == 2024L) {
    # In 2024 the DXA fields are embedded in the all-participant file.
    # The official DXA target population was age 40 years and older.
    dxa <- all_data %>% filter(as_num(age) >= 40)
    merged <- dxa
    dxa_file <- NA_character_
    knhanes_input_files <<- c(knhanes_input_files, all_file)
    weight_name <- "wt_itvex"
    weight_divisor <- 1
    analysis_role <- "contemporary_temporal_validation"
  } else {
    assert_file(dxa_file)
    knhanes_input_files <<- c(knhanes_input_files, all_file, dxa_file)
    dxa <- read_sas(dxa_file)
    by_map <- setNames(all_id, "ID")
    merged <- left_join(dxa, all_data, by = by_map, suffix = c(".dxa", ""))
    weight_name <- dplyr::case_when(
      year == 2008 && "wt_ex1.dxa" %in% names(merged) ~ "wt_ex1.dxa",
      year == 2009 && "wt_dw.dxa" %in% names(merged) ~ "wt_dw.dxa",
      year == 2010 && "wt_itvex.dxa" %in% names(merged) ~ "wt_itvex.dxa",
      year == 2011 && "wt_ex1.dxa" %in% names(merged) ~ "wt_ex1.dxa",
      TRUE ~ NA_character_
    )
    weight_divisor <- 4
    analysis_role <- "historical_geographic_validation"
  }
  if (is.na(weight_name)) stop("No approved DXA survey weight found for KNHANES ", year)

  knhanes_join_qc[[as.character(year)]] <<- tibble(
    cohort = "KNHANES",
    cycle = as.character(year),
    analysis_role = analysis_role,
    dxa_rows = nrow(dxa),
    dxa_unique_ids = n_distinct(dxa[[if ("ID" %in% names(dxa)) "ID" else "id"]]),
    all_rows = nrow(all_data),
    all_unique_ids = n_distinct(all_data[[all_id]]),
    all_data_match_pct = mean(
      dxa[[if ("ID" %in% names(dxa)) "ID" else "id"]] %in% all_data[[all_id]]
    ) * 100,
    rows_after_join = nrow(merged),
    selected_weight = weight_name
  )

  tibble(
    cohort = "KNHANES",
    cycle = as.character(year),
    analysis_role = analysis_role,
    id = clean_id(get_column(merged, c("ID", "id"), default = NA_character_)),
    age = clean_range(get_column(merged, c("age", "age.dxa")), 20, 90),
    female = clean_binary(get_column(merged, c("sex", "sex.dxa")),
                          zero_values = 1, one_values = 2),
    weight_kg = clean_range(get_column(merged, c("HE_wt", "HE_wt.dxa")), 25, 250),
    bmi = clean_range(get_column(merged, c("HE_BMI", "HE_BMI.dxa")), 10, 80),
    survey_weight_raw = clean_positive_weight(merged[[weight_name]]),
    survey_weight_pooled = clean_positive_weight(merged[[weight_name]]) / weight_divisor,
    survey_weight_scope = if_else(
      year <= 2011,
      "KNHANES 2008-2011 pooled",
      "KNHANES 2024"
    ),
    strata = {
      z <- get_column(merged, c("kstrata", "kstrata.dxa"))
      ifelse(is.na(z), NA_character_, paste0(year, "_", z))
    },
    psu = {
      z <- get_column(merged, c("psu", "psu.dxa"))
      ifelse(is.na(z), NA_character_, paste0(year, "_", z))
    },
    fn_bmd = clean_range(get_column(merged, c("DX_NK_BMD", "DX_NK_BMD.dxa")),
                         0.30, 1.80),
    total_hip_bmd = clean_range(get_column(merged, c("DX_TFM_BMD", "DX_TFM_BMD.dxa")),
                                0.30, 1.80),
    lumbar_spine_bmd = clean_range(get_column(merged, c("DX_LSP_BMD", "DX_LSP_BMD.dxa")),
                                   0.30, 2.00),
    fn_t_asian = clean_range(get_column(merged, c("DX_FN_Ts_A", "DX_FN_Ts_A.dxa")),
                             -10, 10),
    total_hip_t_asian = clean_range(get_column(merged, c("DX_F_Ts_A", "DX_F_Ts_A.dxa")),
                                    -10, 10),
    lumbar_spine_t_asian = clean_range(get_column(merged, c("DX_S_Ts_A", "DX_S_Ts_A.dxa")),
                                       -10, 10),
    menopause_code = as_num(get_column(merged, c("DX_Q_MP", "DX_Q_MP.dxa"))),
    natural_postmenopause = case_when(
      menopause_code == 1 ~ 1L,
      menopause_code %in% c(2, 3) ~ 0L,
      TRUE ~ NA_integer_
    ),
    artificial_postmenopause = case_when(
      menopause_code == 3 ~ 1L,
      menopause_code %in% c(1, 2) ~ 0L,
      TRUE ~ NA_integer_
    ),
    postmenopausal = case_when(
      menopause_code %in% c(1, 3) ~ 1L,
      menopause_code == 2 ~ 0L,
      TRUE ~ NA_integer_
    ),
    grip_right_max = pmax(
      as_num(get_column(merged, "GS_mea_r_1")),
      as_num(get_column(merged, "GS_mea_r_2")), na.rm = TRUE
    ),
    grip_left_max = pmax(
      as_num(get_column(merged, "GS_mea_l_1")),
      as_num(get_column(merged, "GS_mea_l_2")), na.rm = TRUE
    ),
    sarcopenia_dxa = as_num(get_column(merged, "GS_SP_DXA"))
  ) %>%
    mutate(
      grip_right_max = if_else(is.infinite(grip_right_max), NA_real_, grip_right_max),
      grip_left_max = if_else(is.infinite(grip_left_max), NA_real_, grip_left_max),
      grip_max = pmax(grip_right_max, grip_left_max, na.rm = TRUE),
      grip_max = if_else(is.infinite(grip_max), NA_real_, grip_max),
      grip_imbalance_pct = if_else(
        !is.na(grip_right_max) & !is.na(grip_left_max) & grip_max > 0,
        abs(grip_right_max - grip_left_max) / grip_max * 100,
        NA_real_
      ),
      fn_t_us = (fn_bmd - 0.858) / 0.120,
      total_hip_t_us = (total_hip_bmd - 0.942) / 0.122,
      lumbar_spine_t_us = (lumbar_spine_bmd - 1.040) / 0.110,
      osteoporosis_fn_asian = if_else(is.na(fn_t_asian), NA_integer_,
                                      as.integer(fn_t_asian <= -2.5)),
      osteoporosis_total_hip_asian = if_else(
        is.na(total_hip_t_asian), NA_integer_,
        as.integer(total_hip_t_asian <= -2.5)
      ),
      osteoporosis_lumbar_asian = if_else(
        is.na(lumbar_spine_t_asian), NA_integer_,
        as.integer(lumbar_spine_t_asian <= -2.5)
      ),
      osteoporosis_any_asian = any_tscore_osteoporosis(
        fn_t_asian, total_hip_t_asian, lumbar_spine_t_asian
      ),
      osteoporosis_any_complete_asian = complete_tscore_osteoporosis(
        fn_t_asian, total_hip_t_asian, lumbar_spine_t_asian
      ),
      osteoporosis_fn_us = if_else(is.na(fn_t_us), NA_integer_,
                                   as.integer(fn_t_us <= -2.5)),
      osteoporosis_total_hip_us = if_else(is.na(total_hip_t_us), NA_integer_,
                                          as.integer(total_hip_t_us <= -2.5)),
      osteoporosis_lumbar_us = if_else(is.na(lumbar_spine_t_us), NA_integer_,
                                       as.integer(lumbar_spine_t_us <= -2.5)),
      osteoporosis_any_us = any_tscore_osteoporosis(
        fn_t_us, total_hip_t_us, lumbar_spine_t_us
      ),
      osteoporosis_any_complete_us = complete_tscore_osteoporosis(
        fn_t_us, total_hip_t_us, lumbar_spine_t_us
      ),
      complete_three_site_bmd = complete.cases(
        fn_bmd, total_hip_bmd, lumbar_spine_bmd
      ),
      osta_raw = if_else(!is.na(age) & !is.na(weight_kg),
                         0.2 * (weight_kg - age), NA_real_),
      osta = trunc(osta_raw),
      eligible_age = age >= age_min & age <= age_max,
      analysis_eligible = female == 1L & eligible_age & !is.na(weight_kg) &
        !is.na(fn_bmd) & !is.na(fn_t_asian) &
        !is.na(survey_weight_pooled) & !is.na(strata) & !is.na(psu)
    )
}

knhanes_clean <- bind_rows(lapply(kn_years, clean_knhanes_year))
knhanes_analysis <- knhanes_clean %>% filter(analysis_eligible)

make_knhanes_flow <- function(data, label) {
  bind_rows(
    flow_row(label, "01 DXA records", data),
    flow_row(label, "02 women", data, data$female == 1L),
    flow_row(label, "03 women age 50-79", data,
             data$female == 1L & data$eligible_age),
    flow_row(label, "04 valid measured weight", data,
             data$female == 1L & data$eligible_age & !is.na(data$weight_kg)),
    flow_row(label, "05 valid femoral-neck BMD and Asian T-score", data,
             data$female == 1L & data$eligible_age & !is.na(data$weight_kg) &
               !is.na(data$fn_bmd) & !is.na(data$fn_t_asian)),
    flow_row(label, "06 positive survey weight and design", data,
             data$analysis_eligible)
  )
}

knhanes_flow <- bind_rows(
  make_knhanes_flow(
    knhanes_clean %>% filter(cycle != "2024"),
    "KNHANES historical validation (2008-2011)"
  ),
  make_knhanes_flow(
    knhanes_clean %>% filter(cycle == "2024"),
    "KNHANES contemporary validation (2024)"
  )
)

# -----------------------------------------------------------------------------
# 5. CHARLS ingestion, official weights, death, and longitudinal outcomes
# -----------------------------------------------------------------------------

assert_file(charls_long_file)
charls_input_files <- charls_long_file

harmonized_file <- file.path(
  charls_original_root, "Harmonized_CHARLS_D", "H_CHARLS_D_Data.dta"
)
assert_file(harmonized_file)
charls_input_files <- c(charls_input_files, harmonized_file)

harmonized <- read_dta(harmonized_file)
id_crosswalk <- harmonized %>%
  transmute(id = clean_id(ID), id_wave1 = clean_id(ID_w1)) %>%
  filter(!is.na(id), !is.na(id_wave1)) %>%
  distinct(id_wave1, .keep_all = TRUE)

weight_specs <- tibble(
  wave = 1:5,
  path = c(
    file.path(charls_raw_root, "2011", "weight.dta"),
    file.path(charls_raw_root, "2013", "Weights.dta"),
    file.path(charls_raw_root, "2015", "Weights.dta"),
    file.path(charls_original_root, "2018", "CHARLS2018r", "Weights.dta"),
    file.path(charls_original_root, "2020", "CHARLS2020r", "Weights.dta")
  )
)
invisible(lapply(weight_specs$path, assert_file))
charls_input_files <- c(charls_input_files, weight_specs$path)

charls_weights <- bind_rows(lapply(seq_len(nrow(weight_specs)), function(i) {
  wave_i <- weight_specs$wave[i]
  path_i <- weight_specs$path[i]
  x <- read_dta(path_i)
  if (wave_i == 1) {
    x <- x %>%
      transmute(
        id_wave1 = clean_id(ID),
        individual_weight = clean_positive_weight(ind_weight_ad2)
      ) %>%
      left_join(id_crosswalk, by = "id_wave1") %>%
      select(id, individual_weight)
  } else {
    x <- x %>%
      transmute(
        id = clean_id(ID),
        individual_weight = clean_positive_weight(INDV_weight_ad2)
      )
  }
  x %>% mutate(wave = wave_i, .before = 1)
})) %>%
  filter(!is.na(id)) %>%
  distinct(wave, id, .keep_all = TRUE)

death_harmonized <- harmonized %>%
  transmute(
    id = clean_id(ID),
    death_year_harmonized = clean_range(radyear, 2011, 2021),
    death_month_harmonized = clean_range(radmonth, 1, 12),
    in_wave1_h = as_num(inw1),
    in_wave2_h = as_num(inw2),
    in_wave3_h = as_num(inw3),
    in_wave4_h = as_num(inw4)
  ) %>%
  filter(!is.na(id)) %>%
  distinct(id, .keep_all = TRUE)

exit_2013_file <- file.path(charls_raw_root, "2013", "Exit_Interview.dta")
exit_2020_file <- file.path(
  charls_original_root, "2020", "CHARLS2020r", "Exit_Module.dta"
)
assert_file(exit_2013_file)
assert_file(exit_2020_file)
charls_input_files <- c(charls_input_files, exit_2013_file, exit_2020_file)

exit_2013 <- read_dta(exit_2013_file) %>%
  transmute(
    id = clean_id(ID),
    source = "2013 exit interview",
    death_year = clean_range(exb001_1, 2011, 2021),
    death_month = clean_range(exb001_2, 1, 12),
    proxy_hip_fracture = clean_binary(exda025, zero_values = 2, one_values = 1)
  )

exit_2020 <- read_dta(exit_2020_file) %>%
  transmute(
    id = clean_id(ID),
    source = "2020 exit module",
    death_year = clean_range(exb001_1, 2011, 2021),
    death_month = clean_range(exb001_2, 1, 12),
    proxy_hip_fracture = clean_binary(exda015, zero_values = 2, one_values = 1)
  )

death_exit <- bind_rows(exit_2013, exit_2020) %>%
  filter(!is.na(id)) %>%
  group_by(id) %>%
  summarise(
    death_year_exit = min_or_na(death_year),
    death_month_exit = {
      if (all(is.na(death_year))) NA_real_ else {
        target <- min(death_year, na.rm = TRUE)
        min_or_na(death_month[death_year == target])
      }
    },
    death_proxy_hip_fracture = any_or_na(proxy_hip_fracture),
    .groups = "drop"
  )

charls_death <- full_join(death_harmonized, death_exit, by = "id") %>%
  rowwise() %>%
  mutate(
    death_year = min_or_na(c(death_year_harmonized, death_year_exit)),
    death_month = case_when(
      !is.na(death_year_exit) & death_year == death_year_exit ~ death_month_exit,
      !is.na(death_year_harmonized) & death_year == death_year_harmonized ~
        death_month_harmonized,
      TRUE ~ NA_real_
    ),
    death_observed = as.integer(!is.na(death_year))
  ) %>%
  ungroup()

charls_raw <- read.csv(charls_long_file, check.names = FALSE)
names(charls_raw) <- make.names(names(charls_raw), unique = TRUE)

required_charls_columns <- c(
  "ID", "wave", "gender", "age", "mweight", "mheight", "bmi", "hip",
  "fall_down", "communityID", "province", "rural", "stroke", "arthre",
  "smoken", "drinkev", "exercise"
)
missing_charls_columns <- setdiff(required_charls_columns, names(charls_raw))
if (length(missing_charls_columns) > 0) {
  stop("CHARLS.csv is missing required fields: ",
       paste(missing_charls_columns, collapse = ", "))
}

charls_clean <- charls_raw %>%
  transmute(
    cohort = "CHARLS",
    id = clean_id(ID),
    wave = as.integer(wave),
    survey_year = as.integer(unname(wave_year[as.character(wave)])),
    household_id = clean_id(householdID),
    community_id = clean_id(communityID),
    province = as.character(province),
    rural = clean_binary(rural, zero_values = 0, one_values = 1),
    age = clean_range(age, 40, 120),
    female = clean_binary(gender, zero_values = 1, one_values = 0),
    weight_kg = clean_range(mweight, 25, 250),
    height_cm = clean_range(mheight, 100, 220),
    bmi_reported = clean_range(bmi, 10, 80),
    hip_fracture = clean_binary(hip, zero_values = 0, one_values = 1),
    fall = clean_binary(fall_down, zero_values = 0, one_values = 1),
    stroke = clean_binary(stroke, zero_values = 0, one_values = 1),
    arthritis = clean_binary(arthre, zero_values = 0, one_values = 1),
    current_smoking = clean_binary(smoken, zero_values = 0, one_values = 1),
    ever_drinking = clean_binary(drinkev, zero_values = 0, one_values = 1),
    exercise = clean_binary(exercise, zero_values = 0, one_values = 1)
  ) %>%
  mutate(
    bmi_calculated = if_else(
      !is.na(weight_kg) & !is.na(height_cm),
      weight_kg / (height_cm / 100)^2,
      NA_real_
    ),
    bmi_calculated = clean_range(bmi_calculated, 10, 80),
    bmi = coalesce(bmi_reported, bmi_calculated),
    osta_raw = if_else(!is.na(age) & !is.na(weight_kg),
                       0.2 * (weight_kg - age), NA_real_),
    osta = trunc(osta_raw)
  ) %>%
  left_join(charls_weights, by = c("wave", "id")) %>%
  left_join(charls_death, by = "id")

charls_duplicate_qc <- duplicate_key_report(
  charls_clean, c("id", "wave"), "CHARLS clean long"
)
if (charls_duplicate_qc$duplicate_key_rows > 0) {
  stop("Duplicate CHARLS id-wave rows detected; cleaning stopped to prevent join inflation.")
}

charls_followup <- charls_clean %>%
  filter(wave > 1) %>%
  group_by(id) %>%
  summarise(
    hip_followup_observed_n = sum(!is.na(hip_fracture)),
    incident_hip = any_or_na(hip_fracture),
    incident_hip_year = min_or_na(survey_year[hip_fracture == 1]),
    fall_followup_observed_n = sum(!is.na(fall)),
    future_fall = any_or_na(fall),
    last_hip_observed_year = max_or_na(survey_year[!is.na(hip_fracture)]),
    last_any_observed_year = max_or_na(survey_year),
    .groups = "drop"
  )

charls_baseline <- charls_clean %>%
  filter(wave == 1) %>%
  select(
    cohort, id, baseline_year = survey_year, age, female, weight_kg,
    height_cm, bmi, osta_raw, osta, baseline_hip = hip_fracture,
    baseline_fall = fall, individual_weight, community_id, province, rural,
    stroke, arthritis, current_smoking, ever_drinking, exercise,
    death_year, death_month, death_observed, death_proxy_hip_fracture
  ) %>%
  left_join(charls_followup, by = "id") %>%
  mutate(
    eligible_age = age >= age_min & age <= age_max,
    incident_hip_augmented = case_when(
      incident_hip == 1L | death_proxy_hip_fracture == 1L ~ 1L,
      incident_hip == 0L & (is.na(death_proxy_hip_fracture) |
                             death_proxy_hip_fracture == 0L) ~ 0L,
      TRUE ~ NA_integer_
    ),
    competing_status = case_when(
      !is.na(incident_hip_year) &
        (is.na(death_year) | incident_hip_year <= death_year) ~ 1L,
      !is.na(death_year) &
        (is.na(incident_hip_year) | death_year < incident_hip_year) ~ 2L,
      TRUE ~ 0L
    ),
    followup_end_year = case_when(
      competing_status == 1L ~ incident_hip_year,
      competing_status == 2L ~ death_year,
      TRUE ~ last_hip_observed_year
    ),
    followup_years = followup_end_year - baseline_year,
    analysis_eligible = female == 1L & eligible_age & baseline_hip == 0L &
      !is.na(weight_kg) & hip_followup_observed_n > 0 & followup_years > 0,
    competing_analysis_eligible = female == 1L & eligible_age &
      baseline_hip == 0L & !is.na(weight_kg) &
      (hip_followup_observed_n > 0 | death_observed == 1L) & followup_years > 0
  )

charls_analysis <- charls_baseline %>% filter(analysis_eligible)
charls_competing_analysis <- charls_baseline %>% filter(competing_analysis_eligible)

charls_flow <- bind_rows(
  flow_row("CHARLS", "01 baseline wave records", charls_baseline),
  flow_row("CHARLS", "02 women", charls_baseline, charls_baseline$female == 1L),
  flow_row("CHARLS", "03 women age 50-79", charls_baseline,
           charls_baseline$female == 1L & charls_baseline$eligible_age),
  flow_row("CHARLS", "04 no baseline hip fracture", charls_baseline,
           charls_baseline$female == 1L & charls_baseline$eligible_age &
             charls_baseline$baseline_hip == 0L),
  flow_row("CHARLS", "05 valid measured weight", charls_baseline,
           charls_baseline$female == 1L & charls_baseline$eligible_age &
             charls_baseline$baseline_hip == 0L & !is.na(charls_baseline$weight_kg)),
  flow_row("CHARLS", "06 at least one observed follow-up hip status", charls_baseline,
           charls_baseline$analysis_eligible)
)

# Validate the already-audited consolidated hip fields against official counts.
official_hip_counts_file <- file.path(
  root, "charls_official", "audit_original", "key_variable_validation.csv"
)
assert_file(official_hip_counts_file)
charls_input_files <- c(charls_input_files, official_hip_counts_file)
official_hip_counts <- read.csv(official_hip_counts_file)
current_hip_counts <- charls_clean %>%
  group_by(wave) %>%
  summarise(
    cleaned_rows = n(),
    cleaned_observed_hip = sum(!is.na(hip_fracture)),
    cleaned_hip_yes_recomputed = sum(hip_fracture == 1, na.rm = TRUE),
    .groups = "drop"
  )
charls_hip_audit <- official_hip_counts %>%
  transmute(
    survey_year = as.integer(wave),
    wave = match(survey_year, c(2011L, 2013L, 2015L, 2018L, 2020L)),
    raw_hip_variables,
    official_raw_hip_yes = raw_hip_yes,
    prior_cleaned_hip_yes = cleaned_hip_yes,
    official_positive_id_set_match = positive_id_set_match
  ) %>%
  left_join(current_hip_counts, by = "wave") %>%
  mutate(count_match_now = official_raw_hip_yes == cleaned_hip_yes_recomputed)

if (any(is.na(charls_hip_audit$count_match_now)) ||
    !all(charls_hip_audit$count_match_now)) {
  stop("CHARLS hip-fracture count audit failed in at least one wave.")
}

# -----------------------------------------------------------------------------
# 6. Harmonized analysis-ready cross-sectional DXA table
# -----------------------------------------------------------------------------

imaging_analysis <- bind_rows(
  nhanes_analysis %>%
    transmute(
      cohort, cycle, analysis_role, id, age, female, weight_kg, bmi,
      postmenopausal = postmenopausal_or_hysterectomy,
      survey_weight = survey_weight_pooled, survey_weight_scope, strata, psu,
      fn_bmd, total_hip_bmd, lumbar_spine_bmd,
      fn_t_local = fn_t_us,
      total_hip_t_local = total_hip_t_us,
      lumbar_spine_t_local = lumbar_spine_t_us,
      fn_t_us, total_hip_t_us, lumbar_spine_t_us,
      osteoporosis_fn_local = osteoporosis_fn_us,
      osteoporosis_total_hip_local = osteoporosis_total_hip_us,
      osteoporosis_lumbar_local = osteoporosis_lumbar_us,
      osteoporosis_any_local = osteoporosis_any_us,
      osteoporosis_any_complete_local = osteoporosis_any_complete_us,
      osteoporosis_fn_us, osteoporosis_total_hip_us,
      osteoporosis_lumbar_us, osteoporosis_any_us,
      osteoporosis_any_complete_us, complete_three_site_bmd,
      osta_raw, osta,
      grip_right_max = NA_real_, grip_left_max = NA_real_,
      grip_max = NA_real_, grip_imbalance_pct = NA_real_,
      sarcopenia_dxa = NA_real_
    ),
  knhanes_analysis %>%
    transmute(
      cohort, cycle, analysis_role, id, age, female, weight_kg, bmi, postmenopausal,
      survey_weight = survey_weight_pooled, survey_weight_scope, strata, psu,
      fn_bmd, total_hip_bmd, lumbar_spine_bmd,
      fn_t_local = fn_t_asian,
      total_hip_t_local = total_hip_t_asian,
      lumbar_spine_t_local = lumbar_spine_t_asian,
      fn_t_us, total_hip_t_us, lumbar_spine_t_us,
      osteoporosis_fn_local = osteoporosis_fn_asian,
      osteoporosis_total_hip_local = osteoporosis_total_hip_asian,
      osteoporosis_lumbar_local = osteoporosis_lumbar_asian,
      osteoporosis_any_local = osteoporosis_any_asian,
      osteoporosis_any_complete_local = osteoporosis_any_complete_asian,
      osteoporosis_fn_us, osteoporosis_total_hip_us,
      osteoporosis_lumbar_us, osteoporosis_any_us,
      osteoporosis_any_complete_us, complete_three_site_bmd,
      osta_raw, osta, grip_right_max, grip_left_max, grip_max,
      grip_imbalance_pct, sarcopenia_dxa
    )
)

# -----------------------------------------------------------------------------
# 7. Quality-control outputs and frozen clean data
# -----------------------------------------------------------------------------

all_flow <- bind_rows(nhanes_flow, knhanes_flow, charls_flow)
join_qc <- bind_rows(bind_rows(nhanes_join_qc), bind_rows(knhanes_join_qc))

duplicate_qc <- bind_rows(
  duplicate_key_report(nhanes_clean, c("cycle", "id"), "NHANES clean"),
  duplicate_key_report(knhanes_clean, c("cycle", "id"), "KNHANES clean"),
  charls_duplicate_qc,
  duplicate_key_report(charls_baseline, "id", "CHARLS baseline")
)

missing_qc <- bind_rows(
  missingness_report(nhanes_analysis, "NHANES analysis-ready"),
  missingness_report(knhanes_analysis, "KNHANES analysis-ready"),
  missingness_report(charls_analysis, "CHARLS analysis-ready")
)

outcome_qc <- bind_rows(
  nhanes_analysis %>%
    group_by(cohort, cycle) %>%
    summarise(
      n = n(),
      events_local_fn = sum(osteoporosis_fn_us == 1, na.rm = TRUE),
      events_local_total_hip = sum(osteoporosis_total_hip_us == 1, na.rm = TRUE),
      events_local_lumbar = sum(osteoporosis_lumbar_us == 1, na.rm = TRUE),
      events_local_any = sum(osteoporosis_any_us == 1, na.rm = TRUE),
      events_local_any_complete = sum(osteoporosis_any_complete_us == 1, na.rm = TRUE),
      events_us_fn = sum(osteoporosis_fn_us == 1, na.rm = TRUE),
      events_us_total_hip = sum(osteoporosis_total_hip_us == 1, na.rm = TRUE),
      events_us_lumbar = sum(osteoporosis_lumbar_us == 1, na.rm = TRUE),
      events_us_any = sum(osteoporosis_any_us == 1, na.rm = TRUE),
      events_us_any_complete = sum(osteoporosis_any_complete_us == 1, na.rm = TRUE),
      complete_three_site_n = sum(complete_three_site_bmd),
      mean_age = mean(age, na.rm = TRUE),
      mean_weight = mean(weight_kg, na.rm = TRUE),
      .groups = "drop"
    ),
  knhanes_analysis %>%
    group_by(cohort, cycle) %>%
    summarise(
      n = n(),
      events_local_fn = sum(osteoporosis_fn_asian == 1, na.rm = TRUE),
      events_local_total_hip = sum(osteoporosis_total_hip_asian == 1, na.rm = TRUE),
      events_local_lumbar = sum(osteoporosis_lumbar_asian == 1, na.rm = TRUE),
      events_local_any = sum(osteoporosis_any_asian == 1, na.rm = TRUE),
      events_local_any_complete = sum(osteoporosis_any_complete_asian == 1, na.rm = TRUE),
      events_us_fn = sum(osteoporosis_fn_us == 1, na.rm = TRUE),
      events_us_total_hip = sum(osteoporosis_total_hip_us == 1, na.rm = TRUE),
      events_us_lumbar = sum(osteoporosis_lumbar_us == 1, na.rm = TRUE),
      events_us_any = sum(osteoporosis_any_us == 1, na.rm = TRUE),
      events_us_any_complete = sum(osteoporosis_any_complete_us == 1, na.rm = TRUE),
      complete_three_site_n = sum(complete_three_site_bmd),
      mean_age = mean(age, na.rm = TRUE),
      mean_weight = mean(weight_kg, na.rm = TRUE),
      .groups = "drop"
    )
)

charls_outcome_qc <- charls_analysis %>%
  summarise(
    cohort = "CHARLS",
    n = n(),
    incident_hip_events = sum(incident_hip == 1, na.rm = TRUE),
    augmented_hip_events = sum(incident_hip_augmented == 1, na.rm = TRUE),
    competing_deaths = sum(competing_status == 2, na.rm = TRUE),
    future_falls = sum(future_fall == 1, na.rm = TRUE),
    median_followup_years = median(followup_years, na.rm = TRUE),
    baseline_weight_available_pct = mean(!is.na(individual_weight)) * 100
  )

postclean_validation <- bind_rows(
  tibble(
    dataset = "NHANES analysis-ready",
    rule = c(
      "age within 50-79", "weight within 25-250 kg",
      "femoral-neck BMD within 0.30-1.80", "binary osteoporosis outcome",
      "positive pooled survey weight", "complete survey design"
    ),
    failures = c(
      sum(is.na(nhanes_analysis$age) | nhanes_analysis$age < age_min |
            nhanes_analysis$age > age_max),
      sum(is.na(nhanes_analysis$weight_kg) | nhanes_analysis$weight_kg < 25 |
            nhanes_analysis$weight_kg > 250),
      sum(is.na(nhanes_analysis$fn_bmd) | nhanes_analysis$fn_bmd < 0.30 |
            nhanes_analysis$fn_bmd > 1.80),
      sum(!nhanes_analysis$osteoporosis_fn_us %in% c(0L, 1L)),
      sum(is.na(nhanes_analysis$survey_weight_pooled) |
            nhanes_analysis$survey_weight_pooled <= 0),
      sum(is.na(nhanes_analysis$strata) | is.na(nhanes_analysis$psu))
    )
  ),
  tibble(
    dataset = "KNHANES analysis-ready",
    rule = c(
      "age within 50-79", "weight within 25-250 kg",
      "femoral-neck BMD within 0.30-1.80", "binary local osteoporosis outcome",
      "positive pooled survey weight", "complete survey design"
    ),
    failures = c(
      sum(is.na(knhanes_analysis$age) | knhanes_analysis$age < age_min |
            knhanes_analysis$age > age_max),
      sum(is.na(knhanes_analysis$weight_kg) | knhanes_analysis$weight_kg < 25 |
            knhanes_analysis$weight_kg > 250),
      sum(is.na(knhanes_analysis$fn_bmd) | knhanes_analysis$fn_bmd < 0.30 |
            knhanes_analysis$fn_bmd > 1.80),
      sum(!knhanes_analysis$osteoporosis_fn_asian %in% c(0L, 1L)),
      sum(is.na(knhanes_analysis$survey_weight_pooled) |
            knhanes_analysis$survey_weight_pooled <= 0),
      sum(is.na(knhanes_analysis$strata) | is.na(knhanes_analysis$psu))
    )
  ),
  tibble(
    dataset = "CHARLS analysis-ready",
    rule = c(
      "unique baseline ID", "age within 50-79", "weight within 25-250 kg",
      "no baseline hip fracture", "at least one follow-up hip observation",
      "positive follow-up time"
    ),
    failures = c(
      sum(duplicated(charls_analysis$id)),
      sum(is.na(charls_analysis$age) | charls_analysis$age < age_min |
            charls_analysis$age > age_max),
      sum(is.na(charls_analysis$weight_kg) | charls_analysis$weight_kg < 25 |
            charls_analysis$weight_kg > 250),
      sum(is.na(charls_analysis$baseline_hip) | charls_analysis$baseline_hip != 0L),
      sum(is.na(charls_analysis$hip_followup_observed_n) |
            charls_analysis$hip_followup_observed_n <= 0),
      sum(is.na(charls_analysis$followup_years) | charls_analysis$followup_years <= 0)
    )
  )
)

if (any(postclean_validation$failures > 0)) {
  stop("At least one post-clean validation rule failed. See postclean_validation.")
}

input_manifest <- file_manifest(c(
  nhanes_input_files, knhanes_input_files, charls_input_files
))

write_both(nhanes_clean, "nhanes_clean_all")
write_both(nhanes_analysis, "nhanes_analysis_ready")
write_both(knhanes_clean, "knhanes_clean_all")
write_both(knhanes_analysis, "knhanes_analysis_ready")
write_both(charls_clean, "charls_long_clean")
write_both(charls_baseline, "charls_baseline_all")
write_both(charls_analysis, "charls_analysis_ready")
write_both(charls_competing_analysis, "charls_competing_risk_analysis_ready")
write_both(imaging_analysis, "imaging_cross_country_analysis_ready")

write.csv(input_manifest, file.path(qc_dir, "01_input_file_manifest.csv"),
          row.names = FALSE, na = "")
write.csv(all_flow, file.path(qc_dir, "02_sample_flow.csv"),
          row.names = FALSE, na = "")
write.csv(join_qc, file.path(qc_dir, "03_join_coverage.csv"),
          row.names = FALSE, na = "")
write.csv(duplicate_qc, file.path(qc_dir, "04_duplicate_keys.csv"),
          row.names = FALSE, na = "")
write.csv(missing_qc, file.path(qc_dir, "05_missingness.csv"),
          row.names = FALSE, na = "")
write.csv(outcome_qc, file.path(qc_dir, "06_imaging_outcome_counts.csv"),
          row.names = FALSE, na = "")
write.csv(charls_outcome_qc, file.path(qc_dir, "07_charls_outcome_counts.csv"),
          row.names = FALSE, na = "")
write.csv(charls_hip_audit, file.path(qc_dir, "08_charls_hip_audit.csv"),
          row.names = FALSE, na = "")
write.csv(postclean_validation, file.path(qc_dir, "09_postclean_validation.csv"),
          row.names = FALSE, na = "")
writeLines(capture.output(sessionInfo()), file.path(qc_dir, "10_sessionInfo.txt"))

cat("\nSTEP 1 DATA CLEANING COMPLETED\n")
cat("Clean data: ", normalizePath(clean_dir, winslash = "/"), "\n", sep = "")
cat("Quality control: ", normalizePath(qc_dir, winslash = "/"), "\n", sep = "")
cat("NHANES analysis-ready n: ", nrow(nhanes_analysis), "\n", sep = "")
cat("KNHANES analysis-ready n: ", nrow(knhanes_analysis), "\n", sep = "")
cat("CHARLS analysis-ready n: ", nrow(charls_analysis), "\n", sep = "")
