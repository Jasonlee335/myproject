# STEP 6: prospective clinical validation in CHARLS
# Primary exposure: within-CHARLS continuous OSTA risk ranking.
# Primary outcome: first self-reported hip fracture after the 2011 baseline.

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "survival", "sandwich", "lmtest",
  "cmprsk", "splines", "tibble"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages before running Step 6: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(sandwich)
  library(lmtest)
  library(cmprsk)
  library(splines)
  library(tibble)
})

set.seed(20260805)

project_dir <- Sys.getenv("DXA_PROJECT_DIR", unset = file.path(getwd(), "DXA_transportability_results"))
input_dir <- file.path(project_dir, "data_clean")
output_dir <- file.path(project_dir, "step6_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_replicates <- as.integer(Sys.getenv("STEP6_BOOTSTRAP_REPS", unset = "200"))
if (!is.finite(bootstrap_replicates) || bootstrap_replicates < 20) bootstrap_replicates <- 200L
capacities <- c(0.10, 0.20, 0.30, 0.40)
primary_capacity <- 0.20

weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

weighted_quantile <- function(x, w, probabilities) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) return(rep(NA_real_, length(probabilities)))
  order_index <- order(x)
  x <- x[order_index]
  w <- w[order_index]
  cumulative <- cumsum(w) / sum(w)
  vapply(probabilities, function(p) x[which(cumulative >= p)[1]], numeric(1))
}

capacity_rule <- function(score, w, capacity) {
  keep <- is.finite(score) & is.finite(w) & w > 0
  referral <- rep(NA_real_, length(score))
  if (!any(keep)) return(list(cutoff = NA_real_, tie_fraction = NA_real_, referral = referral))
  valid_score <- score[keep]
  valid_weight <- w[keep]
  cutoff <- weighted_quantile(valid_score, valid_weight, 1 - capacity)
  weight_above <- sum(valid_weight[valid_score > cutoff])
  weight_at <- sum(valid_weight[valid_score == cutoff])
  target_weight <- capacity * sum(valid_weight)
  tie_fraction <- if (weight_at > 0) {
    min(1, max(0, (target_weight - weight_above) / weight_at))
  } else {
    0
  }
  referral[keep] <- ifelse(valid_score > cutoff, 1, ifelse(valid_score == cutoff, tie_fraction, 0))
  list(cutoff = cutoff, tie_fraction = tie_fraction, referral = referral)
}

exact_rank_referral <- function(data, score, capacity) {
  keep <- is.finite(score)
  referral <- rep(NA_integer_, nrow(data))
  referral[keep] <- 0L
  eligible_index <- which(keep)
  ordered <- eligible_index[order(-score[eligible_index], as.character(data$id[eligible_index]))]
  n_refer <- max(1L, min(length(ordered), round(capacity * length(ordered))))
  referral[ordered[seq_len(n_refer)]] <- 1L
  referral
}

enrichment_metrics <- function(event, referral, w = rep(1, length(event))) {
  keep <- event %in% c(0, 1) & is.finite(referral) & is.finite(w) & w > 0
  event <- event[keep]
  referral <- pmin(pmax(referral[keep], 0), 1)
  w <- w[keep]
  tp <- sum(w * referral * event)
  fp <- sum(w * referral * (1 - event))
  fn <- sum(w * (1 - referral) * event)
  tn <- sum(w * (1 - referral) * (1 - event))
  referred <- tp + fp
  nonreferred <- fn + tn
  total_events <- tp + fn
  risk_referred <- ifelse(referred > 0, tp / referred, NA_real_)
  risk_nonreferred <- ifelse(nonreferred > 0, fn / nonreferred, NA_real_)
  c(
    referral_rate = referred / (referred + nonreferred),
    event_capture = ifelse(total_events > 0, tp / total_events, NA_real_),
    event_risk_referred = risk_referred,
    event_risk_nonreferred = risk_nonreferred,
    risk_ratio = ifelse(risk_nonreferred > 0, risk_referred / risk_nonreferred, NA_real_),
    risk_difference = risk_referred - risk_nonreferred,
    number_needed_to_screen = ifelse(tp > 0, referred / tp, NA_real_),
    specificity = ifelse((fp + tn) > 0, tn / (fp + tn), NA_real_)
  )
}

bootstrap_interval <- function(x) {
  valid <- x[is.finite(x)]
  if (length(valid) < 10) return(c(NA_real_, NA_real_))
  unname(quantile(valid, c(0.025, 0.975), type = 6, names = FALSE))
}

extract_cox_term <- function(fit, term, analysis, exposure_label) {
  coefficient_table <- summary(fit)$coefficients
  if (!term %in% rownames(coefficient_table)) stop("Cox term not found: ", term)
  log_estimate <- coefficient_table[term, "coef"]
  standard_error_column <- if ("robust se" %in% colnames(coefficient_table)) "robust se" else "se(coef)"
  standard_error <- coefficient_table[term, standard_error_column]
  tibble(
    analysis = analysis,
    exposure = exposure_label,
    estimate_scale = "Hazard ratio",
    estimate = exp(log_estimate),
    lower = exp(log_estimate - 1.96 * standard_error),
    upper = exp(log_estimate + 1.96 * standard_error),
    p_value = 2 * pnorm(-abs(log_estimate / standard_error)),
    model_n = fit$n,
    model_events = fit$nevent
  )
}

extract_glm_term <- function(fit, term, cluster, analysis, exposure_label, estimate_scale) {
  robust_table <- lmtest::coeftest(fit, vcov. = sandwich::vcovCL(fit, cluster = cluster, type = "HC1"))
  if (!term %in% rownames(robust_table)) stop("GLM term not found: ", term)
  log_estimate <- robust_table[term, "Estimate"]
  standard_error <- robust_table[term, "Std. Error"]
  tibble(
    analysis = analysis,
    exposure = exposure_label,
    estimate_scale = estimate_scale,
    estimate = exp(log_estimate),
    lower = exp(log_estimate - 1.96 * standard_error),
    upper = exp(log_estimate + 1.96 * standard_error),
    p_value = robust_table[term, "Pr(>|z|)"],
    model_n = nobs(fit),
    model_events = sum(model.response(model.frame(fit)) == 1)
  )
}

build_person_period <- function(data, outcome = "incident_hip") {
  interval_definition <- tibble(
    interval_index = 1:4,
    interval_start = c(2011, 2013, 2015, 2018),
    interval_end = c(2013, 2015, 2018, 2020)
  )
  data %>%
    mutate(model_event = .data[[outcome]]) %>%
    tidyr::crossing(interval_definition) %>%
    filter(interval_end <= followup_end_year) %>%
    mutate(
      event_interval = as.integer(model_event == 1L & incident_hip_year == interval_end),
      censor_interval = as.integer(
        model_event == 0L & followup_end_year == interval_end & interval_end < 2020
      )
    ) %>%
    group_by(id) %>%
    arrange(interval_index, .by_group = TRUE) %>%
    mutate(prior_event = lag(cumsum(event_interval), default = 0L)) %>%
    filter(prior_event == 0L) %>%
    ungroup()
}

build_landmarks <- function(long_data) {
  landmark_rows <- list()
  landmark_id <- 1L
  for (baseline_wave in 1:3) {
    baseline <- long_data %>%
      filter(
        wave == baseline_wave,
        female == 1L,
        age >= 50,
        age <= 79,
        hip_fracture == 0L,
        complete.cases(weight_kg, osta_raw, community_id)
      ) %>%
      select(
        id, landmark_year = survey_year, age, weight_kg, bmi, osta_raw, osta,
        community_id, rural, stroke, arthritis, current_smoking, ever_drinking,
        baseline_fall = fall
      )
    future <- long_data %>%
      filter(wave > baseline_wave) %>%
      group_by(id) %>%
      summarise(
        observed_n = sum(!is.na(hip_fracture)),
        event = as.integer(any(hip_fracture == 1L, na.rm = TRUE)),
        event_year = ifelse(
          any(hip_fracture == 1L, na.rm = TRUE),
          min(survey_year[hip_fracture == 1L], na.rm = TRUE),
          NA_real_
        ),
        last_observed_year = ifelse(
          any(!is.na(hip_fracture)),
          max(survey_year[!is.na(hip_fracture)], na.rm = TRUE),
          NA_real_
        ),
        .groups = "drop"
      )
    landmark <- baseline %>%
      inner_join(future, by = "id") %>%
      filter(observed_n > 0) %>%
      mutate(
        followup_end_year = if_else(event == 1L, event_year, last_observed_year),
        followup_years = followup_end_year - landmark_year,
        risk_score = -osta_raw,
        risk_z = as.numeric(scale(risk_score)),
        landmark_wave = baseline_wave
      )
    landmark$high20 <- exact_rank_referral(landmark, landmark$risk_score, 0.20)
    landmark <- landmark %>% filter(followup_years > 0)
    landmark_rows[[landmark_id]] <- landmark
    landmark_id <- landmark_id + 1L
  }
  bind_rows(landmark_rows)
}

baseline_file <- file.path(input_dir, "charls_analysis_ready.rds")
competing_file <- file.path(input_dir, "charls_competing_risk_analysis_ready.rds")
baseline_all_file <- file.path(input_dir, "charls_baseline_all.rds")
long_file <- file.path(input_dir, "charls_long_clean.rds")
for (path in c(baseline_file, competing_file, baseline_all_file, long_file)) {
  if (!file.exists(path)) stop("Required Step 1 file is missing: ", path)
}

baseline <- readRDS(baseline_file) %>%
  mutate(
    risk_score = -osta_raw,
    risk_z = as.numeric(scale(risk_score)),
    clinical_osta = as.integer(osta <= -1)
  )
baseline$high20 <- exact_rank_referral(baseline, baseline$risk_score, primary_capacity)
competing <- readRDS(competing_file) %>%
  mutate(
    risk_score = -osta_raw,
    risk_z = (risk_score - mean(baseline$risk_score)) / sd(baseline$risk_score),
    clinical_osta = as.integer(osta <= -1)
  )
competing$high20 <- exact_rank_referral(competing, competing$risk_score, primary_capacity)
baseline_all <- readRDS(baseline_all_file)
long_data <- readRDS(long_file)

# Capacity-constrained event enrichment with community-cluster bootstrap.
strategy_definitions <- tribble(
  ~strategy, ~capacity,
  "Clinical OSTA <= -1", NA_real_,
  "Target-percentile OSTA 10%", 0.10,
  "Target-percentile OSTA 20%", 0.20,
  "Target-percentile OSTA 30%", 0.30,
  "Target-percentile OSTA 40%", 0.40
)

communities <- unique(baseline$community_id)
cluster_index <- match(baseline$community_id, communities)
bootstrap_cluster_counts <- replicate(
  bootstrap_replicates,
  tabulate(sample(seq_along(communities), length(communities), replace = TRUE), nbins = length(communities))
)

enrichment_rows <- list()
bootstrap_rows <- list()
for (strategy_index in seq_len(nrow(strategy_definitions))) {
  definition <- strategy_definitions[strategy_index, ]
  if (definition$strategy == "Clinical OSTA <= -1") {
    point_referral <- baseline$clinical_osta
    point_cutoff <- -1
  } else {
    point_referral <- exact_rank_referral(baseline, baseline$risk_score, definition$capacity)
    point_cutoff <- min(baseline$risk_score[point_referral == 1], na.rm = TRUE)
  }
  point <- enrichment_metrics(baseline$incident_hip, point_referral)
  bootstrap_metrics <- matrix(
    NA_real_, nrow = bootstrap_replicates, ncol = length(point),
    dimnames = list(NULL, names(point))
  )
  for (replicate_index in seq_len(bootstrap_replicates)) {
    replicate_weight <- bootstrap_cluster_counts[cluster_index, replicate_index]
    if (definition$strategy == "Clinical OSTA <= -1") {
      replicate_referral <- baseline$clinical_osta
    } else {
      replicate_referral <- capacity_rule(
        baseline$risk_score, replicate_weight, definition$capacity
      )$referral
    }
    bootstrap_metrics[replicate_index, ] <- enrichment_metrics(
      baseline$incident_hip, replicate_referral, replicate_weight
    )
  }
  interval_list <- lapply(seq_along(point), function(metric_index) {
    bootstrap_interval(bootstrap_metrics[, metric_index])
  })
  names(interval_list) <- names(point)
  enrichment_rows[[strategy_index]] <- tibble(
    strategy = definition$strategy,
    intended_capacity = definition$capacity,
    cutoff = point_cutoff,
    n = nrow(baseline),
    events = sum(baseline$incident_hip == 1L),
    bootstrap_replicates = bootstrap_replicates,
    referral_rate = point[["referral_rate"]],
    referral_rate_lower = interval_list$referral_rate[1],
    referral_rate_upper = interval_list$referral_rate[2],
    event_capture = point[["event_capture"]],
    event_capture_lower = interval_list$event_capture[1],
    event_capture_upper = interval_list$event_capture[2],
    event_risk_referred = point[["event_risk_referred"]],
    event_risk_nonreferred = point[["event_risk_nonreferred"]],
    risk_ratio = point[["risk_ratio"]],
    risk_ratio_lower = interval_list$risk_ratio[1],
    risk_ratio_upper = interval_list$risk_ratio[2],
    risk_difference = point[["risk_difference"]],
    number_needed_to_screen = point[["number_needed_to_screen"]],
    specificity = point[["specificity"]]
  )
  bootstrap_rows[[strategy_index]] <- as_tibble(bootstrap_metrics) %>%
    mutate(replicate = row_number(), strategy = definition$strategy, .before = 1)
}
enrichment <- bind_rows(enrichment_rows)
enrichment_bootstrap <- bind_rows(bootstrap_rows)

# Primary complete-case regression set. Age and weight are not adjusted because they define OSTA.
model_data <- baseline %>%
  filter(complete.cases(
    followup_years, incident_hip, risk_z, high20, clinical_osta,
    baseline_fall, rural, stroke, arthritis, current_smoking, ever_drinking,
    community_id, bmi
  ))

cox_continuous <- coxph(
  Surv(followup_years, incident_hip) ~ risk_z + baseline_fall + rural + stroke + arthritis,
  data = model_data,
  ties = "efron",
  cluster = community_id
)
cox_high20 <- coxph(
  Surv(followup_years, incident_hip) ~ high20 + baseline_fall + rural + stroke + arthritis,
  data = model_data,
  ties = "efron",
  cluster = community_id
)
cox_lifestyle <- coxph(
  Surv(followup_years, incident_hip) ~ risk_z + bmi + baseline_fall + rural + stroke +
    arthritis + current_smoking + ever_drinking,
  data = model_data,
  ties = "efron",
  cluster = community_id
)
cox_clinical <- coxph(
  Surv(followup_years, incident_hip) ~ clinical_osta + baseline_fall + rural + stroke + arthritis,
  data = model_data,
  ties = "efron",
  cluster = community_id
)

regression_rows <- list(
  extract_cox_term(cox_continuous, "risk_z", "Cause-specific Cox", "OSTA risk per SD"),
  extract_cox_term(cox_high20, "high20", "Cause-specific Cox", "Top 20% OSTA risk"),
  extract_cox_term(cox_lifestyle, "risk_z", "Cox, lifestyle/BMI sensitivity", "OSTA risk per SD"),
  extract_cox_term(cox_clinical, "clinical_osta", "Cause-specific Cox", "Clinical OSTA <= -1")
)

poisson_fit <- glm(
  incident_hip ~ high20 + baseline_fall + rural + stroke + arthritis + offset(log(followup_years)),
  data = model_data,
  family = poisson(link = "log")
)
regression_rows[[length(regression_rows) + 1L]] <- extract_glm_term(
  poisson_fit,
  "high20",
  model_data$community_id,
  "Person-time Poisson",
  "Top 20% OSTA risk",
  "Incidence-rate ratio"
)

person_period <- build_person_period(model_data)
cloglog_fit <- glm(
  event_interval ~ risk_z + factor(interval_index) + baseline_fall + rural + stroke + arthritis,
  data = person_period,
  family = binomial(link = "cloglog")
)
regression_rows[[length(regression_rows) + 1L]] <- extract_glm_term(
  cloglog_fit,
  "risk_z",
  person_period$id,
  "Discrete-time complementary log-log",
  "OSTA risk per SD",
  "Discrete-time hazard ratio"
)

# Incremental concordance beyond baseline clinical factors.
concordance_data <- model_data %>%
  filter(complete.cases(followup_years, incident_hip, baseline_fall, rural, stroke, arthritis, risk_z))
cox_base <- coxph(
  Surv(followup_years, incident_hip) ~ baseline_fall + rural + stroke + arthritis,
  data = concordance_data,
  ties = "efron"
)
cox_base_osta <- coxph(
  Surv(followup_years, incident_hip) ~ risk_z + baseline_fall + rural + stroke + arthritis,
  data = concordance_data,
  ties = "efron"
)
concordance_summary <- tibble(
  model = c("Baseline clinical factors", "Baseline factors plus OSTA risk"),
  harrell_c = c(summary(cox_base)$concordance[1], summary(cox_base_osta)$concordance[1]),
  harrell_c_se = c(summary(cox_base)$concordance[2], summary(cox_base_osta)$concordance[2])
) %>%
  mutate(delta_c_vs_baseline = harrell_c - first(harrell_c))

ph_test <- cox.zph(cox_continuous, transform = "km")
ph_diagnostics <- tibble(
  term = rownames(ph_test$table),
  chisq = ph_test$table[, "chisq"],
  p_value = ph_test$table[, "p"]
)

# Competing mortality: cause-specific Cox and Fine-Gray.
competing_data <- competing %>%
  filter(complete.cases(
    followup_years, competing_status, risk_z, baseline_fall, rural, stroke,
    arthritis, community_id
  ))
cause_specific_fit <- coxph(
  Surv(followup_years, competing_status == 1L) ~ risk_z + baseline_fall + rural + stroke + arthritis,
  data = competing_data,
  ties = "efron",
  cluster = community_id
)
competing_rows <- list(
  extract_cox_term(
    cause_specific_fit, "risk_z", "Cause-specific Cox with competing deaths",
    "OSTA risk per SD"
  )
)
fine_gray_matrix <- model.matrix(
  ~ risk_z + baseline_fall + rural + stroke + arthritis,
  data = competing_data
)[, -1, drop = FALSE]
fine_gray_fit <- cmprsk::crr(
  ftime = competing_data$followup_years,
  fstatus = competing_data$competing_status,
  cov1 = fine_gray_matrix,
  failcode = 1,
  cencode = 0
)
fg_index <- which(names(fine_gray_fit$coef) == "risk_z")
fg_estimate <- fine_gray_fit$coef[fg_index]
fg_se <- sqrt(diag(fine_gray_fit$var))[fg_index]
competing_rows[[2]] <- tibble(
  analysis = "Fine-Gray subdistribution model",
  exposure = "OSTA risk per SD",
  estimate_scale = "Subdistribution hazard ratio",
  estimate = exp(fg_estimate),
  lower = exp(fg_estimate - 1.96 * fg_se),
  upper = exp(fg_estimate + 1.96 * fg_se),
  p_value = 2 * pnorm(-abs(fg_estimate / fg_se)),
  model_n = nrow(competing_data),
  model_events = sum(competing_data$competing_status == 1L),
  competing_deaths = sum(competing_data$competing_status == 2L)
)

# Time-varying inverse probability of censoring weights in the discrete-time model.
ipcw_data <- person_period %>%
  filter(complete.cases(
    censor_interval, interval_index, risk_z, baseline_fall, rural, stroke, arthritis
  ))
censor_marginal <- glm(
  censor_interval ~ factor(interval_index),
  data = ipcw_data,
  family = binomial()
)
censor_conditional <- glm(
  censor_interval ~ factor(interval_index) + risk_z + baseline_fall + rural + stroke + arthritis,
  data = ipcw_data,
  family = binomial()
)
ipcw_data <- ipcw_data %>%
  mutate(
    p_censor_marginal = pmin(pmax(predict(censor_marginal, type = "response"), 1e-4), 1 - 1e-4),
    p_censor_conditional = pmin(pmax(predict(censor_conditional, type = "response"), 1e-4), 1 - 1e-4)
  ) %>%
  group_by(id) %>%
  arrange(interval_index, .by_group = TRUE) %>%
  mutate(
    numerator_survival = lag(cumprod(1 - p_censor_marginal), default = 1),
    denominator_survival = lag(cumprod(1 - p_censor_conditional), default = 1),
    ipcw_raw = numerator_survival / denominator_survival
  ) %>%
  ungroup()
ipcw_limits <- quantile(ipcw_data$ipcw_raw, c(0.01, 0.99), na.rm = TRUE)
ipcw_data <- ipcw_data %>%
  mutate(ipcw = pmin(pmax(ipcw_raw, ipcw_limits[1]), ipcw_limits[2]))
ipcw_fit <- glm(
  event_interval ~ risk_z + factor(interval_index) + baseline_fall + rural + stroke + arthritis,
  data = ipcw_data,
  weights = ipcw,
  family = binomial(link = "cloglog")
)
competing_rows[[3]] <- extract_glm_term(
  ipcw_fit,
  "risk_z",
  ipcw_data$id,
  "Discrete-time model with stabilized IPCW",
  "OSTA risk per SD",
  "Discrete-time hazard ratio"
) %>%
  mutate(
    competing_deaths = NA_integer_,
    ipcw_lower = ipcw_limits[1],
    ipcw_upper = ipcw_limits[2]
  )

# Proxy-augmented hip-fracture outcome.
augmented_data <- baseline %>%
  filter(incident_hip_augmented %in% c(0, 1)) %>%
  mutate(
    augmented_event_year = case_when(
      incident_hip == 1L ~ incident_hip_year,
      incident_hip_augmented == 1L & death_proxy_hip_fracture == 1L ~ death_year,
      TRUE ~ NA_real_
    ),
    augmented_end_year = if_else(
      incident_hip_augmented == 1L,
      augmented_event_year,
      followup_end_year
    ),
    augmented_followup = augmented_end_year - baseline_year
  ) %>%
  filter(
    augmented_followup > 0,
    complete.cases(risk_z, baseline_fall, rural, stroke, arthritis, community_id)
  )
augmented_fit <- coxph(
  Surv(augmented_followup, incident_hip_augmented) ~ risk_z + baseline_fall + rural + stroke + arthritis,
  data = augmented_data,
  ties = "efron",
  cluster = community_id
)
competing_rows[[4]] <- extract_cox_term(
  augmented_fit, "risk_z", "Proxy-augmented hip-fracture outcome", "OSTA risk per SD"
) %>% mutate(competing_deaths = NA_integer_)

regression_results <- bind_rows(regression_rows)
robustness_results <- bind_rows(competing_rows)

# Sequential landmarks at waves 1-3; waves 4-5 have no usable OSTA.
landmarks <- build_landmarks(long_data)
landmark_rows <- list()
for (landmark_year in sort(unique(landmarks$landmark_year))) {
  landmark_data <- landmarks %>%
    filter(
      landmark_year == .env$landmark_year,
      complete.cases(followup_years, event, high20, baseline_fall, rural, stroke, arthritis, community_id)
    )
  fit <- coxph(
    Surv(followup_years, event) ~ high20 + baseline_fall + rural + stroke + arthritis,
    data = landmark_data,
    ties = "efron",
    cluster = community_id
  )
  landmark_rows[[length(landmark_rows) + 1L]] <- extract_cox_term(
    fit, "high20", paste0("Sequential landmark ", landmark_year), "Top 20% OSTA risk"
  ) %>% mutate(landmark_year = landmark_year)
}
landmark_complete <- landmarks %>%
  filter(complete.cases(followup_years, event, high20, baseline_fall, rural, stroke, arthritis, id))
landmark_pooled_fit <- coxph(
  Surv(followup_years, event) ~ high20 + strata(landmark_year) + baseline_fall + rural + stroke + arthritis,
  data = landmark_complete,
  ties = "efron",
  cluster = id
)
landmark_rows[[length(landmark_rows) + 1L]] <- extract_cox_term(
  landmark_pooled_fit,
  "high20",
  "Pooled sequential landmarks, stratified by year",
  "Top 20% OSTA risk"
) %>% mutate(landmark_year = NA_real_)
landmark_results <- bind_rows(landmark_rows)

# Supporting future-fall outcome.
fall_data <- baseline %>%
  filter(
    future_fall %in% c(0, 1),
    complete.cases(risk_z, baseline_fall, rural, stroke, arthritis, community_id)
  )
fall_fit <- glm(
  future_fall ~ risk_z + baseline_fall + rural + stroke + arthritis,
  data = fall_data,
  family = binomial()
)
fall_result <- extract_glm_term(
  fall_fit,
  "risk_z",
  fall_data$community_id,
  "Supporting future-fall outcome",
  "OSTA risk per SD",
  "Odds ratio"
)

# Figures.
capture_plot_data <- enrichment %>%
  filter(!is.na(intended_capacity))
capture_plot <- ggplot(
  capture_plot_data,
  aes(x = intended_capacity * 100, y = event_capture * 100)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#7A7F85") +
  geom_ribbon(
    aes(ymin = event_capture_lower * 100, ymax = event_capture_upper * 100),
    fill = "#B7D7D3", alpha = 0.55
  ) +
  geom_line(color = "#19726B", linewidth = 1) +
  geom_point(color = "#19726B", size = 3) +
  scale_x_continuous(breaks = capacities * 100) +
  labs(
    title = "Future hip fractures captured by OSTA-ranked referral",
    subtitle = "CHARLS 2011-2020; ribbons are community-cluster bootstrap 95% intervals",
    x = "Intended DXA referral capacity (%)",
    y = "Hip-fracture events captured (%)"
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    plot.title = element_text(face = "bold", color = "#1F2428"),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(output_dir, "Figure9_charls_event_capture.png"), capture_plot, width = 8.5, height = 5.8, dpi = 320)

forest_data <- bind_rows(
  regression_results %>% filter(exposure != "Clinical OSTA <= -1"),
  robustness_results,
  landmark_results %>% filter(is.na(landmark_year)),
  fall_result
) %>%
  mutate(analysis = factor(analysis, levels = rev(unique(analysis))))
forest_plot <- ggplot(forest_data, aes(x = estimate, y = analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#7A7F85") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.14, orientation = "y", color = "#245E9A") +
  geom_point(size = 2.8, color = "#245E9A") +
  scale_x_log10() +
  labs(
    title = "Prospective association of OSTA risk with clinical outcomes",
    subtitle = "Community-cluster robust estimates unless otherwise stated",
    x = "Relative effect estimate (log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 10.8) +
  theme(
    plot.title = element_text(face = "bold", color = "#1F2428"),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(output_dir, "Figure10_charls_forest.png"), forest_plot, width = 9.2, height = 7.0, dpi = 320)

landmark_plot_data <- landmark_results %>% filter(!is.na(landmark_year))
landmark_plot <- ggplot(landmark_plot_data, aes(x = estimate, y = factor(landmark_year))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#7A7F85") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.16, orientation = "y", color = "#9A3E50") +
  geom_point(size = 3, color = "#9A3E50") +
  scale_x_log10() +
  labs(
    title = "Sequential-landmark validation of the top 20% OSTA rank",
    subtitle = "Separate baselines in 2011, 2013, and 2015",
    x = "Hazard ratio (log scale)",
    y = "Landmark year"
  ) +
  theme_minimal(base_size = 11.5) +
  theme(plot.title = element_text(face = "bold", color = "#1F2428"), panel.grid.minor = element_blank())
ggsave(file.path(output_dir, "Figure11_charls_landmarks.png"), landmark_plot, width = 8, height = 5.5, dpi = 320)

qc <- bind_rows(
  enrichment %>% transmute(
    check = paste(strategy, "finite enrichment metrics"),
    passed = is.finite(event_capture) & is.finite(risk_ratio) & is.finite(number_needed_to_screen),
    severity = "error",
    value = event_capture
  ),
  regression_results %>% transmute(
    check = paste(analysis, exposure, "finite regression estimate"),
    passed = is.finite(estimate) & is.finite(lower) & is.finite(upper) & lower > 0,
    severity = "error",
    value = estimate
  ),
  robustness_results %>% transmute(
    check = paste(analysis, "finite robustness estimate"),
    passed = is.finite(estimate) & is.finite(lower) & is.finite(upper) & lower > 0,
    severity = "error",
    value = estimate
  ),
  landmark_results %>% transmute(
    check = paste(analysis, "finite landmark estimate"),
    passed = is.finite(estimate) & is.finite(lower) & is.finite(upper) & lower > 0,
    severity = "error",
    value = estimate
  ),
  tibble(
    check = "Proportional-hazards test for continuous OSTA risk",
    passed = ph_diagnostics$p_value[ph_diagnostics$term == "risk_z"] >= 0.05,
    severity = "warning",
    value = ph_diagnostics$p_value[ph_diagnostics$term == "risk_z"]
  )
)

write.csv(enrichment, file.path(output_dir, "Table10_charls_capacity_enrichment.csv"), row.names = FALSE, na = "")
write.csv(regression_results, file.path(output_dir, "Table11_charls_regression_models.csv"), row.names = FALSE, na = "")
write.csv(robustness_results, file.path(output_dir, "Table12_charls_competing_ipcw_proxy.csv"), row.names = FALSE, na = "")
write.csv(concordance_summary, file.path(output_dir, "TableS9_incremental_concordance.csv"), row.names = FALSE, na = "")
write.csv(ph_diagnostics, file.path(output_dir, "TableS10_proportional_hazards.csv"), row.names = FALSE, na = "")
write.csv(landmark_results, file.path(output_dir, "TableS11_sequential_landmarks.csv"), row.names = FALSE, na = "")
write.csv(fall_result, file.path(output_dir, "TableS12_future_falls.csv"), row.names = FALSE, na = "")
write.csv(enrichment_bootstrap, file.path(output_dir, "TableS13_enrichment_bootstrap_replicates.csv"), row.names = FALSE, na = "")
write.csv(qc, file.path(output_dir, "Step6_quality_control.csv"), row.names = FALSE, na = "")
saveRDS(
  list(
    baseline = baseline,
    model_data = model_data,
    person_period = person_period,
    competing_data = competing_data,
    landmarks = landmarks,
    fitted_models = list(
      cox_continuous = cox_continuous,
      cox_high20 = cox_high20,
      cloglog = cloglog_fit,
      cause_specific = cause_specific_fit,
      fine_gray = fine_gray_fit,
      ipcw = ipcw_fit,
      landmark_pooled = landmark_pooled_fit
    )
  ),
  file.path(output_dir, "Step6_analysis_objects.rds")
)
capture.output(sessionInfo(), file = file.path(output_dir, "Step6_sessionInfo.txt"))

cat("\nSTEP 6 COMPLETED\n")
cat("Output:", output_dir, "\n")
print(enrichment)
print(regression_results)
print(robustness_results)
print(landmark_results)
