# STEP 7: robustness analyses, cross-cohort synthesis, and final go/no-go assessment

required_packages <- c(
  "dplyr", "tidyr", "ggplot2", "survey", "survival", "mice",
  "splines", "tibble"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages before running Step 7: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survey)
  library(survival)
  library(mice)
  library(splines)
  library(tibble)
})

options(survey.lonely.psu = "adjust")
set.seed(20260805)

project_dir <- Sys.getenv("DXA_PROJECT_DIR", unset = file.path(getwd(), "DXA_transportability_results"))
input_dir <- file.path(project_dir, "data_clean")
output_dir <- file.path(project_dir, "step7_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_replicates <- as.integer(Sys.getenv("STEP7_BOOTSTRAP_REPS", unset = "200"))
if (!is.finite(bootstrap_replicates) || bootstrap_replicates < 20) bootstrap_replicates <- 200L
imputation_count <- as.integer(Sys.getenv("STEP7_IMPUTATIONS", unset = "20"))
if (!is.finite(imputation_count) || imputation_count < 5) imputation_count <- 20L

weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

weighted_auc <- function(y, score, w) {
  keep <- y %in% c(0, 1) & is.finite(score) & is.finite(w) & w > 0
  y <- y[keep]
  score <- score[keep]
  w <- w[keep]
  if (length(unique(y)) < 2) return(NA_real_)
  ordered <- order(score)
  y <- y[ordered]
  score <- score[ordered]
  w <- w[ordered]
  groups <- split(seq_along(score), cumsum(c(TRUE, diff(score) != 0)))
  negative_before <- 0
  concordant <- 0
  for (index in groups) {
    positive_weight <- sum(w[index][y[index] == 1])
    negative_weight <- sum(w[index][y[index] == 0])
    concordant <- concordant + positive_weight * (negative_before + 0.5 * negative_weight)
    negative_before <- negative_before + negative_weight
  }
  concordant / (sum(w[y == 1]) * sum(w[y == 0]))
}

make_design <- function(data) {
  design_data <- data %>% mutate(survey_strata = factor(strata), survey_psu = factor(psu))
  svydesign(
    ids = ~survey_psu,
    strata = ~survey_strata,
    weights = ~survey_weight,
    nest = TRUE,
    data = design_data
  )
}

make_replicate_weights <- function(data, replicates) {
  replicate_design <- as.svrepdesign(
    make_design(data),
    type = "bootstrap",
    replicates = replicates,
    mse = TRUE
  )
  as.matrix(weights(replicate_design, type = "analysis"))
}

bootstrap_auc_interval <- function(y, score, point_weight, replicate_weight) {
  point <- weighted_auc(y, score, point_weight)
  replicate_auc <- apply(
    replicate_weight,
    2,
    function(w) weighted_auc(y, score, w)
  )
  valid <- replicate_auc[is.finite(replicate_auc)]
  interval <- if (length(valid) >= 10) {
    quantile(valid, c(0.025, 0.975), names = FALSE, type = 6)
  } else {
    c(NA_real_, NA_real_)
  }
  c(auc = point, lower = interval[1], upper = interval[2], valid_replicates = length(valid))
}

extract_cox_term <- function(fit, term, analysis, subgroup = "All eligible") {
  coefficient_table <- summary(fit)$coefficients
  if (!term %in% rownames(coefficient_table)) stop("Term not found: ", term)
  log_estimate <- coefficient_table[term, "coef"]
  standard_error_column <- if ("robust se" %in% colnames(coefficient_table)) "robust se" else "se(coef)"
  standard_error <- coefficient_table[term, standard_error_column]
  tibble(
    analysis = analysis,
    subgroup = subgroup,
    estimate = exp(log_estimate),
    lower = exp(log_estimate - 1.96 * standard_error),
    upper = exp(log_estimate + 1.96 * standard_error),
    p_value = 2 * pnorm(-abs(log_estimate / standard_error)),
    model_n = fit$n,
    model_events = fit$nevent
  )
}

rubin_pool_scalar <- function(estimates, variances) {
  m <- length(estimates)
  q_bar <- mean(estimates)
  u_bar <- mean(variances)
  between <- var(estimates)
  total <- u_bar + (1 + 1 / m) * between
  if (!is.finite(between) || between <= 0) {
    degrees_freedom <- Inf
  } else {
    degrees_freedom <- (m - 1) * (1 + u_bar / ((1 + 1 / m) * between))^2
  }
  critical <- if (is.finite(degrees_freedom)) qt(0.975, degrees_freedom) else qnorm(0.975)
  tibble(
    log_estimate = q_bar,
    standard_error = sqrt(total),
    degrees_freedom = degrees_freedom,
    estimate = exp(q_bar),
    lower = exp(q_bar - critical * sqrt(total)),
    upper = exp(q_bar + critical * sqrt(total)),
    p_value = 2 * pt(-abs(q_bar / sqrt(total)), df = degrees_freedom)
  )
}

required_files <- c(
  file.path(input_dir, "imaging_cross_country_analysis_ready.rds"),
  file.path(input_dir, "charls_analysis_ready.rds"),
  file.path(project_dir, "step3_results", "Table4_model_performance.csv"),
  file.path(project_dir, "step4_results", "Table5_risk_scale_decomposition.csv"),
  file.path(project_dir, "step5_results", "Table7_capacity_strategy_performance.csv"),
  file.path(project_dir, "step6_results", "Table10_charls_capacity_enrichment.csv"),
  file.path(project_dir, "step6_results", "Table11_charls_regression_models.csv"),
  file.path(project_dir, "step6_results", "Table12_charls_competing_ipcw_proxy.csv"),
  file.path(project_dir, "step6_results", "TableS11_sequential_landmarks.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) stop("Required prior-stage files are missing: ", paste(missing_files, collapse = ", "))

imaging <- readRDS(file.path(input_dir, "imaging_cross_country_analysis_ready.rds")) %>%
  mutate(
    period = case_when(
      analysis_role == "development" ~ "NHANES 2005-2010",
      analysis_role == "temporal_validation_1" ~ "NHANES 2013-2014",
      analysis_role == "temporal_validation_2" ~ "NHANES 2017-Mar 2020",
      analysis_role == "historical_geographic_validation" ~ "KNHANES 2008-2011",
      analysis_role == "contemporary_temporal_validation" ~ "KNHANES 2024",
      TRUE ~ NA_character_
    ),
    risk_score = -osta_raw
  ) %>%
  filter(!is.na(period), complete.cases(risk_score, survey_weight, strata, psu), survey_weight > 0)
period_data <- split(imaging, imaging$period)
period_replicates <- lapply(period_data, make_replicate_weights, replicates = bootstrap_replicates)

outcome_sites <- tribble(
  ~site, ~local_suffix,
  "Femoral neck", "osteoporosis_fn",
  "Total hip", "osteoporosis_total_hip",
  "Lumbar spine", "osteoporosis_lumbar",
  "Any complete site", "osteoporosis_any_complete"
)
subgroup_definitions <- c(
  "All eligible", "Age 50-64", "Age 65-79", "BMI <25 kg/m2",
  "BMI >=25 kg/m2", "Postmenopausal observed"
)

imaging_rows <- list()
for (period_name in names(period_data)) {
  data <- period_data[[period_name]]
  replicate_weight <- period_replicates[[period_name]]
  reference_definitions <- if (grepl("KNHANES", period_name, fixed = TRUE)) {
    c("local", "us")
  } else {
    "us"
  }
  for (reference in reference_definitions) {
    reference_label <- ifelse(reference == "local", "Asian local", "US common")
    for (site_index in seq_len(nrow(outcome_sites))) {
      outcome_variable <- paste0(outcome_sites$local_suffix[site_index], "_", reference)
      for (subgroup in subgroup_definitions) {
        keep <- switch(
          subgroup,
          "All eligible" = rep(TRUE, nrow(data)),
          "Age 50-64" = data$age < 65,
          "Age 65-79" = data$age >= 65,
          "BMI <25 kg/m2" = data$bmi < 25,
          "BMI >=25 kg/m2" = data$bmi >= 25,
          "Postmenopausal observed" = data$postmenopausal == 1,
          rep(FALSE, nrow(data))
        )
        keep <- keep & !is.na(keep) & data[[outcome_variable]] %in% c(0, 1)
        y <- data[[outcome_variable]][keep]
        if (length(y) < 100 || sum(y == 1) < 10 || length(unique(y)) < 2) next
        metrics <- bootstrap_auc_interval(
          y,
          data$risk_score[keep],
          data$survey_weight[keep],
          replicate_weight[keep, , drop = FALSE]
        )
        imaging_rows[[length(imaging_rows) + 1L]] <- tibble(
          period = period_name,
          cohort = first(data$cohort),
          outcome_reference = reference_label,
          site = outcome_sites$site[site_index],
          outcome = outcome_variable,
          subgroup = subgroup,
          n = sum(keep),
          events = sum(y == 1),
          weighted_prevalence = weighted_mean(y, data$survey_weight[keep]),
          auc = metrics[["auc"]],
          auc_lower = metrics[["lower"]],
          auc_upper = metrics[["upper"]],
          valid_bootstrap_replicates = metrics[["valid_replicates"]],
          bootstrap_replicates = bootstrap_replicates
        )
      }
    }
  }
}
imaging_robustness <- bind_rows(imaging_rows)

# Cycle-specific primary-outcome stability.
cycle_rows <- list()
cycle_groups <- split(imaging, interaction(imaging$period, imaging$cycle, drop = TRUE))
for (cycle_data in cycle_groups) {
  period_name <- first(cycle_data$period)
  outcome_variable <- ifelse(
    grepl("KNHANES", period_name, fixed = TRUE),
    "osteoporosis_fn_local",
    "osteoporosis_fn_us"
  )
  keep <- cycle_data[[outcome_variable]] %in% c(0, 1)
  if (sum(keep) < 100 || sum(cycle_data[[outcome_variable]][keep] == 1) < 8) next
  cycle_rows[[length(cycle_rows) + 1L]] <- tibble(
    period = period_name,
    cycle = first(cycle_data$cycle),
    n = sum(keep),
    events = sum(cycle_data[[outcome_variable]][keep] == 1),
    weighted_prevalence = weighted_mean(
      cycle_data[[outcome_variable]][keep], cycle_data$survey_weight[keep]
    ),
    auc = weighted_auc(
      cycle_data[[outcome_variable]][keep], cycle_data$risk_score[keep], cycle_data$survey_weight[keep]
    )
  )
}
cycle_stability <- bind_rows(cycle_rows)

# CHARLS robustness: complete-case, sampling-weight, BMI restriction, subgroups, and interactions.
charls <- readRDS(file.path(input_dir, "charls_analysis_ready.rds")) %>%
  mutate(
    risk_score = -osta_raw,
    risk_z = as.numeric(scale(risk_score)),
    age_group = factor(if_else(age >= 65, "65-79", "50-64")),
    bmi_group = factor(if_else(bmi >= 25, ">=25", "<25"))
  )

charls_complete <- charls %>%
  filter(complete.cases(
    followup_years, incident_hip, risk_z, baseline_fall, rural, stroke, arthritis,
    current_smoking, ever_drinking, bmi, community_id
  ))

charls_sensitivity_rows <- list()
complete_fit <- coxph(
  Surv(followup_years, incident_hip) ~ risk_z + baseline_fall + rural + stroke + arthritis,
  data = charls_complete,
  ties = "efron",
  cluster = community_id
)
charls_sensitivity_rows[[1]] <- extract_cox_term(
  complete_fit, "risk_z", "Complete-case primary model"
)

weighted_data <- charls_complete %>%
  filter(is.finite(individual_weight), individual_weight > 0) %>%
  mutate(analysis_weight = individual_weight / mean(individual_weight))
weighted_fit <- coxph(
  Surv(followup_years, incident_hip) ~ risk_z + baseline_fall + rural + stroke + arthritis,
  data = weighted_data,
  weights = analysis_weight,
  ties = "efron",
  robust = TRUE,
  cluster = community_id
)
charls_sensitivity_rows[[2]] <- extract_cox_term(
  weighted_fit, "risk_z", "Individual-weight complete-case sensitivity"
)

bmi_restricted <- charls_complete %>% filter(bmi >= 15, bmi <= 40)
bmi_restricted_fit <- coxph(
  Surv(followup_years, incident_hip) ~ risk_z + baseline_fall + rural + stroke + arthritis,
  data = bmi_restricted,
  ties = "efron",
  cluster = community_id
)
charls_sensitivity_rows[[3]] <- extract_cox_term(
  bmi_restricted_fit, "risk_z", "BMI 15-40 kg/m2 sensitivity"
)

subgroups <- list(
  "Age 50-64" = charls_complete$age < 65,
  "Age 65-79" = charls_complete$age >= 65,
  "BMI <25 kg/m2" = charls_complete$bmi < 25,
  "BMI >=25 kg/m2" = charls_complete$bmi >= 25,
  "Rural" = charls_complete$rural == 1,
  "Urban" = charls_complete$rural == 0
)
for (subgroup_name in names(subgroups)) {
  subgroup_data <- charls_complete[subgroups[[subgroup_name]], , drop = FALSE]
  fit <- coxph(
    Surv(followup_years, incident_hip) ~ risk_z + baseline_fall + rural + stroke + arthritis,
    data = subgroup_data,
    ties = "efron",
    cluster = community_id
  )
  charls_sensitivity_rows[[length(charls_sensitivity_rows) + 1L]] <- extract_cox_term(
    fit, "risk_z", "Prespecified subgroup", subgroup_name
  )
}
charls_sensitivity <- bind_rows(charls_sensitivity_rows)

age_interaction_fit <- coxph(
  Surv(followup_years, incident_hip) ~ risk_z * age_group + baseline_fall + rural + stroke + arthritis,
  data = charls_complete,
  ties = "efron",
  cluster = community_id
)
bmi_interaction_fit <- coxph(
  Surv(followup_years, incident_hip) ~ risk_z * bmi_group + baseline_fall + rural + stroke + arthritis,
  data = charls_complete,
  ties = "efron",
  cluster = community_id
)
interaction_rows <- list()
for (definition in list(
  list(fit = age_interaction_fit, pattern = "risk_z:age_group", modifier = "Age group"),
  list(fit = bmi_interaction_fit, pattern = "risk_z:bmi_group", modifier = "BMI group")
)) {
  coefficient_table <- summary(definition$fit)$coefficients
  term <- grep(definition$pattern, rownames(coefficient_table), value = TRUE, fixed = TRUE)[1]
  log_estimate <- coefficient_table[term, "coef"]
  standard_error <- coefficient_table[term, "robust se"]
  interaction_rows[[length(interaction_rows) + 1L]] <- tibble(
    modifier = definition$modifier,
    interaction_term = term,
    interaction_ratio = exp(log_estimate),
    lower = exp(log_estimate - 1.96 * standard_error),
    upper = exp(log_estimate + 1.96 * standard_error),
    p_interaction = 2 * pnorm(-abs(log_estimate / standard_error))
  )
}
interactions <- bind_rows(interaction_rows)

# Multiple imputation of the small amount of missing baseline covariate information.
imputation_data <- charls %>%
  select(
    id, community_id, followup_years, incident_hip, risk_z, bmi,
    baseline_fall, rural, stroke, arthritis, current_smoking, ever_drinking
  )
methods <- mice::make.method(imputation_data)
methods[c("id", "community_id", "followup_years", "incident_hip", "risk_z")] <- ""
predictor_matrix <- mice::make.predictorMatrix(imputation_data)
diag(predictor_matrix) <- 0
predictor_matrix[, c("id", "community_id")] <- 0
imputed <- mice::mice(
  imputation_data,
  m = imputation_count,
  maxit = 10,
  method = methods,
  predictorMatrix = predictor_matrix,
  seed = 20260805,
  printFlag = FALSE
)
mi_estimates <- numeric(imputation_count)
mi_variances <- numeric(imputation_count)
for (imputation_index in seq_len(imputation_count)) {
  completed <- mice::complete(imputed, imputation_index)
  fit <- coxph(
    Surv(followup_years, incident_hip) ~ risk_z + bmi + baseline_fall + rural + stroke +
      arthritis + current_smoking + ever_drinking,
    data = completed,
    ties = "efron",
    cluster = community_id
  )
  mi_estimates[imputation_index] <- coef(fit)["risk_z"]
  mi_variances[imputation_index] <- vcov(fit)["risk_z", "risk_z"]
}
mi_result <- rubin_pool_scalar(mi_estimates, mi_variances) %>%
  mutate(
    analysis = "Multiple imputation with Rubin pooling",
    subgroup = "All eligible",
    imputations = imputation_count,
    model_n = nrow(imputation_data),
    model_events = sum(imputation_data$incident_hip == 1L),
    .before = 1
  )

# Prior-stage evidence synthesis and final go/no-go decision.
step3 <- read.csv(file.path(project_dir, "step3_results", "Table4_model_performance.csv"), check.names = FALSE)
step4 <- read.csv(file.path(project_dir, "step4_results", "Table5_risk_scale_decomposition.csv"), check.names = FALSE)
step5 <- read.csv(file.path(project_dir, "step5_results", "Table7_capacity_strategy_performance.csv"), check.names = FALSE)
step6_enrichment <- read.csv(file.path(project_dir, "step6_results", "Table10_charls_capacity_enrichment.csv"), check.names = FALSE)
step6_regression <- read.csv(file.path(project_dir, "step6_results", "Table11_charls_regression_models.csv"), check.names = FALSE)
step6_robust <- read.csv(file.path(project_dir, "step6_results", "Table12_charls_competing_ipcw_proxy.csv"), check.names = FALSE)
step6_landmark <- read.csv(file.path(project_dir, "step6_results", "TableS11_sequential_landmarks.csv"), check.names = FALSE)

validation_scenarios <- c(
  "NH temporal validation 1", "NH temporal validation 2",
  "KN historical geographic", "KN contemporary temporal-geographic",
  "Reverse KN-to-NH common reference", "Korean temporal validation"
)
validation_osta <- step3 %>% filter(model == "OSTA", scenario %in% validation_scenarios)
primary20 <- step6_enrichment %>% filter(strategy == "Target-percentile OSTA 20%")
primary_hr <- step6_regression %>% filter(
  analysis == "Cause-specific Cox", exposure == "Top 20% OSTA risk"
)
target_rank20 <- step5 %>% filter(
  outcome_scope == "Primary definition",
  intended_capacity == 0.20,
  strategy == "Target-percentile OSTA"
)

go_no_go <- tribble(
  ~criterion, ~threshold, ~observed, ~passed,
  "Bidirectional/temporal OSTA discrimination", "All prespecified validation AUC >=0.75",
  sprintf("Minimum AUC %.3f", min(validation_osta$auc, na.rm = TRUE)),
  min(validation_osta$auc, na.rm = TRUE) >= 0.75,
  "Clinically meaningful calibration or threshold drift", "Absolute calibration intercept >0.20 or target referral deviation >10 pp",
  sprintf("Maximum |intercept| %.3f; fixed-OSTA Korean deviation 30.9-33.4 pp",
          max(abs(validation_osta$calibration_intercept), na.rm = TRUE)),
  max(abs(validation_osta$calibration_intercept), na.rm = TRUE) > 0.20,
  "Capacity control", "Target-rank referral deviation <0.1 pp",
  sprintf("Maximum deviation %.3f pp", 100 * max(abs(target_rank20$capacity_deviation), na.rm = TRUE)),
  max(abs(target_rank20$capacity_deviation), na.rm = TRUE) < 0.001,
  "Prospective CHARLS clinical enrichment", "Top 20% captures >20% events and adjusted HR lower CI >1",
  sprintf("Captured %.1f%%; HR %.2f (%.2f-%.2f)",
          100 * primary20$event_capture, primary_hr$estimate, primary_hr$lower, primary_hr$upper),
  primary20$event_capture > 0.20 & primary_hr$lower > 1,
  "Robustness to competing risk, IPCW, proxy outcome, and landmarks", "All prespecified OSTA estimates remain directionally positive",
  sprintf("Minimum lower CI %.2f", min(c(step6_robust$lower, step6_landmark$lower), na.rm = TRUE)),
  min(c(step6_robust$lower, step6_landmark$lower), na.rm = TRUE) > 1,
  "Worst-case same-capacity regret", "OSTA maximum regret <=0.10",
  sprintf("Maximum regret %.3f", max(target_rank20$sensitivity_regret, na.rm = TRUE)),
  max(target_rank20$sensitivity_regret, na.rm = TRUE) <= 0.10
) %>%
  mutate(decision = if_else(passed, "GO", "NO-GO"))

final_decision <- if (all(go_no_go$passed)) "GO: proceed to full manuscript and submission preparation" else
  "CONDITIONAL/NO-GO: resolve failed criteria before submission"

# Figures.
imaging_plot_data <- imaging_robustness %>%
  filter(subgroup == "All eligible") %>%
  mutate(
    period = factor(period, levels = rev(c(
      "NHANES 2005-2010", "NHANES 2013-2014", "NHANES 2017-Mar 2020",
      "KNHANES 2008-2011", "KNHANES 2024"
    ))),
    site = factor(site, levels = c("Femoral neck", "Total hip", "Lumbar spine", "Any complete site"))
  )
imaging_plot <- ggplot(imaging_plot_data, aes(x = auc, y = period, color = site)) +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "#7A7F85") +
  geom_errorbar(
    aes(xmin = auc_lower, xmax = auc_upper),
    width = 0.12,
    orientation = "y",
    position = position_dodge(width = 0.62)
  ) +
  geom_point(size = 2.3, position = position_dodge(width = 0.62)) +
  facet_wrap(~outcome_reference) +
  scale_color_manual(values = c(
    "Femoral neck" = "#245E9A", "Total hip" = "#19726B",
    "Lumbar spine" = "#C07A32", "Any complete site" = "#9A3E50"
  )) +
  coord_cartesian(xlim = c(0.55, 0.92)) +
  labs(
    title = "OSTA discrimination across DXA sites and reference standards",
    subtitle = "Survey-weighted AUC with 200-replicate 95% intervals",
    x = "Area under the ROC curve",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 10.8) +
  theme(
    plot.title = element_text(face = "bold", color = "#1F2428"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
ggsave(file.path(output_dir, "Figure12_imaging_robustness.png"), imaging_plot, width = 10.5, height = 7.4, dpi = 320)

charls_plot_data <- bind_rows(
  charls_sensitivity %>% select(analysis, subgroup, estimate, lower, upper),
  mi_result %>% select(analysis, subgroup, estimate, lower, upper)
) %>%
  mutate(label = if_else(subgroup == "All eligible", analysis, subgroup), label = factor(label, levels = rev(unique(label))))
charls_plot <- ggplot(charls_plot_data, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#7A7F85") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.14, orientation = "y", color = "#19726B") +
  geom_point(size = 2.7, color = "#19726B") +
  scale_x_log10() +
  labs(
    title = "Robustness of the CHARLS OSTA association",
    subtitle = "Hazard ratio per 1-SD increase in OSTA risk",
    x = "Hazard ratio (log scale)",
    y = NULL
  ) +
  theme_minimal(base_size = 10.8) +
  theme(plot.title = element_text(face = "bold", color = "#1F2428"), panel.grid.minor = element_blank())
ggsave(file.path(output_dir, "Figure13_charls_robustness.png"), charls_plot, width = 8.8, height = 6.6, dpi = 320)

qc <- bind_rows(
  imaging_robustness %>% transmute(
    check = paste(period, outcome_reference, site, subgroup, "finite AUC", sep = ": "),
    passed = is.finite(auc) & between(auc, 0, 1) & is.finite(auc_lower) & is.finite(auc_upper),
    severity = "error",
    value = auc
  ),
  charls_sensitivity %>% transmute(
    check = paste(analysis, subgroup, "finite CHARLS estimate", sep = ": "),
    passed = is.finite(estimate) & is.finite(lower) & is.finite(upper) & lower > 0,
    severity = "error",
    value = estimate
  ),
  go_no_go %>% transmute(
    check = paste("Go/no-go", criterion),
    passed = passed,
    severity = "decision",
    value = as.numeric(passed)
  )
)

write.csv(imaging_robustness, file.path(output_dir, "Table13_imaging_site_subgroup_robustness.csv"), row.names = FALSE, na = "")
write.csv(cycle_stability, file.path(output_dir, "TableS14_cycle_stability.csv"), row.names = FALSE, na = "")
write.csv(charls_sensitivity, file.path(output_dir, "Table14_charls_weight_subgroup_sensitivity.csv"), row.names = FALSE, na = "")
write.csv(interactions, file.path(output_dir, "TableS15_charls_interactions.csv"), row.names = FALSE, na = "")
write.csv(mi_result, file.path(output_dir, "TableS16_multiple_imputation.csv"), row.names = FALSE, na = "")
write.csv(go_no_go, file.path(output_dir, "Table15_final_go_no_go.csv"), row.names = FALSE, na = "")
writeLines(final_decision, file.path(output_dir, "FINAL_DECISION.txt"))
write.csv(qc, file.path(output_dir, "Step7_quality_control.csv"), row.names = FALSE, na = "")
saveRDS(
  list(
    imputed_data = imputed,
    age_interaction_fit = age_interaction_fit,
    bmi_interaction_fit = bmi_interaction_fit,
    imaging_robustness = imaging_robustness,
    go_no_go = go_no_go
  ),
  file.path(output_dir, "Step7_analysis_objects.rds")
)
capture.output(sessionInfo(), file = file.path(output_dir, "Step7_sessionInfo.txt"))

cat("\nSTEP 7 COMPLETED\n")
cat("Output:", output_dir, "\n")
cat("Decision:", final_decision, "\n")
print(go_no_go)
print(charls_sensitivity)
print(mi_result)
