# STEP 4: calibration-drift decomposition and target-cohort recalibration
# OSTA is primary; the age-weight spline is a prespecified sensitivity model.

required_packages <- c("dplyr", "tidyr", "ggplot2", "survey", "splines", "glmnet", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages before running Step 4: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survey)
  library(splines)
  library(glmnet)
  library(tibble)
})

options(survey.lonely.psu = "adjust")
set.seed(20260805)

project_dir <- Sys.getenv("DXA_PROJECT_DIR", unset = file.path(getwd(), "DXA_transportability_results"))
input_file <- file.path(project_dir, "data_clean", "imaging_cross_country_analysis_ready.rds")
model_file <- file.path(project_dir, "step3_results", "Step3_fitted_models.rds")
output_dir <- file.path(project_dir, "step4_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

weighted_var <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (sum(keep) < 2) return(NA_real_)
  x <- x[keep]
  w <- w[keep]
  mu <- sum(x * w) / sum(w)
  sum(w * (x - mu)^2) / sum(w)
}

weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cumulative <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cumulative >= p)[1]], numeric(1))
}

kish_n <- function(w) {
  keep <- is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(w[keep])^2 / sum(w[keep]^2)
}

weighted_auc <- function(y, score, w) {
  keep <- y %in% c(0, 1) & is.finite(score) & is.finite(w) & w > 0
  y <- y[keep]
  score <- score[keep]
  w <- w[keep]
  if (length(unique(y)) < 2) return(NA_real_)
  ord <- order(score)
  y <- y[ord]
  score <- score[ord]
  w <- w[ord]
  groups <- split(seq_along(score), cumsum(c(TRUE, diff(score) != 0)))
  negative_before <- 0
  concordant <- 0
  for (idx in groups) {
    positive_weight <- sum(w[idx][y[idx] == 1])
    negative_weight <- sum(w[idx][y[idx] == 0])
    concordant <- concordant + positive_weight * (negative_before + 0.5 * negative_weight)
    negative_before <- negative_before + negative_weight
  }
  concordant / (sum(w[y == 1]) * sum(w[y == 0]))
}

weighted_pr_auc <- function(y, score, w) {
  keep <- y %in% c(0, 1) & is.finite(score) & is.finite(w) & w > 0
  y <- y[keep]
  score <- score[keep]
  w <- w[keep]
  if (length(unique(y)) < 2) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  score <- score[ord]
  w <- w[ord]
  groups <- split(seq_along(score), cumsum(c(TRUE, diff(score) != 0)))
  total_positive <- sum(w[y == 1])
  true_positive <- false_positive <- recall <- area <- 0
  for (idx in groups) {
    true_positive <- true_positive + sum(w[idx][y[idx] == 1])
    false_positive <- false_positive + sum(w[idx][y[idx] == 0])
    new_recall <- true_positive / total_positive
    precision <- true_positive / (true_positive + false_positive)
    area <- area + (new_recall - recall) * precision
    recall <- new_recall
  }
  area
}

make_design <- function(data, weight_variable = "survey_weight") {
  data <- data %>%
    mutate(survey_strata = factor(strata), survey_psu = factor(psu))
  formula_weight <- as.formula(paste0("~", weight_variable))
  svydesign(
    ids = ~survey_psu,
    strata = ~survey_strata,
    weights = formula_weight,
    nest = TRUE,
    data = data
  )
}

prevalence_summary <- function(data, outcome, weight_variable = "survey_weight") {
  analysis <- data %>%
    mutate(model_outcome = .data[[outcome]]) %>%
    filter(model_outcome %in% c(0, 1), is.finite(.data[[weight_variable]]), .data[[weight_variable]] > 0)
  estimate <- svymean(~model_outcome, make_design(analysis, weight_variable))
  value <- as.numeric(coef(estimate)[1])
  standard_error <- as.numeric(SE(estimate)[1])
  tibble(
    prevalence = value,
    prevalence_se = standard_error,
    prevalence_lower = max(0, value - 1.96 * standard_error),
    prevalence_upper = min(1, value + 1.96 * standard_error)
  )
}

standardized_difference <- function(source_x, source_w, target_x, target_w) {
  difference <- weighted_mean(source_x, source_w) - weighted_mean(target_x, target_w)
  denominator <- sqrt((weighted_var(source_x, source_w) + weighted_var(target_x, target_w)) / 2)
  if (!is.finite(denominator) || denominator == 0) return(NA_real_)
  difference / denominator
}

ridge_matrix <- function(data, bmi_median) {
  prepared <- data %>% mutate(bmi_imputed = if_else(is.finite(bmi), bmi, bmi_median))
  model.matrix(
    ~ age + I(age^2) + weight_kg + I(weight_kg^2) + bmi_imputed + age:weight_kg,
    data = prepared
  )[, -1, drop = FALSE]
}

predict_selected_models <- function(models, new_data) {
  list(
    OSTA = as.numeric(predict(models$osta, newdata = new_data, type = "response")),
    `Age-weight spline` = as.numeric(predict(models$spline, newdata = new_data, type = "response"))
  )
}

performance_metrics <- function(y, prediction, w) {
  keep <- y %in% c(0, 1) & is.finite(prediction) & is.finite(w) & w > 0
  y <- y[keep]
  prediction <- pmin(pmax(prediction[keep], 1e-6), 1 - 1e-6)
  w <- w[keep]
  tibble(
    auc = weighted_auc(y, prediction, w),
    pr_auc = weighted_pr_auc(y, prediction, w),
    brier = weighted_mean((y - prediction)^2, w)
  )
}

calibration_statistics <- function(data, prediction) {
  analysis <- data %>%
    mutate(
      prediction_internal = pmin(pmax(prediction, 1e-6), 1 - 1e-6),
      linear_predictor = qlogis(prediction_internal)
    ) %>%
    filter(model_outcome %in% c(0, 1), complete.cases(prediction_internal, linear_predictor))
  design <- make_design(analysis)
  intercept_fit <- svyglm(
    model_outcome ~ 1,
    offset = linear_predictor,
    design = design,
    family = quasibinomial()
  )
  slope_fit <- svyglm(
    model_outcome ~ linear_predictor,
    design = design,
    family = quasibinomial()
  )
  c(
    calibration_intercept = unname(coef(intercept_fit)[1]),
    calibration_intercept_se = unname(sqrt(vcov(intercept_fit)[1, 1])),
    calibration_slope = unname(coef(slope_fit)[2]),
    calibration_slope_se = unname(sqrt(vcov(slope_fit)[2, 2]))
  )
}

recalibrate_predictions <- function(data, prediction) {
  analysis <- data %>%
    mutate(
      prediction_direct = pmin(pmax(prediction, 1e-6), 1 - 1e-6),
      linear_predictor = qlogis(prediction_direct)
    )
  design <- make_design(analysis)
  intercept_fit <- svyglm(
    model_outcome ~ 1,
    offset = linear_predictor,
    design = design,
    family = quasibinomial()
  )
  slope_fit <- svyglm(
    model_outcome ~ linear_predictor,
    design = design,
    family = quasibinomial()
  )
  intercept <- unname(coef(intercept_fit)[1])
  updated <- coef(slope_fit)
  list(
    Direct = analysis$prediction_direct,
    `Intercept recalibration` = plogis(intercept + analysis$linear_predictor),
    `Intercept-slope recalibration` = plogis(updated[1] + updated[2] * analysis$linear_predictor),
    intercept_update = intercept,
    slope_update_intercept = unname(updated[1]),
    slope_update = unname(updated[2])
  )
}

if (!file.exists(input_file)) stop("Step 1 input is missing: ", input_file)
if (!file.exists(model_file)) stop("Step 3 model object is missing: ", model_file)

imaging <- readRDS(input_file) %>%
  mutate(
    period = case_when(
      analysis_role == "development" ~ "NHANES 2005-2010",
      analysis_role == "temporal_validation_1" ~ "NHANES 2013-2014",
      analysis_role == "temporal_validation_2" ~ "NHANES 2017-Mar 2020",
      analysis_role == "historical_geographic_validation" ~ "KNHANES 2008-2011",
      analysis_role == "contemporary_temporal_validation" ~ "KNHANES 2024",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(period))

period_data <- split(imaging, imaging$period)
models <- readRDS(model_file)$NH_US
source <- period_data[["NHANES 2005-2010"]] %>%
  filter(osteoporosis_fn_us %in% c(0, 1), complete.cases(age, weight_kg, osta_raw, survey_weight, strata, psu))

target_definitions <- tribble(
  ~target_period, ~target_label,
  "KNHANES 2008-2011", "Historical Korea",
  "KNHANES 2024", "Contemporary Korea"
)

decomposition_rows <- list()
balance_rows <- list()
recalibration_rows <- list()
drift_performance_rows <- list()

for (target_index in seq_len(nrow(target_definitions))) {
  target_period <- target_definitions$target_period[target_index]
  target_label <- target_definitions$target_label[target_index]
  target <- period_data[[target_period]] %>%
    filter(
      osteoporosis_fn_us %in% c(0, 1),
      osteoporosis_fn_local %in% c(0, 1),
      complete.cases(age, weight_kg, osta_raw, survey_weight, strata, psu)
    )

  domain_data <- bind_rows(
    source %>% transmute(age, weight_kg, bmi, domain = 0L, base_weight = survey_weight),
    target %>% transmute(age, weight_kg, bmi, domain = 1L, base_weight = survey_weight)
  ) %>%
    group_by(domain) %>%
    mutate(domain_weight = base_weight / sum(base_weight) * n()) %>%
    ungroup()
  domain_fit <- glm(
    domain ~ ns(age, df = 3) + ns(weight_kg, df = 3) + ns(bmi, df = 3) + age:weight_kg,
    data = domain_data,
    family = quasibinomial(),
    weights = domain_weight
  )
  source_domain_probability <- pmin(pmax(
    predict(domain_fit, newdata = source, type = "response"),
    1e-4
  ), 1 - 1e-4)
  raw_density_ratio <- source_domain_probability / (1 - source_domain_probability)
  trim_limits <- weighted_quantile(raw_density_ratio, source$survey_weight, c(0.01, 0.99))
  trimmed_density_ratio <- pmin(pmax(raw_density_ratio, trim_limits[1]), trim_limits[2])
  capped_density_ratio <- pmin(pmax(raw_density_ratio, 0.05), 20)

  source_weighted <- source %>%
    mutate(
      density_ratio_raw = raw_density_ratio,
      density_ratio_trimmed = trimmed_density_ratio,
      density_ratio_capped = capped_density_ratio,
      transport_weight = survey_weight * density_ratio_trimmed,
      transport_weight_capped = survey_weight * density_ratio_capped
    )

  source_prevalence <- prevalence_summary(source, "osteoporosis_fn_us")
  reweighted_prevalence <- prevalence_summary(source_weighted, "osteoporosis_fn_us", "transport_weight")
  reweighted_capped <- prevalence_summary(source_weighted, "osteoporosis_fn_us", "transport_weight_capped")
  target_us_summary <- prevalence_summary(target, "osteoporosis_fn_us")
  target_local_summary <- prevalence_summary(target, "osteoporosis_fn_local")

  decomposition_rows[[target_index]] <- tibble(
    target_period = target_period,
    target_label = target_label,
    source_n = nrow(source),
    target_n = nrow(target),
    source_kish_n = kish_n(source$survey_weight),
    transported_source_kish_n = kish_n(source_weighted$transport_weight),
    transported_source_kish_fraction = kish_n(source_weighted$transport_weight) / kish_n(source$survey_weight),
    density_ratio_trim_lower = trim_limits[1],
    density_ratio_trim_upper = trim_limits[2],
    source_us_prevalence = source_prevalence$prevalence,
    source_us_lower = source_prevalence$prevalence_lower,
    source_us_upper = source_prevalence$prevalence_upper,
    transported_source_us_prevalence = reweighted_prevalence$prevalence,
    transported_source_us_lower = reweighted_prevalence$prevalence_lower,
    transported_source_us_upper = reweighted_prevalence$prevalence_upper,
    target_us_prevalence = target_us_summary$prevalence,
    target_us_lower = target_us_summary$prevalence_lower,
    target_us_upper = target_us_summary$prevalence_upper,
    target_local_prevalence = target_local_summary$prevalence,
    target_local_lower = target_local_summary$prevalence_lower,
    target_local_upper = target_local_summary$prevalence_upper,
    population_composition_component = reweighted_prevalence$prevalence - source_prevalence$prevalence,
    residual_common_reference_component = target_us_summary$prevalence - reweighted_prevalence$prevalence,
    reference_standard_component = target_local_summary$prevalence - target_us_summary$prevalence,
    total_local_difference = target_local_summary$prevalence - source_prevalence$prevalence,
    capped_transport_prevalence = reweighted_capped$prevalence,
    limited_support_warning =
      kish_n(source_weighted$transport_weight) / kish_n(source$survey_weight) < 0.50 |
      trim_limits[2] > 20,
    decomposition_identity_error = abs(
      (reweighted_prevalence$prevalence - source_prevalence$prevalence) +
        (target_us_summary$prevalence - reweighted_prevalence$prevalence) +
        (target_local_summary$prevalence - target_us_summary$prevalence) -
        (target_local_summary$prevalence - source_prevalence$prevalence)
    )
  )

  for (variable in c("age", "weight_kg", "bmi", "osta_raw")) {
    balance_rows[[length(balance_rows) + 1L]] <- tibble(
      target_period = target_period,
      variable = variable,
      source_mean = weighted_mean(source[[variable]], source$survey_weight),
      target_mean = weighted_mean(target[[variable]], target$survey_weight),
      transported_source_mean = weighted_mean(source_weighted[[variable]], source_weighted$transport_weight),
      smd_before = standardized_difference(
        source[[variable]], source$survey_weight, target[[variable]], target$survey_weight
      ),
      smd_after = standardized_difference(
        source_weighted[[variable]], source_weighted$transport_weight, target[[variable]], target$survey_weight
      )
    )
  }

  source_predictions <- predict_selected_models(models, source_weighted)
  for (model_name in names(source_predictions)) {
    direct_prediction <- source_predictions[[model_name]]
    original_metrics <- performance_metrics(
      source_weighted$osteoporosis_fn_us, direct_prediction, source_weighted$survey_weight
    )
    transported_metrics <- performance_metrics(
      source_weighted$osteoporosis_fn_us, direct_prediction, source_weighted$transport_weight
    )
    drift_performance_rows[[length(drift_performance_rows) + 1L]] <- bind_cols(
      tibble(target_period = target_period, model = model_name),
      original_metrics %>% rename_with(~paste0("source_", .x)),
      transported_metrics %>% rename_with(~paste0("transported_source_", .x))
    )
  }

  target_predictions <- predict_selected_models(models, target)
  for (outcome_definition in c("osteoporosis_fn_local", "osteoporosis_fn_us")) {
    outcome_reference <- ifelse(outcome_definition == "osteoporosis_fn_local", "Asian local", "US common")
    target_analysis <- target %>% mutate(model_outcome = .data[[outcome_definition]])
    for (model_name in names(target_predictions)) {
      updated_predictions <- recalibrate_predictions(target_analysis, target_predictions[[model_name]])
      for (update_name in c("Direct", "Intercept recalibration", "Intercept-slope recalibration")) {
        prediction <- updated_predictions[[update_name]]
        metrics <- performance_metrics(target_analysis$model_outcome, prediction, target_analysis$survey_weight)
        calibration <- calibration_statistics(target_analysis, prediction)
        recalibration_rows[[length(recalibration_rows) + 1L]] <- bind_cols(
          tibble(
            target_period = target_period,
            target_label = target_label,
            outcome = outcome_definition,
            outcome_reference = outcome_reference,
            model = model_name,
            updating_method = update_name,
            target_n = nrow(target_analysis),
            target_events = sum(target_analysis$model_outcome == 1),
            intercept_update = updated_predictions$intercept_update,
            slope_update_intercept = updated_predictions$slope_update_intercept,
            slope_update = updated_predictions$slope_update,
            validation_status = ifelse(
              update_name == "Direct",
              "External validation without updating",
              "Apparent target-cohort model updating; not external validation"
            )
          ),
          metrics,
          as_tibble_row(calibration)
        )
      }
    }
  }
}

decomposition <- bind_rows(decomposition_rows)
balance <- bind_rows(balance_rows)
recalibration <- bind_rows(recalibration_rows)
drift_performance <- bind_rows(drift_performance_rows)

decomposition_plot_data <- decomposition %>%
  select(
    target_period,
    `Population composition` = population_composition_component,
    `Residual cross-country difference` = residual_common_reference_component,
    `T-score reference standard` = reference_standard_component
  ) %>%
  pivot_longer(-target_period, names_to = "component", values_to = "risk_difference") %>%
  mutate(component = factor(
    component,
    levels = c("Population composition", "Residual cross-country difference", "T-score reference standard")
  ))

decomposition_plot <- ggplot(
  decomposition_plot_data,
  aes(x = target_period, y = risk_difference * 100, fill = component)
) +
  geom_col(width = 0.62) +
  geom_hline(yintercept = 0, color = "#4B5258", linewidth = 0.4) +
  scale_fill_manual(values = c(
    "Population composition" = "#245E9A",
    "Residual cross-country difference" = "#19726B",
    "T-score reference standard" = "#C07A32"
  )) +
  labs(
    title = "Decomposition of the osteoporosis prevalence difference",
    subtitle = "Additive risk-scale components relative to NHANES 2005-2010",
    x = NULL,
    y = "Absolute prevalence difference (percentage points)",
    fill = NULL
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    plot.title = element_text(face = "bold", color = "#1F2428"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
ggsave(file.path(output_dir, "Figure5_drift_decomposition.png"), decomposition_plot, width = 9, height = 5.8, dpi = 320)

recalibration_plot_data <- recalibration %>%
  filter(model == "OSTA", outcome_reference == "Asian local") %>%
  mutate(updating_method = factor(
    updating_method,
    levels = c("Direct", "Intercept recalibration", "Intercept-slope recalibration")
  ))
recalibration_plot <- ggplot(
  recalibration_plot_data,
  aes(x = updating_method, y = brier, color = target_period, group = target_period)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c("KNHANES 2008-2011" = "#245E9A", "KNHANES 2024" = "#C07A32")) +
  labs(
    title = "Effect of target-cohort recalibration on OSTA accuracy",
    subtitle = "Lower survey-weighted Brier score indicates better absolute-risk accuracy",
    x = NULL,
    y = "Brier score",
    color = NULL
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    plot.title = element_text(face = "bold", color = "#1F2428"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 15, hjust = 1)
  )
ggsave(file.path(output_dir, "Figure6_recalibration.png"), recalibration_plot, width = 9, height = 5.8, dpi = 320)

qc <- bind_rows(
  decomposition %>% transmute(
    check = paste0(target_period, ": additive decomposition identity"),
    passed = decomposition_identity_error < 1e-10,
    value = decomposition_identity_error
  ),
  balance %>% group_by(target_period) %>% summarise(
    check = paste0(first(target_period), ": mean absolute SMD improved"),
    passed = mean(abs(smd_after), na.rm = TRUE) < mean(abs(smd_before), na.rm = TRUE),
    value = mean(abs(smd_after), na.rm = TRUE),
    .groups = "drop"
  ) %>% select(check, passed, value),
  recalibration %>% filter(updating_method == "Direct") %>% transmute(
    check = paste(target_period, outcome_reference, model, "direct metrics finite", sep = ": "),
    passed = is.finite(auc) & is.finite(pr_auc) & is.finite(brier) &
      is.finite(calibration_intercept) & is.finite(calibration_slope),
    value = brier
  ),
  decomposition %>% transmute(
    check = paste0(target_period, ": transport common-support diagnostic"),
    passed = !limited_support_warning,
    value = transported_source_kish_fraction
  )
) %>%
  mutate(severity = if_else(grepl("common-support diagnostic", check), "warning", "error"))

write.csv(decomposition, file.path(output_dir, "Table5_risk_scale_decomposition.csv"), row.names = FALSE, na = "")
write.csv(recalibration, file.path(output_dir, "Table6_recalibration_performance.csv"), row.names = FALSE, na = "")
write.csv(balance, file.path(output_dir, "TableS5_density_weight_balance.csv"), row.names = FALSE, na = "")
write.csv(drift_performance, file.path(output_dir, "TableS6_transported_source_performance.csv"), row.names = FALSE, na = "")
write.csv(qc, file.path(output_dir, "Step4_quality_control.csv"), row.names = FALSE, na = "")
capture.output(sessionInfo(), file = file.path(output_dir, "Step4_sessionInfo.txt"))

cat("\nSTEP 4 COMPLETED\n")
cat("Output:", output_dir, "\n")
print(decomposition %>% select(
  target_period, source_us_prevalence, transported_source_us_prevalence,
  target_us_prevalence, target_local_prevalence,
  population_composition_component, residual_common_reference_component,
  reference_standard_component
))
print(recalibration %>%
  filter(model == "OSTA", outcome_reference == "Asian local") %>%
  select(target_period, updating_method, auc, brier, calibration_intercept, calibration_slope))
