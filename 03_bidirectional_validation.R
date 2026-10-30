# STEP 3: bidirectional temporal and geographic external validation
# Models are developed in prespecified source periods and transported without recalibration.

required_packages <- c("dplyr", "tidyr", "ggplot2", "survey", "splines", "glmnet", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages before running Step 3: ", paste(missing_packages, collapse = ", "))
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
output_dir <- file.path(project_dir, "step3_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_replicates <- as.integer(Sys.getenv("STEP3_BOOTSTRAP_REPS", unset = "200"))
if (!is.finite(bootstrap_replicates) || bootstrap_replicates < 20) bootstrap_replicates <- 200L

weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
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
  vapply(probs, function(probability) x[which(cumulative >= probability)[1]], numeric(1))
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
  group_id <- cumsum(c(TRUE, diff(score) != 0))
  groups <- split(seq_along(score), group_id)
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
  group_id <- cumsum(c(TRUE, diff(score) != 0))
  groups <- split(seq_along(score), group_id)
  total_positive <- sum(w[y == 1])
  true_positive <- 0
  false_positive <- 0
  recall <- 0
  area <- 0
  for (index in groups) {
    true_positive <- true_positive + sum(w[index][y[index] == 1])
    false_positive <- false_positive + sum(w[index][y[index] == 0])
    new_recall <- true_positive / total_positive
    precision <- true_positive / (true_positive + false_positive)
    area <- area + (new_recall - recall) * precision
    recall <- new_recall
  }
  area
}

basic_metrics <- function(y, prediction, w) {
  keep <- y %in% c(0, 1) & is.finite(prediction) & is.finite(w) & w > 0
  y <- y[keep]
  prediction <- pmin(pmax(prediction[keep], 1e-6), 1 - 1e-6)
  w <- w[keep]
  c(
    auc = weighted_auc(y, prediction, w),
    pr_auc = weighted_pr_auc(y, prediction, w),
    brier = weighted_mean((y - prediction)^2, w)
  )
}

make_design <- function(data) {
  data <- data %>%
    mutate(
      survey_strata = factor(strata),
      survey_psu = factor(psu)
    )
  svydesign(
    ids = ~survey_psu,
    strata = ~survey_strata,
    weights = ~survey_weight,
    nest = TRUE,
    data = data
  )
}

ridge_matrix <- function(data, bmi_median) {
  prepared <- data %>%
    mutate(bmi_imputed = if_else(is.finite(bmi), bmi, bmi_median))
  model.matrix(
    ~ age + I(age^2) + weight_kg + I(weight_kg^2) + bmi_imputed + age:weight_kg,
    data = prepared
  )[, -1, drop = FALSE]
}

fit_model_suite <- function(data, outcome_variable, source_key) {
  fit_data <- data %>%
    mutate(model_outcome = .data[[outcome_variable]]) %>%
    filter(
      model_outcome %in% c(0, 1),
      complete.cases(age, weight_kg, osta_raw, survey_weight, strata, psu),
      survey_weight > 0
    )
  design <- make_design(fit_data)
  age_model <- svyglm(model_outcome ~ age, design = design, family = quasibinomial())
  osta_model <- svyglm(model_outcome ~ osta_raw, design = design, family = quasibinomial())
  spline_model <- svyglm(
    model_outcome ~ ns(age, df = 3) + ns(weight_kg, df = 3),
    design = design,
    family = quasibinomial()
  )

  bmi_median <- weighted_quantile(fit_data$bmi, fit_data$survey_weight, 0.5)
  x_ridge <- ridge_matrix(fit_data, bmi_median)
  cluster <- interaction(fit_data$strata, fit_data$psu, drop = TRUE)
  unique_clusters <- unique(cluster)
  cluster_fold <- sample(rep(seq_len(10), length.out = length(unique_clusters)))
  fold_id <- cluster_fold[match(cluster, unique_clusters)]
  ridge_weights <- fit_data$survey_weight / mean(fit_data$survey_weight)
  ridge_model <- cv.glmnet(
    x = x_ridge,
    y = fit_data$model_outcome,
    family = "binomial",
    alpha = 0,
    weights = ridge_weights,
    foldid = fold_id,
    type.measure = "deviance",
    standardize = TRUE,
    nlambda = 100
  )

  list(
    source_key = source_key,
    outcome_variable = outcome_variable,
    fit_n = nrow(fit_data),
    fit_events = sum(fit_data$model_outcome == 1),
    age = age_model,
    osta = osta_model,
    spline = spline_model,
    ridge = ridge_model,
    ridge_bmi_median = bmi_median,
    ridge_lambda = ridge_model$lambda.1se
  )
}

predict_model_suite <- function(models, new_data) {
  ridge_x <- ridge_matrix(new_data, models$ridge_bmi_median)
  list(
    "Age only" = as.numeric(predict(models$age, newdata = new_data, type = "response")),
    "OSTA" = as.numeric(predict(models$osta, newdata = new_data, type = "response")),
    "Age-weight spline" = as.numeric(predict(models$spline, newdata = new_data, type = "response")),
    "Ridge" = as.numeric(predict(models$ridge, newx = ridge_x, s = "lambda.1se", type = "response"))
  )
}

calibration_statistics <- function(data, prediction) {
  analysis_data <- data %>%
    mutate(
      prediction = pmin(pmax(prediction, 1e-6), 1 - 1e-6),
      linear_predictor = qlogis(prediction)
    ) %>%
    filter(model_outcome %in% c(0, 1), complete.cases(prediction, linear_predictor))
  design <- make_design(analysis_data)
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
  intercept_estimate <- coef(intercept_fit)[1]
  intercept_se <- sqrt(vcov(intercept_fit)[1, 1])
  slope_estimate <- coef(slope_fit)[2]
  slope_se <- sqrt(vcov(slope_fit)[2, 2])
  tibble(
    calibration_intercept = unname(intercept_estimate),
    calibration_intercept_se = unname(intercept_se),
    calibration_intercept_lower = unname(intercept_estimate - 1.96 * intercept_se),
    calibration_intercept_upper = unname(intercept_estimate + 1.96 * intercept_se),
    calibration_slope = unname(slope_estimate),
    calibration_slope_se = unname(slope_se),
    calibration_slope_lower = unname(slope_estimate - 1.96 * slope_se),
    calibration_slope_upper = unname(slope_estimate + 1.96 * slope_se)
  )
}

replicate_metric_statistics <- function(data, prediction, replicates) {
  analysis_data <- data %>%
    mutate(prediction = prediction) %>%
    filter(model_outcome %in% c(0, 1), complete.cases(prediction, survey_weight, strata, psu))
  design <- make_design(analysis_data)
  point <- basic_metrics(
    analysis_data$model_outcome,
    analysis_data$prediction,
    analysis_data$survey_weight
  )
  replicate_design <- as.svrepdesign(
    design,
    type = "bootstrap",
    replicates = replicates,
    mse = TRUE
  )
  replicate_estimate <- withReplicates(
    replicate_design,
    function(weights, data) basic_metrics(data$model_outcome, data$prediction, weights)
  )
  standard_error <- SE(replicate_estimate)
  lower <- point - 1.96 * standard_error
  upper <- point + 1.96 * standard_error
  tibble(
    auc = point[["auc"]],
    auc_se = standard_error[["auc"]],
    auc_lower = max(0, lower[["auc"]]),
    auc_upper = min(1, upper[["auc"]]),
    pr_auc = point[["pr_auc"]],
    pr_auc_se = standard_error[["pr_auc"]],
    pr_auc_lower = max(0, lower[["pr_auc"]]),
    pr_auc_upper = min(1, upper[["pr_auc"]]),
    brier = point[["brier"]],
    brier_se = standard_error[["brier"]],
    brier_lower = max(0, lower[["brier"]]),
    brier_upper = min(1, upper[["brier"]])
  )
}

calibration_curve <- function(data, prediction, groups = 10) {
  keep <- data$model_outcome %in% c(0, 1) & is.finite(prediction) &
    is.finite(data$survey_weight) & data$survey_weight > 0
  curve_data <- data[keep, , drop = FALSE]
  curve_data$prediction <- pmin(pmax(prediction[keep], 1e-6), 1 - 1e-6)
  boundaries <- unique(weighted_quantile(
    curve_data$prediction,
    curve_data$survey_weight,
    seq(0, 1, length.out = groups + 1)
  ))
  if (length(boundaries) < 3) return(tibble())
  boundaries[1] <- -Inf
  boundaries[length(boundaries)] <- Inf
  curve_data$group <- cut(curve_data$prediction, breaks = boundaries, include.lowest = TRUE, labels = FALSE)
  curve_data %>%
    group_by(group) %>%
    summarise(
      predicted = weighted_mean(prediction, survey_weight),
      observed = weighted_mean(model_outcome, survey_weight),
      unweighted_n = n(),
      events = sum(model_outcome == 1),
      .groups = "drop"
    )
}

if (!file.exists(input_file)) stop("Step 1 input is missing: ", input_file)

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

source_definitions <- tribble(
  ~source_key, ~source_period, ~source_outcome, ~source_description,
  "NH_US", "NHANES 2005-2010", "osteoporosis_fn_us", "NHANES development, US reference",
  "KN_US", "KNHANES 2008-2011", "osteoporosis_fn_us", "KNHANES historical, common US reference",
  "KN_LOCAL", "KNHANES 2008-2011", "osteoporosis_fn_local", "KNHANES historical, local Asian reference"
)

model_sets <- list()
for (index in seq_len(nrow(source_definitions))) {
  definition <- source_definitions[index, ]
  model_sets[[definition$source_key]] <- fit_model_suite(
    period_data[[definition$source_period]],
    definition$source_outcome,
    definition$source_key
  )
}

scenarios <- tribble(
  ~scenario, ~direction, ~source_key, ~target_period, ~target_outcome, ~outcome_reference, ~role,
  "NH development apparent", "Forward", "NH_US", "NHANES 2005-2010", "osteoporosis_fn_us", "US", "Development apparent",
  "NH temporal validation 1", "Forward", "NH_US", "NHANES 2013-2014", "osteoporosis_fn_us", "US", "Temporal validation",
  "NH temporal validation 2", "Forward", "NH_US", "NHANES 2017-Mar 2020", "osteoporosis_fn_us", "US", "Temporal validation",
  "KN historical geographic", "Forward", "NH_US", "KNHANES 2008-2011", "osteoporosis_fn_local", "Asian local", "Geographic validation",
  "KN contemporary temporal-geographic", "Forward", "NH_US", "KNHANES 2024", "osteoporosis_fn_local", "Asian local", "Temporal-geographic validation",
  "KN historical common reference", "Forward sensitivity", "NH_US", "KNHANES 2008-2011", "osteoporosis_fn_us", "US", "Reference-standard sensitivity",
  "KN contemporary common reference", "Forward sensitivity", "NH_US", "KNHANES 2024", "osteoporosis_fn_us", "US", "Reference-standard sensitivity",
  "Reverse KN-to-NH common reference", "Reverse", "KN_US", "NHANES 2005-2010", "osteoporosis_fn_us", "US", "Reverse geographic validation",
  "Korean temporal validation", "Korean temporal", "KN_LOCAL", "KNHANES 2024", "osteoporosis_fn_local", "Asian local", "Within-country temporal validation"
)

performance_rows <- list()
calibration_rows <- list()
row_id <- 1L
curve_id <- 1L

for (scenario_index in seq_len(nrow(scenarios))) {
  scenario_row <- scenarios[scenario_index, ]
  models <- model_sets[[scenario_row$source_key]]
  target_data <- period_data[[scenario_row$target_period]] %>%
    mutate(model_outcome = .data[[scenario_row$target_outcome]]) %>%
    filter(
      model_outcome %in% c(0, 1),
      complete.cases(age, weight_kg, osta_raw, survey_weight, strata, psu),
      survey_weight > 0
    )
  predictions <- predict_model_suite(models, target_data)

  prevalence_design <- make_design(target_data)
  prevalence_estimate <- svymean(~model_outcome, prevalence_design)
  prevalence <- as.numeric(prevalence_estimate[1])

  for (model_name in names(predictions)) {
    prediction <- predictions[[model_name]]
    replicate_metrics <- replicate_metric_statistics(target_data, prediction, bootstrap_replicates)
    calibration <- calibration_statistics(target_data, prediction)
    performance_rows[[row_id]] <- bind_cols(
      scenario_row,
      tibble(
        model = model_name,
        target_n = nrow(target_data),
        target_events = sum(target_data$model_outcome == 1),
        weighted_prevalence = prevalence,
        bootstrap_replicates = bootstrap_replicates,
        confidence_interval_scope = "Target survey-design uncertainty; transported source model fixed"
      ),
      replicate_metrics,
      calibration
    )
    row_id <- row_id + 1L

    calibration_rows[[curve_id]] <- calibration_curve(target_data, prediction) %>%
      mutate(
        scenario = scenario_row$scenario,
        direction = scenario_row$direction,
        model = model_name,
        outcome_reference = scenario_row$outcome_reference,
        .before = 1
      )
    curve_id <- curve_id + 1L
  }
}

performance <- bind_rows(performance_rows)
calibration_points <- bind_rows(calibration_rows)

coefficient_rows <- list()
coefficient_id <- 1L
for (index in seq_len(nrow(source_definitions))) {
  definition <- source_definitions[index, ]
  models <- model_sets[[definition$source_key]]
  for (model_name in c("age", "osta", "spline")) {
    estimates <- coef(models[[model_name]])
    coefficient_rows[[coefficient_id]] <- tibble(
      source_key = definition$source_key,
      source_period = definition$source_period,
      source_outcome = definition$source_outcome,
      model = recode(model_name, age = "Age only", osta = "OSTA", spline = "Age-weight spline"),
      term = names(estimates),
      estimate = as.numeric(estimates),
      lambda = NA_real_,
      fit_n = models$fit_n,
      fit_events = models$fit_events
    )
    coefficient_id <- coefficient_id + 1L
  }
  ridge_coefficients <- as.matrix(coef(models$ridge, s = "lambda.1se"))[, 1]
  coefficient_rows[[coefficient_id]] <- tibble(
    source_key = definition$source_key,
    source_period = definition$source_period,
    source_outcome = definition$source_outcome,
    model = "Ridge",
    term = names(ridge_coefficients),
    estimate = as.numeric(ridge_coefficients),
    lambda = models$ridge_lambda,
    fit_n = models$fit_n,
    fit_events = models$fit_events
  )
  coefficient_id <- coefficient_id + 1L
}
coefficients <- bind_rows(coefficient_rows)

core_scenarios <- c(
  "NH development apparent", "NH temporal validation 1", "NH temporal validation 2",
  "KN historical geographic", "KN contemporary temporal-geographic"
)
performance_plot_data <- performance %>%
  filter(scenario %in% core_scenarios) %>%
  mutate(
    scenario = factor(scenario, levels = rev(core_scenarios)),
    model = factor(model, levels = c("Age only", "OSTA", "Age-weight spline", "Ridge"))
  )

performance_plot <- ggplot(performance_plot_data, aes(x = auc, y = scenario, color = model)) +
  geom_vline(xintercept = 0.75, color = "#8A9299", linetype = "dashed") +
  geom_errorbar(
    aes(xmin = auc_lower, xmax = auc_upper),
    width = 0.14,
    orientation = "y",
    position = position_dodge(width = 0.56)
  ) +
  geom_point(size = 2.7, position = position_dodge(width = 0.56)) +
  scale_color_manual(values = c(
    "Age only" = "#7A7F85", "OSTA" = "#245E9A",
    "Age-weight spline" = "#19726B", "Ridge" = "#C07A32"
  )) +
  coord_cartesian(xlim = c(0.55, 0.90)) +
  labs(
    title = "Discrimination across development and validation settings",
    subtitle = "Survey-weighted AUC with 95% confidence intervals from bootstrap replicate designs",
    x = "Area under the ROC curve", y = NULL, color = NULL
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    plot.title = element_text(face = "bold", color = "#1F2428"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
ggsave(file.path(output_dir, "Figure3_discrimination.png"), performance_plot, width = 10, height = 6.4, dpi = 320)

calibration_plot_data <- calibration_points %>%
  filter(
    scenario %in% core_scenarios,
    model == "Age-weight spline"
  ) %>%
  mutate(scenario = factor(scenario, levels = core_scenarios))
calibration_plot <- ggplot(calibration_plot_data, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, color = "#70777D", linetype = "dashed") +
  geom_line(color = "#245E9A", linewidth = 0.9) +
  geom_point(aes(size = unweighted_n), color = "#245E9A", fill = "white", shape = 21, stroke = 0.9) +
  facet_wrap(~scenario, ncol = 3) +
  coord_equal(xlim = c(0, 0.6), ylim = c(0, 0.6)) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  labs(
    title = "Calibration of the transported age-weight spline model",
    subtitle = "Survey-weighted observed versus predicted femoral-neck osteoporosis by risk group",
    x = "Mean predicted probability", y = "Observed prevalence"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", color = "#1F2428"),
    panel.grid.minor = element_blank()
  )
ggsave(file.path(output_dir, "Figure4_calibration.png"), calibration_plot, width = 10, height = 7.2, dpi = 320)

qc <- performance %>%
  transmute(
    scenario,
    model,
    target_n,
    target_events,
    prediction_metrics_finite = is.finite(auc) & is.finite(pr_auc) & is.finite(brier),
    auc_in_range = between(auc, 0, 1),
    pr_auc_in_range = between(pr_auc, 0, 1),
    brier_in_range = between(brier, 0, 1),
    calibration_finite = is.finite(calibration_intercept) & is.finite(calibration_slope)
  )

write.csv(performance, file.path(output_dir, "Table4_model_performance.csv"), row.names = FALSE, na = "")
write.csv(coefficients, file.path(output_dir, "TableS3_model_coefficients.csv"), row.names = FALSE, na = "")
write.csv(calibration_points, file.path(output_dir, "TableS4_calibration_points.csv"), row.names = FALSE, na = "")
write.csv(qc, file.path(output_dir, "Step3_quality_control.csv"), row.names = FALSE, na = "")
saveRDS(model_sets, file.path(output_dir, "Step3_fitted_models.rds"))
capture.output(sessionInfo(), file = file.path(output_dir, "Step3_sessionInfo.txt"))

cat("\nSTEP 3 COMPLETED\n")
cat("Output:", output_dir, "\n")
print(performance %>%
  filter(model %in% c("OSTA", "Age-weight spline"), scenario %in% core_scenarios) %>%
  select(scenario, model, auc, auc_lower, auc_upper, brier, calibration_intercept, calibration_slope))
