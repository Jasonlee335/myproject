# STEP 5: capacity-constrained referral analysis
# Primary question: can a portable ranking improve case finding when DXA capacity is limited?

required_packages <- c("dplyr", "tidyr", "ggplot2", "survey", "splines", "glmnet", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages before running Step 5: ", paste(missing_packages, collapse = ", "))
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
output_dir <- file.path(project_dir, "step5_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_replicates <- as.integer(Sys.getenv("STEP5_BOOTSTRAP_REPS", unset = "200"))
if (!is.finite(bootstrap_replicates) || bootstrap_replicates < 20) bootstrap_replicates <- 200L
capacities <- c(0.10, 0.20, 0.30, 0.40)
primary_capacity <- 0.20

weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
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

ridge_matrix <- function(data, bmi_median) {
  prepared <- data %>% mutate(bmi_imputed = if_else(is.finite(bmi), bmi, bmi_median))
  model.matrix(
    ~ age + I(age^2) + weight_kg + I(weight_kg^2) + bmi_imputed + age:weight_kg,
    data = prepared
  )[, -1, drop = FALSE]
}

predict_risk_scores <- function(models, data) {
  ridge_x <- ridge_matrix(data, models$ridge_bmi_median)
  list(
    osta = -data$osta_raw,
    spline = as.numeric(predict(models$spline, newdata = data, type = "response")),
    ridge = as.numeric(predict(models$ridge, newx = ridge_x, s = "lambda.1se", type = "response"))
  )
}

capacity_rule <- function(score, w, capacity) {
  keep <- is.finite(score) & is.finite(w) & w > 0
  if (!any(keep)) {
    return(list(cutoff = NA_real_, tie_fraction = NA_real_, referral = rep(NA_real_, length(score))))
  }
  score_valid <- score[keep]
  weight_valid <- w[keep]
  ord <- order(score_valid, decreasing = TRUE)
  score_ordered <- score_valid[ord]
  weight_ordered <- weight_valid[ord]
  target_weight <- capacity * sum(weight_ordered)
  cumulative <- cumsum(weight_ordered)
  cutoff_index <- which(cumulative >= target_weight)[1]
  cutoff <- score_ordered[cutoff_index]
  weight_above <- sum(weight_valid[score_valid > cutoff])
  weight_at <- sum(weight_valid[score_valid == cutoff])
  tie_fraction <- if (weight_at > 0) {
    min(1, max(0, (target_weight - weight_above) / weight_at))
  } else {
    0
  }
  referral <- rep(NA_real_, length(score))
  referral[keep] <- ifelse(score_valid > cutoff, 1, ifelse(score_valid == cutoff, tie_fraction, 0))
  list(cutoff = cutoff, tie_fraction = tie_fraction, referral = referral)
}

apply_rule <- function(score, cutoff, tie_fraction) {
  if (!is.finite(cutoff) || !is.finite(tie_fraction)) return(rep(NA_real_, length(score)))
  ifelse(!is.finite(score), NA_real_, ifelse(score > cutoff, 1, ifelse(score == cutoff, tie_fraction, 0)))
}

screening_metrics <- function(y, referral, w) {
  keep <- y %in% c(0, 1) & is.finite(referral) & is.finite(w) & w > 0
  y <- y[keep]
  referral <- pmin(pmax(referral[keep], 0), 1)
  w <- w[keep]
  tp <- sum(w * referral * y)
  fp <- sum(w * referral * (1 - y))
  fn <- sum(w * (1 - referral) * y)
  tn <- sum(w * (1 - referral) * (1 - y))
  referred <- tp + fp
  total <- tp + fp + fn + tn
  cases <- tp + fn
  noncases <- fp + tn
  c(
    referral_rate = referred / total,
    sensitivity = ifelse(cases > 0, tp / cases, NA_real_),
    specificity = ifelse(noncases > 0, tn / noncases, NA_real_),
    ppv = ifelse(referred > 0, tp / referred, NA_real_),
    npv = ifelse((tn + fn) > 0, tn / (tn + fn), NA_real_),
    number_needed_to_screen = ifelse(tp > 0, referred / tp, NA_real_),
    missed_case_rate = ifelse(cases > 0, fn / cases, NA_real_)
  )
}

strategy_allocation <- function(strategy, capacity, target_scores, target_data, target_w,
                                source_scores, source_data, source_w) {
  if (strategy == "Clinical OSTA <= -1") {
    return(list(
      referral = as.numeric(target_data$osta <= -1),
      cutoff = 1,
      tie_fraction = 1,
      threshold_origin = "Published clinical OSTA cutoff"
    ))
  }

  score_key <- case_when(
    grepl("OSTA", strategy, fixed = TRUE) ~ "osta",
    grepl("spline", strategy, fixed = TRUE) ~ "spline",
    grepl("ridge", strategy, fixed = TRUE) ~ "ridge",
    TRUE ~ NA_character_
  )
  if (is.na(score_key)) stop("Unknown strategy: ", strategy)

  if (grepl("Source-fixed", strategy, fixed = TRUE)) {
    source_rule <- capacity_rule(source_scores[[score_key]], source_w, capacity)
    return(list(
      referral = apply_rule(target_scores[[score_key]], source_rule$cutoff, source_rule$tie_fraction),
      cutoff = source_rule$cutoff,
      tie_fraction = source_rule$tie_fraction,
      threshold_origin = "NHANES 2005-2010 capacity threshold"
    ))
  }

  target_rule <- capacity_rule(target_scores[[score_key]], target_w, capacity)
  list(
    referral = target_rule$referral,
    cutoff = target_rule$cutoff,
    tie_fraction = target_rule$tie_fraction,
    threshold_origin = "Unlabelled within-target weighted percentile"
  )
}

bootstrap_interval <- function(values) {
  valid <- values[is.finite(values)]
  if (length(valid) < 10) return(c(lower = NA_real_, upper = NA_real_))
  unname(quantile(valid, c(0.025, 0.975), names = FALSE, type = 6))
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
  filter(
    !is.na(period),
    complete.cases(age, weight_kg, osta_raw, osta, survey_weight, strata, psu),
    survey_weight > 0
  )

period_data <- split(imaging, imaging$period)
models <- readRDS(model_file)$NH_US
source_period <- "NHANES 2005-2010"
source_data <- period_data[[source_period]]
source_scores <- predict_risk_scores(models, source_data)

replicate_weights <- lapply(period_data, make_replicate_weights, replicates = bootstrap_replicates)

scenario_definitions <- tribble(
  ~target_period, ~outcome, ~outcome_reference, ~outcome_scope,
  "NHANES 2005-2010", "osteoporosis_fn_us", "US", "Primary definition",
  "NHANES 2013-2014", "osteoporosis_fn_us", "US", "Primary definition",
  "NHANES 2017-Mar 2020", "osteoporosis_fn_us", "US", "Primary definition",
  "KNHANES 2008-2011", "osteoporosis_fn_local", "Asian local", "Primary definition",
  "KNHANES 2024", "osteoporosis_fn_local", "Asian local", "Primary definition",
  "KNHANES 2008-2011", "osteoporosis_fn_us", "US common", "Reference-standard sensitivity",
  "KNHANES 2024", "osteoporosis_fn_us", "US common", "Reference-standard sensitivity"
)

strategies <- c(
  "Clinical OSTA <= -1",
  "Source-fixed OSTA",
  "Source-fixed spline",
  "Target-percentile OSTA",
  "Target-percentile spline",
  "Target-percentile ridge"
)

performance_rows <- list()
threshold_rows <- list()
subgroup_rows <- list()

for (scenario_index in seq_len(nrow(scenario_definitions))) {
  scenario <- scenario_definitions[scenario_index, ]
  target_data <- period_data[[scenario$target_period]]
  target_scores <- predict_risk_scores(models, target_data)
  y <- target_data[[scenario$outcome]]
  target_point_w <- target_data$survey_weight
  target_rep_w <- replicate_weights[[scenario$target_period]]
  source_rep_w <- replicate_weights[[source_period]]

  for (capacity in capacities) {
    for (strategy in strategies) {
      point_allocation <- strategy_allocation(
        strategy, capacity, target_scores, target_data, target_point_w,
        source_scores, source_data, source_data$survey_weight
      )
      point_metrics <- screening_metrics(y, point_allocation$referral, target_point_w)

      replicate_metrics <- matrix(
        NA_real_,
        nrow = bootstrap_replicates,
        ncol = length(point_metrics),
        dimnames = list(NULL, names(point_metrics))
      )
      for (replicate_index in seq_len(bootstrap_replicates)) {
        replicate_allocation <- strategy_allocation(
          strategy,
          capacity,
          target_scores,
          target_data,
          target_rep_w[, replicate_index],
          source_scores,
          source_data,
          source_rep_w[, replicate_index]
        )
        replicate_metrics[replicate_index, ] <- screening_metrics(
          y,
          replicate_allocation$referral,
          target_rep_w[, replicate_index]
        )
      }

      sensitivity_ci <- bootstrap_interval(replicate_metrics[, "sensitivity"])
      referral_ci <- bootstrap_interval(replicate_metrics[, "referral_rate"])
      ppv_ci <- bootstrap_interval(replicate_metrics[, "ppv"])
      nns_ci <- bootstrap_interval(replicate_metrics[, "number_needed_to_screen"])

      performance_rows[[length(performance_rows) + 1L]] <- bind_cols(
        scenario,
        tibble(
          strategy = strategy,
          intended_capacity = capacity,
          target_n = nrow(target_data),
          target_events = sum(y == 1, na.rm = TRUE),
          cutoff = point_allocation$cutoff,
          tie_fraction = point_allocation$tie_fraction,
          threshold_origin = point_allocation$threshold_origin,
          bootstrap_replicates = bootstrap_replicates,
          confidence_interval_scope = paste0(
            "Survey bootstrap with capacity threshold re-selected; ",
            "NHANES model coefficients fixed"
          ),
          referral_rate = point_metrics[["referral_rate"]],
          referral_rate_lower = referral_ci[1],
          referral_rate_upper = referral_ci[2],
          capacity_deviation = point_metrics[["referral_rate"]] - capacity,
          sensitivity = point_metrics[["sensitivity"]],
          sensitivity_lower = sensitivity_ci[1],
          sensitivity_upper = sensitivity_ci[2],
          specificity = point_metrics[["specificity"]],
          ppv = point_metrics[["ppv"]],
          ppv_lower = ppv_ci[1],
          ppv_upper = ppv_ci[2],
          npv = point_metrics[["npv"]],
          number_needed_to_screen = point_metrics[["number_needed_to_screen"]],
          number_needed_to_screen_lower = nns_ci[1],
          number_needed_to_screen_upper = nns_ci[2],
          missed_case_rate = point_metrics[["missed_case_rate"]]
        )
      )

      threshold_rows[[length(threshold_rows) + 1L]] <- tibble(
        target_period = scenario$target_period,
        outcome_reference = scenario$outcome_reference,
        strategy = strategy,
        intended_capacity = capacity,
        cutoff = point_allocation$cutoff,
        tie_fraction = point_allocation$tie_fraction,
        threshold_origin = point_allocation$threshold_origin
      )

      if (scenario$outcome_scope == "Primary definition" && capacity == primary_capacity) {
        subgroup_definitions <- list(
          `Age 50-64` = target_data$age < 65,
          `Age 65-79` = target_data$age >= 65,
          `BMI <25 kg/m2` = target_data$bmi < 25,
          `BMI >=25 kg/m2` = target_data$bmi >= 25
        )
        for (subgroup_name in names(subgroup_definitions)) {
          subgroup_keep <- subgroup_definitions[[subgroup_name]] & !is.na(subgroup_definitions[[subgroup_name]])
          subgroup_metrics <- screening_metrics(
            y[subgroup_keep],
            point_allocation$referral[subgroup_keep],
            target_point_w[subgroup_keep]
          )
          subgroup_rows[[length(subgroup_rows) + 1L]] <- tibble(
            target_period = scenario$target_period,
            outcome_reference = scenario$outcome_reference,
            subgroup = subgroup_name,
            subgroup_n = sum(subgroup_keep),
            subgroup_events = sum(y[subgroup_keep] == 1, na.rm = TRUE),
            strategy = strategy,
            intended_capacity = capacity,
            referral_rate = subgroup_metrics[["referral_rate"]],
            sensitivity = subgroup_metrics[["sensitivity"]],
            specificity = subgroup_metrics[["specificity"]],
            ppv = subgroup_metrics[["ppv"]],
            number_needed_to_screen = subgroup_metrics[["number_needed_to_screen"]]
          )
        }
      }
    }
  }
}

performance <- bind_rows(performance_rows) %>%
  group_by(target_period, outcome_reference, outcome_scope, intended_capacity) %>%
  mutate(
    best_capacity_constrained_sensitivity = max(
      sensitivity[grepl("Target-percentile", strategy, fixed = TRUE)],
      na.rm = TRUE
    ),
    sensitivity_regret = if_else(
      !grepl("Target-percentile", strategy, fixed = TRUE),
      NA_real_,
      best_capacity_constrained_sensitivity - sensitivity
    ),
    strategy_role = case_when(
      strategy == "Clinical OSTA <= -1" ~ "Published clinical benchmark",
      grepl("Source-fixed", strategy, fixed = TRUE) ~ "Fixed-threshold transportability diagnostic",
      TRUE ~ "Capacity-compliant target-rank implementation"
    )
  ) %>%
  ungroup()

thresholds <- bind_rows(threshold_rows) %>% distinct()
subgroups <- bind_rows(subgroup_rows)

worst_case <- performance %>%
  filter(
    outcome_scope == "Primary definition",
    strategy_role == "Capacity-compliant target-rank implementation"
  ) %>%
  group_by(strategy, intended_capacity) %>%
  summarise(
    minimum_sensitivity = min(sensitivity, na.rm = TRUE),
    median_sensitivity = median(sensitivity, na.rm = TRUE),
    maximum_sensitivity_regret = max(sensitivity_regret, na.rm = TRUE),
    maximum_absolute_capacity_deviation = max(abs(capacity_deviation), na.rm = TRUE),
    median_ppv = median(ppv, na.rm = TRUE),
    median_number_needed_to_screen = median(number_needed_to_screen, na.rm = TRUE),
    .groups = "drop"
  )

fixed_threshold_drift <- performance %>%
  filter(
    outcome_scope == "Primary definition",
    strategy_role != "Capacity-compliant target-rank implementation"
  ) %>%
  select(
    target_period, strategy, strategy_role, intended_capacity,
    referral_rate, capacity_deviation, sensitivity, ppv, number_needed_to_screen
  )

primary_plot_data <- performance %>%
  filter(outcome_scope == "Primary definition", intended_capacity == primary_capacity) %>%
  mutate(
    target_period = factor(target_period, levels = rev(c(
      "NHANES 2005-2010", "NHANES 2013-2014", "NHANES 2017-Mar 2020",
      "KNHANES 2008-2011", "KNHANES 2024"
    ))),
    strategy = factor(strategy, levels = strategies)
  )

capacity_plot <- ggplot(primary_plot_data, aes(x = sensitivity, y = target_period, color = strategy)) +
  geom_errorbar(
    aes(xmin = sensitivity_lower, xmax = sensitivity_upper),
    width = 0.12,
    orientation = "y",
    position = position_dodge(width = 0.62)
  ) +
  geom_point(size = 2.4, position = position_dodge(width = 0.62)) +
  scale_color_manual(values = c(
    "Clinical OSTA <= -1" = "#7A7F85",
    "Source-fixed OSTA" = "#245E9A",
    "Source-fixed spline" = "#19726B",
    "Target-percentile OSTA" = "#9A3E50",
    "Target-percentile spline" = "#C07A32",
    "Target-percentile ridge" = "#5F5A93"
  )) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(
    title = "Case detection under a 20% intended DXA referral capacity",
    subtitle = "Survey-weighted sensitivity; local T-score definitions are primary in each country",
    x = "Sensitivity",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", color = "#1F2428"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
ggsave(file.path(output_dir, "Figure7_capacity_performance.png"), capacity_plot, width = 10.5, height = 7.0, dpi = 320)

tradeoff_plot <- ggplot(
  worst_case,
  aes(x = intended_capacity * 100, y = minimum_sensitivity, color = strategy)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c(
    "Target-percentile OSTA" = "#9A3E50",
    "Target-percentile spline" = "#C07A32",
    "Target-percentile ridge" = "#5F5A93"
  )) +
  labs(
    title = "Worst-case sensitivity at fixed DXA referral capacities",
    subtitle = "Minimum survey-weighted sensitivity across five survey periods",
    x = "Intended referral capacity (%)",
    y = "Minimum sensitivity across five survey periods",
    color = NULL
  ) +
  theme_minimal(base_size = 11.5) +
  theme(
    plot.title = element_text(face = "bold", color = "#1F2428"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
ggsave(file.path(output_dir, "Figure8_worst_case_tradeoff.png"), tradeoff_plot, width = 9, height = 6.2, dpi = 320)

qc <- bind_rows(
  performance %>% filter(strategy != "Clinical OSTA <= -1") %>% transmute(
    check = paste(target_period, outcome_reference, strategy, intended_capacity, "finite metrics", sep = ": "),
    passed = is.finite(referral_rate) & is.finite(sensitivity) & is.finite(specificity) &
      is.finite(ppv) & is.finite(npv) & is.finite(number_needed_to_screen),
    value = sensitivity
  ),
  performance %>% filter(grepl("Target-percentile", strategy, fixed = TRUE)) %>% transmute(
    check = paste(target_period, outcome_reference, strategy, intended_capacity, "capacity exact", sep = ": "),
    passed = abs(capacity_deviation) < 1e-10,
    value = capacity_deviation
  )
)

write.csv(performance, file.path(output_dir, "Table7_capacity_strategy_performance.csv"), row.names = FALSE, na = "")
write.csv(worst_case, file.path(output_dir, "Table8_worst_case_summary.csv"), row.names = FALSE, na = "")
write.csv(fixed_threshold_drift, file.path(output_dir, "Table9_fixed_threshold_capacity_drift.csv"), row.names = FALSE, na = "")
write.csv(subgroups, file.path(output_dir, "TableS7_subgroup_performance.csv"), row.names = FALSE, na = "")
write.csv(thresholds, file.path(output_dir, "TableS8_strategy_thresholds.csv"), row.names = FALSE, na = "")
write.csv(qc, file.path(output_dir, "Step5_quality_control.csv"), row.names = FALSE, na = "")
capture.output(sessionInfo(), file = file.path(output_dir, "Step5_sessionInfo.txt"))

cat("\nSTEP 5 COMPLETED\n")
cat("Output:", output_dir, "\n")
print(performance %>%
  filter(outcome_scope == "Primary definition", intended_capacity == primary_capacity) %>%
  select(target_period, strategy, referral_rate, sensitivity, sensitivity_lower, sensitivity_upper, ppv))
cat("\nWorst-case summary at 20% intended capacity:\n")
print(worst_case %>% filter(intended_capacity == primary_capacity))
