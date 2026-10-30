# STEP 2: survey-weighted description and distribution-shift diagnostics
# Inputs are frozen Step 1 RDS files. No eligibility or outcome definitions are changed here.

required_packages <- c("dplyr", "tidyr", "ggplot2", "survey", "splines", "tibble")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages before running Step 2: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survey)
  library(splines)
  library(tibble)
})

options(survey.lonely.psu = "adjust")
set.seed(20260805)

project_dir <- Sys.getenv("DXA_PROJECT_DIR", unset = file.path(getwd(), "DXA_transportability_results"))
input_file <- file.path(project_dir, "data_clean", "imaging_cross_country_analysis_ready.rds")
output_dir <- file.path(project_dir, "step2_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Step 1 input is missing: ", input_file)

weighted_quantile <- function(x, w, probs = c(0.01, 0.99)) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (length(x) == 0) return(rep(NA_real_, length(probs)))
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

weighted_mean <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  if (!any(keep)) return(NA_real_)
  sum(x[keep] * w[keep]) / sum(w[keep])
}

weighted_var <- function(x, w) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  if (length(x) < 2) return(NA_real_)
  mu <- sum(w * x) / sum(w)
  sum(w * (x - mu)^2) / (sum(w) - sum(w^2) / sum(w))
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
  groups <- split(seq_along(score), score)
  neg_before <- 0
  concordant <- 0
  for (idx in groups) {
    pos_w <- sum(w[idx][y[idx] == 1])
    neg_w <- sum(w[idx][y[idx] == 0])
    concordant <- concordant + pos_w * (neg_before + 0.5 * neg_w)
    neg_before <- neg_before + neg_w
  }
  concordant / (sum(w[y == 1]) * sum(w[y == 0]))
}

weighted_overlap_1d <- function(x1, w1, x2, w2, n_grid = 512) {
  keep1 <- is.finite(x1) & is.finite(w1) & w1 > 0
  keep2 <- is.finite(x2) & is.finite(w2) & w2 > 0
  x1 <- x1[keep1]
  w1 <- w1[keep1]
  x2 <- x2[keep2]
  w2 <- w2[keep2]
  limits <- range(c(x1, x2))
  common_bandwidth <- bw.nrd0(c(x1, x2))
  d1 <- density(
    x1, weights = w1 / sum(w1), bw = common_bandwidth,
    from = limits[1], to = limits[2], n = n_grid
  )
  d2 <- density(
    x2, weights = w2 / sum(w2), bw = common_bandwidth,
    from = limits[1], to = limits[2], n = n_grid
  )
  sum(pmin(d1$y, d2$y)) * (d1$x[2] - d1$x[1])
}

weighted_joint_overlap <- function(source, target) {
  age_breaks <- seq(50, 80, by = 2)
  weight_breaks <- c(seq(25, 120, by = 5), Inf)
  cell_prob <- function(data) {
    keep <- is.finite(data$age) & is.finite(data$weight_kg) &
      is.finite(data$survey_weight) & data$survey_weight > 0
    age_bin <- cut(data$age[keep], age_breaks, right = FALSE, include.lowest = TRUE)
    weight_bin <- cut(data$weight_kg[keep], weight_breaks, right = FALSE, include.lowest = TRUE)
    cells <- interaction(age_bin, weight_bin, drop = FALSE)
    totals <- tapply(data$survey_weight[keep], cells, sum, default = 0)
    as.numeric(totals / sum(totals))
  }
  sum(pmin(cell_prob(source), cell_prob(target)))
}

imaging <- readRDS(input_file) %>%
  mutate(
    period = case_when(
      analysis_role == "development" ~ "NHANES 2005-2010",
      analysis_role == "temporal_validation_1" ~ "NHANES 2013-2014",
      analysis_role == "temporal_validation_2" ~ "NHANES 2017-Mar 2020",
      analysis_role == "historical_geographic_validation" ~ "KNHANES 2008-2011",
      analysis_role == "contemporary_temporal_validation" ~ "KNHANES 2024",
      TRUE ~ NA_character_
    ),
    primary_fn_op = if_else(cohort == "KNHANES", osteoporosis_fn_local, osteoporosis_fn_us),
    period = factor(
      period,
      levels = c(
        "NHANES 2005-2010", "NHANES 2013-2014", "NHANES 2017-Mar 2020",
        "KNHANES 2008-2011", "KNHANES 2024"
      )
    )
  ) %>%
  filter(!is.na(period))

period_levels <- levels(imaging$period)
period_data <- setNames(lapply(period_levels, function(label) {
  imaging %>% filter(period == label) %>% droplevels()
}), period_levels)

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
period_designs <- lapply(period_data, make_design)

continuous_variables <- c(
  age = "Age, years",
  weight_kg = "Weight, kg",
  bmi = "BMI, kg/m2",
  osta_raw = "OSTA, continuous",
  fn_bmd = "Femoral-neck BMD, g/cm2",
  total_hip_bmd = "Total-hip BMD, g/cm2",
  lumbar_spine_bmd = "Lumbar-spine BMD, g/cm2"
)
binary_variables <- c(
  primary_fn_op = "Femoral-neck osteoporosis, primary reference",
  osteoporosis_fn_us = "Femoral-neck osteoporosis, US reference",
  osteoporosis_any_complete_us = "Any-site osteoporosis, US reference, complete three-site"
)

survey_continuous <- function(design, variable, label, period) {
  estimate <- svymean(as.formula(paste0("~", variable)), design, na.rm = TRUE)
  variance <- svyvar(as.formula(paste0("~", variable)), design, na.rm = TRUE)
  standard_error <- sqrt(as.numeric(attr(estimate, "var")[1, 1]))
  point_estimate <- as.numeric(estimate[1])
  tibble(
    period = period,
    variable = variable,
    label = label,
    type = "continuous",
    estimate = point_estimate,
    standard_error = standard_error,
    ci_lower = point_estimate - 1.96 * standard_error,
    ci_upper = point_estimate + 1.96 * standard_error,
    weighted_sd = sqrt(as.numeric(variance[1]))
  )
}

survey_binary <- function(design, variable, label, period) {
  estimate <- svymean(as.formula(paste0("~", variable)), design, na.rm = TRUE)
  standard_error <- sqrt(as.numeric(attr(estimate, "var")[1, 1]))
  point_estimate <- as.numeric(estimate[1])
  tibble(
    period = period,
    variable = variable,
    label = label,
    type = "binary",
    estimate = point_estimate,
    standard_error = standard_error,
    ci_lower = max(0, point_estimate - 1.96 * standard_error),
    ci_upper = min(1, point_estimate + 1.96 * standard_error),
    weighted_sd = NA_real_
  )
}

descriptive_rows <- list()
row_id <- 1L
for (period in period_levels) {
  design <- period_designs[[period]]
  for (variable in names(continuous_variables)) {
    descriptive_rows[[row_id]] <- survey_continuous(
      design, variable, continuous_variables[[variable]], period
    )
    row_id <- row_id + 1L
  }
  for (variable in names(binary_variables)) {
    descriptive_rows[[row_id]] <- survey_binary(
      design, variable, binary_variables[[variable]], period
    )
    row_id <- row_id + 1L
  }
}
weighted_descriptives <- bind_rows(descriptive_rows)

sample_summary <- bind_rows(lapply(period_levels, function(period) {
  data <- period_data[[period]]
  w <- data$survey_weight
  tibble(
    period = period,
    unweighted_n = nrow(data),
    primary_fn_events = sum(data$primary_fn_op == 1, na.rm = TRUE),
    us_reference_fn_events = sum(data$osteoporosis_fn_us == 1, na.rm = TRUE),
    weighted_population = sum(w, na.rm = TRUE),
    kish_effective_n = sum(w, na.rm = TRUE)^2 / sum(w^2, na.rm = TRUE),
    strata_n = n_distinct(data$strata),
    psu_n = n_distinct(interaction(data$strata, data$psu, drop = TRUE))
  )
}))

source_label <- "NHANES 2005-2010"
source_data <- period_data[[source_label]]
shift_variables <- c("age", "weight_kg", "bmi", "osta_raw", "fn_bmd")

smd_rows <- list()
shift_rows <- list()
for (target_label in setdiff(period_levels, source_label)) {
  target_data <- period_data[[target_label]]

  for (variable in shift_variables) {
    m_source <- weighted_mean(source_data[[variable]], source_data$survey_weight)
    m_target <- weighted_mean(target_data[[variable]], target_data$survey_weight)
    sd_source <- sqrt(weighted_var(source_data[[variable]], source_data$survey_weight))
    sd_target <- sqrt(weighted_var(target_data[[variable]], target_data$survey_weight))
    smd_rows[[length(smd_rows) + 1L]] <- tibble(
      source = source_label,
      target = target_label,
      variable = variable,
      source_mean = m_source,
      target_mean = m_target,
      source_sd = sd_source,
      target_sd = sd_target,
      standardized_mean_difference = (m_target - m_source) / sqrt((sd_source^2 + sd_target^2) / 2)
    )
  }

  domain_data <- bind_rows(
    source_data %>% mutate(domain = 0L),
    target_data %>% mutate(domain = 1L)
  ) %>%
    filter(complete.cases(age, weight_kg, survey_weight)) %>%
    group_by(domain) %>%
    mutate(domain_weight = survey_weight / sum(survey_weight)) %>%
    ungroup()

  domain_model <- glm(
    domain ~ ns(age, df = 3) + ns(weight_kg, df = 3) + age:weight_kg,
    data = domain_data,
    weights = domain_weight,
    family = quasibinomial()
  )
  domain_data$domain_probability <- predict(domain_model, type = "response")
  source_ps <- domain_data$domain_probability[domain_data$domain == 0]
  source_ps_w <- domain_data$domain_weight[domain_data$domain == 0]
  target_ps <- domain_data$domain_probability[domain_data$domain == 1]
  target_ps_w <- domain_data$domain_weight[domain_data$domain == 1]
  source_limits <- weighted_quantile(source_ps, source_ps_w, c(0.01, 0.99))
  target_limits <- weighted_quantile(target_ps, target_ps_w, c(0.01, 0.99))
  common_lower <- max(source_limits[1], target_limits[1])
  common_upper <- min(source_limits[2], target_limits[2])

  age_limits <- weighted_quantile(source_data$age, source_data$survey_weight, c(0.01, 0.99))
  weight_limits <- weighted_quantile(source_data$weight_kg, source_data$survey_weight, c(0.01, 0.99))
  target_rect <- target_data$age >= age_limits[1] & target_data$age <= age_limits[2] &
    target_data$weight_kg >= weight_limits[1] & target_data$weight_kg <= weight_limits[2]

  shift_rows[[length(shift_rows) + 1L]] <- tibble(
    source = source_label,
    target = target_label,
    osta_overlap_coefficient = weighted_overlap_1d(
      source_data$osta_raw, source_data$survey_weight,
      target_data$osta_raw, target_data$survey_weight
    ),
    age_weight_joint_overlap = weighted_joint_overlap(source_data, target_data),
    domain_auc = weighted_auc(
      domain_data$domain, domain_data$domain_probability, domain_data$domain_weight
    ),
    propensity_common_lower = common_lower,
    propensity_common_upper = common_upper,
    source_in_propensity_common_support = weighted_mean(
      as.numeric(source_ps >= common_lower & source_ps <= common_upper), source_ps_w
    ),
    target_in_propensity_common_support = weighted_mean(
      as.numeric(target_ps >= common_lower & target_ps <= common_upper), target_ps_w
    ),
    target_in_source_age_weight_1_99_rectangle = weighted_mean(
      as.numeric(target_rect), target_data$survey_weight
    )
  )
}

pairwise_smd <- bind_rows(smd_rows)
distribution_shift <- bind_rows(shift_rows)

plot_data <- imaging %>%
  group_by(period) %>%
  mutate(plot_weight = survey_weight / sum(survey_weight, na.rm = TRUE)) %>%
  ungroup()

palette <- c(
  "NHANES 2005-2010" = "#245E9A",
  "NHANES 2013-2014" = "#6B93B8",
  "NHANES 2017-Mar 2020" = "#A5B8C8",
  "KNHANES 2008-2011" = "#19726B",
  "KNHANES 2024" = "#C07A32"
)

osta_plot <- ggplot(plot_data, aes(x = osta_raw, weight = plot_weight, color = period)) +
  geom_density(linewidth = 1.05, adjust = 1.05, na.rm = TRUE) +
  scale_color_manual(values = palette) +
  coord_cartesian(xlim = c(-12, 10)) +
  labs(
    title = "Survey-weighted OSTA distributions",
    subtitle = "Women aged 50-79 years; weights normalized within each analytical period",
    x = "OSTA, continuous", y = "Density", color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", color = "#1F2428"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
ggsave(file.path(output_dir, "Figure2A_OSTA_distributions.png"), osta_plot, width = 9.5, height = 5.8, dpi = 320)

smd_plot_data <- pairwise_smd %>%
  mutate(
    variable = recode(
      variable,
      age = "Age", weight_kg = "Weight", bmi = "BMI", osta_raw = "OSTA", fn_bmd = "FN BMD"
    ),
    target = factor(target, levels = rev(setdiff(period_levels, source_label)))
  )
smd_plot <- ggplot(smd_plot_data, aes(x = target, y = variable, fill = standardized_mean_difference)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", standardized_mean_difference)), size = 3.4) +
  scale_fill_gradient2(low = "#245E9A", mid = "#F5F6F7", high = "#C07A32", midpoint = 0, limits = c(-1.5, 1.5), oob = scales::squish) +
  labs(
    title = "Distribution shift relative to NHANES 2005-2010",
    subtitle = "Survey-weighted standardized mean differences; absolute values above 0.10 indicate meaningful imbalance",
    x = NULL, y = NULL, fill = "SMD"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", color = "#1F2428"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 24, hjust = 1),
    legend.position = "right"
  )
ggsave(file.path(output_dir, "Figure2B_standardized_differences.png"), smd_plot, width = 9.5, height = 5.8, dpi = 320)

age_weight_bins <- plot_data %>%
  filter(age >= 50, age < 80, weight_kg >= 25, weight_kg < 120) %>%
  mutate(
    age_bin = 50 + 2 * floor((age - 50) / 2),
    weight_bin = 25 + 5 * floor((weight_kg - 25) / 5)
  ) %>%
  group_by(period, age_bin, weight_bin) %>%
  summarise(weighted_mass = sum(plot_weight, na.rm = TRUE), .groups = "drop")

age_weight_plot <- ggplot(age_weight_bins, aes(x = age_bin + 1, y = weight_bin + 2.5, fill = weighted_mass)) +
  geom_tile(width = 2, height = 5) +
  facet_wrap(~period, ncol = 3) +
  scale_fill_gradient(low = "#F2F4F7", high = "#245E9A", trans = "sqrt") +
  coord_cartesian(xlim = c(50, 79), ylim = c(30, 110)) +
  labs(
    title = "Survey-weighted joint distribution of age and body weight",
    subtitle = "Color represents normalized weighted mass within each period",
    x = "Age, years", y = "Weight, kg", fill = "Weighted mass"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", color = "#1F2428"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
ggsave(file.path(output_dir, "Figure2C_age_weight_joint_distribution.png"), age_weight_plot, width = 10.5, height = 6.8, dpi = 320)

qc <- bind_rows(lapply(period_levels, function(period) {
  data <- period_data[[period]]
  summary_row <- sample_summary %>% filter(.data$period == .env$period)
  tibble(
    period = period,
    unweighted_n = nrow(data),
    positive_weights = all(data$survey_weight > 0),
    complete_design = all(complete.cases(data[, c("survey_weight", "strata", "psu")])),
    both_outcome_classes = n_distinct(data$primary_fn_op) == 2,
    strata_n = summary_row$strata_n,
    psu_n = summary_row$psu_n
  )
}))

write.csv(weighted_descriptives, file.path(output_dir, "Table1_weighted_descriptives.csv"), row.names = FALSE, na = "")
write.csv(sample_summary, file.path(output_dir, "TableS2_period_sample_summary.csv"), row.names = FALSE, na = "")
write.csv(pairwise_smd, file.path(output_dir, "Table2_pairwise_standardized_differences.csv"), row.names = FALSE, na = "")
write.csv(distribution_shift, file.path(output_dir, "Table3_distribution_shift_diagnostics.csv"), row.names = FALSE, na = "")
write.csv(qc, file.path(output_dir, "Step2_quality_control.csv"), row.names = FALSE, na = "")
capture.output(sessionInfo(), file = file.path(output_dir, "Step2_sessionInfo.txt"))

cat("\nSTEP 2 COMPLETED\n")
cat("Output:", output_dir, "\n")
print(sample_summary %>% select(period, unweighted_n, primary_fn_events, kish_effective_n))
print(distribution_shift %>% select(target, osta_overlap_coefficient, age_weight_joint_overlap, domain_auc, target_in_propensity_common_support))
