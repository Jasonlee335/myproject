# Paired allocation contrasts and final figures from the prespecified OSTA ranking.
options(stringsAsFactors = FALSE, scipen = 999)
stopifnot(all(vapply(c("survey", "haven", "ggplot2", "patchwork", "svglite"),
                    requireNamespace, logical(1), quietly = TRUE)))
suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
options(survey.lonely.psu = "adjust")
set.seed(20261028)
project <- Sys.getenv("DXA_PROJECT_DIR", file.path(getwd(), "DXA_transportability_results"))
results <- file.path(project, "reviewer1_revision", "analysis_results")
figures <- file.path(project, "reviewer1_revision", "figures")
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
B <- as.integer(Sys.getenv("REVISION_BOOTSTRAP_REPS", "1000"))
stopifnot(is.finite(B), B >= 100)
periods <- c(development = "NHANES 2005-2010", temporal_validation_1 = "NHANES 2013-2014",
             temporal_validation_2 = "NHANES 2017-Mar 2020",
             historical_geographic_validation = "KNHANES 2008-2011",
             contemporary_temporal_validation = "KNHANES 2024")
bands <- c("50-59", "60-69", "70-79")
capacity <- function(score, weight) {
  order_id <- order(score, decreasing = TRUE)
  cut <- score[order_id[which(cumsum(weight[order_id]) >= .2 * sum(weight))[1]]]
  fraction <- (.2 * sum(weight) - sum(weight[score > cut])) / sum(weight[score == cut])
  ifelse(score > cut, 1, ifelse(score == cut, min(1, max(0, fraction)), 0))
}
sensitivity <- function(y, selected, weight) sum(y * selected * weight) / sum(y * weight)
allocation <- function(d, weight) {
  within <- numeric(nrow(d))
  for (band in bands) {
    k <- d$age_band == band
    within[k] <- capacity(d$risk[k], weight[k])
  }
  list(pooled = capacity(d$risk, weight), within = within)
}
contrast <- function(d, y, weight) {
  a <- allocation(d, weight)
  answer <- c(overall = sensitivity(y, a$within, weight) - sensitivity(y, a$pooled, weight))
  for (band in bands) {
    k <- d$age_band == band
    answer[band] <- sensitivity(y[k], a$within[k], weight[k]) -
      sensitivity(y[k], a$pooled[k], weight[k])
  }
  answer
}
interval <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < .9 * B) return(c(NA_real_, NA_real_))
  unname(quantile(x, c(.025, .975), type = 6))
}
im <- readRDS(file.path(project, "data_clean", "imaging_cross_country_analysis_ready.rds"))
im$period <- unname(periods[im$analysis_role])
im$age_band <- cut(im$age, c(50, 60, 70, 80), right = FALSE, labels = bands)
im$risk <- -im$osta_raw
sets <- split(im, factor(im$period, levels = unname(periods)))
# Shared survey replicates make within-minus-pooled contrasts paired.
replicates <- lapply(sets, function(d) {
  design <- svydesign(ids = ~psu, strata = ~strata, weights = ~survey_weight, nest = TRUE, data = d)
  as.matrix(weights(as.svrepdesign(design, type = "bootstrap", replicates = B, mse = TRUE),
                    type = "analysis"))
})
rows <- list()
for (period in names(sets)) {
  message("Paired allocation contrasts: ", period)
  d <- sets[[period]]
  for (reference in c("common_US", if (grepl("KNHANES", period)) "local")) {
    y <- if (reference == "common_US") d$osteoporosis_fn_us else d$osteoporosis_fn_local
    point <- contrast(d, y, d$survey_weight)
    boot <- vapply(seq_len(B), function(b) contrast(d, y, replicates[[period]][, b]), numeric(4))
    for (band in names(point)) {
      limits <- interval(boot[band, ])
      k <- if (band == "overall") rep(TRUE, nrow(d)) else d$age_band == band
      rows[[length(rows) + 1]] <- data.frame(
        period = period, reference = reference, age_band = band,
        contrast = "within_age20_minus_pooled20_sensitivity",
        estimate = unname(point[band]), lower = limits[1], upper = limits[2],
        n = sum(k), cases = sum(y[k] == 1), bootstrap_replicates = B,
        successful_replicates = sum(is.finite(boot[band, ])))
    }
  }
}
paired <- do.call(rbind, rows)
write.csv(paired, file.path(results, "R1_paired_allocation_contrasts.csv"), row.names = FALSE)
age <- read.csv(file.path(results, "R1_age_stratified_performance.csv"))
age_auc <- subset(age, strategy == "pooled20" & metric == "auc" &
                    (reference == "common_US" | grepl("KNHANES", period)))
age_auc <- age_auc[, setdiff(names(age_auc), c("strategy", "metric", "population_share", "case_share"))]
write.csv(age_auc, file.path(results, "R1_age_auc_summary.csv"), row.names = FALSE)

raw_root <- file.path(Sys.getenv("DXA_DATA_ROOT", file.path(getwd(), "data_raw")), "charls_official/raw_extracted/2011")
audit <- list()
modules <- list.files(raw_root, pattern = "[.]dta$", full.names = TRUE)
for (file in modules) {
  d <- haven::read_dta(file, n_max = 0)
  labels <- vapply(d, function(x) { z <- attr(x, "label"); if (is.null(z)) "" else as.character(z) }, character(1))
  hits <- grepl("fracture|steroid|prednis|glucocort|rheuma|arthritis|alcohol|liquor|wine|beer|parent|mother|father",
                paste(names(d), labels), ignore.case = TRUE)
  if (any(hits)) audit[[length(audit) + 1]] <- data.frame(module = basename(file), variable = names(d)[hits], label = labels[hits])
}
write.csv(do.call(rbind, audit), file.path(results, "R1_FRAX_all_baseline_modules_audit.csv"), row.names = FALSE)
write.csv(data.frame(module = basename(modules), inspected = TRUE, records_read = 0L),
          file.path(results, "R1_FRAX_audit_module_manifest.csv"), row.names = FALSE)

performance <- read.csv(file.path(project, "step3_results", "Table4_model_performance.csv"))
pooled_auc <- subset(performance, model == "OSTA" & source_key == "NH_US" &
                      target_outcome == "osteoporosis_fn_us")
stopifnot(nrow(pooled_auc) == 5)
pooled_auc$period <- pooled_auc$target_period
short <- c("US 2005-2010", "US 2013-2014", "US 2017-2020", "Korea 2008-2011", "Korea 2024")
labels <- setNames(short, unname(periods))
y_factor <- function(period) factor(labels[period], levels = rev(short))
pooled_auc$period_label <- y_factor(pooled_auc$period)
joint <- subset(read.csv(file.path(results, "R1_joint_source_target_bootstrap.csv")),
                bootstrap == "source_and_target_joint" & metric == "referral_share")
joint <- rbind(data.frame(period = unname(periods[1]), bootstrap = "source_and_target_joint",
                          metric = "referral_share", estimate = .2, lower = .2, upper = .2,
                          successful_replicates = B), joint)
joint$period_label <- y_factor(joint$period)
caps <- subset(read.csv(file.path(results, "R1_capacity_tradeoffs.csv")),
               reference == "common_US" & metric == "sensitivity" &
                 strategy %in% c("pooled20", "age_stratified20"))
caps$period_label <- y_factor(caps$period)
age_sens <- subset(age, reference == "common_US" & grepl("KNHANES", period) &
                     metric == "sensitivity" & strategy %in% c("pooled20", "age_stratified20"))
strategy_labels <- c(pooled20 = "Pooled 20%", age_stratified20 = "20% within age band")
for (name in c("caps", "age_sens")) {
  d <- get(name)
  d$strategy_label <- factor(strategy_labels[d$strategy], levels = unname(strategy_labels))
  assign(name, d)
}
age_sens$age_band <- factor(age_sens$age_band, levels = rev(bands))
age_sens$period_label <- factor(ifelse(age_sens$period == "KNHANES 2008-2011", "Korea 2008-2011", "Korea 2024"),
                              levels = c("Korea 2008-2011", "Korea 2024"))
palette <- c("Pooled 20%" = "#0072B2", "20% within age band" = "#D55E00")
symbols <- c("Pooled 20%" = 16, "20% within age band" = 17)
theme_set(theme_classic(base_size = 8.5, base_family = "Arial"))
common <- theme(axis.text = element_text(colour = "black", size = 8),
                axis.title = element_text(size = 8.5),
                plot.tag = element_text(face = "bold", size = 11),
                plot.tag.position = "topleft", plot.margin = margin(5, 6, 4, 4),
                legend.position = "bottom", legend.title = element_blank(),
                legend.text = element_text(size = 8), strip.background = element_blank(),
                strip.text = element_text(size = 8, face = "bold"),
                panel.grid.major.x = element_line(colour = "grey90", linewidth = .25))
p_a <- ggplot(pooled_auc, aes(x = auc, y = period_label)) +
  geom_errorbar(aes(xmin = auc_lower, xmax = auc_upper), orientation = "y", width = .12, linewidth = .4) +
  geom_point(size = 1.7, colour = "#0072B2") +
  scale_x_continuous(limits = c(.65, .9), breaks = seq(.65, .9, .05)) +
  labs(x = "Pooled AUC (95% CI)", y = NULL, tag = "a") + common
p_b <- ggplot(joint, aes(x = 100 * estimate, y = period_label)) +
  geom_vline(xintercept = 20, linetype = "dashed", colour = "grey35", linewidth = .35) +
  geom_errorbar(aes(xmin = 100 * lower, xmax = 100 * upper), orientation = "y", width = .12, linewidth = .4) +
  geom_point(size = 1.8, shape = 18, colour = "#0072B2") +
  scale_x_continuous(limits = c(10, 60), breaks = seq(10, 60, 10)) +
  labs(x = "Source-cutoff referral (%)", y = NULL, tag = "b") + common
dodge <- position_dodge(width = .42)
p_c <- ggplot(age_sens, aes(x = 100 * estimate, y = age_band, colour = strategy_label,
                           shape = strategy_label, group = strategy_label)) +
  geom_errorbar(aes(xmin = 100 * lower, xmax = 100 * upper), orientation = "y",
                position = dodge, width = .1, linewidth = .4) +
  geom_point(position = dodge, size = 1.8) +
  facet_wrap(~period_label, ncol = 1) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
  scale_colour_manual(values = palette) + scale_shape_manual(values = symbols) +
  labs(x = "Age-specific sensitivity (%)", y = "Age (years)", tag = "c") + common +
  theme(strip.text = element_text(size = 8))
p_d <- ggplot(caps, aes(x = 100 * estimate, y = period_label, colour = strategy_label,
                       shape = strategy_label, group = strategy_label)) +
  geom_errorbar(aes(xmin = 100 * lower, xmax = 100 * upper), orientation = "y",
                position = dodge, width = .12, linewidth = .4) +
  geom_point(position = dodge, size = 1.8) +
  scale_x_continuous(limits = c(20, 80), breaks = c(20, 40, 60, 80)) +
  scale_colour_manual(values = palette) + scale_shape_manual(values = symbols) +
  labs(x = "Overall sensitivity (%)", y = NULL, tag = "d") + common
composite <- ((p_a | p_b) / (p_c | p_d)) + plot_layout(heights = c(1, 1.2), guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(figures, "Figure2_Transportability_Capacity.png"), composite,
       width = 6.3, height = 6.45, units = "in", dpi = 600, bg = "white")
ggsave(file.path(figures, "Figure2_Transportability_Capacity.tiff"), composite,
       width = 6.3, height = 6.45, units = "in", dpi = 600, compression = "lzw", bg = "white")
ggsave(file.path(figures, "Figure2_Transportability_Capacity.svg"), composite,
       width = 6.3, height = 6.45, units = "in", device = svglite::svglite, bg = "white")
ggsave(file.path(figures, "Figure2_Transportability_Capacity.pdf"), composite,
       width = 6.3, height = 6.45, units = "in", device = cairo_pdf, bg = "white")
gray_a <- p_a
gray_a$layers[[2]]$aes_params$colour <- "grey15"
gray_b <- p_b
gray_b$layers[[3]]$aes_params$colour <- "grey15"
gray_values <- c("Pooled 20%" = "grey15", "20% within age band" = "grey55")
gray_c <- p_c + scale_colour_manual(values = gray_values)
gray_d <- p_d + scale_colour_manual(values = gray_values)
gray <- ((gray_a | gray_b) / (gray_c | gray_d)) +
  plot_layout(heights = c(1, 1.2), guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(figures, "Figure2_Transportability_Capacity_grayscale.png"), gray,
       width = 6.3, height = 6.45, units = "in", dpi = 300, bg = "white")
benchmark <- read.csv(file.path(project, "step5_results", "Table7_capacity_strategy_performance.csv"))
benchmark <- subset(benchmark, abs(intended_capacity - .2) < 1e-8 &
                      strategy %in% c("Target-percentile OSTA", "Target-percentile spline", "Target-percentile ridge"))
benchmark$model <- factor(benchmark$strategy,
                          levels = c("Target-percentile OSTA", "Target-percentile spline", "Target-percentile ridge"),
                          labels = c("OSTA", "Age-weight spline", "Ridge"))
benchmark$period_label <- y_factor(benchmark$target_period)
model_colours <- c("OSTA" = "#0072B2", "Age-weight spline" = "#D55E00", "Ridge" = "#009E73")
model_shapes <- c("OSTA" = 16, "Age-weight spline" = 17, "Ridge" = 15)
benchmark_plot <- function(data, tag) {
  ggplot(data, aes(x = 100 * sensitivity, y = period_label, colour = model, shape = model, group = model)) +
    geom_errorbar(aes(xmin = 100 * sensitivity_lower, xmax = 100 * sensitivity_upper),
                  orientation = "y", position = position_dodge(.55), width = .1, linewidth = .4) +
    geom_point(position = position_dodge(.55), size = 1.7) +
    scale_x_continuous(limits = c(25, 80), breaks = c(30, 50, 70)) +
    scale_colour_manual(values = model_colours) + scale_shape_manual(values = model_shapes) +
    labs(x = "Sensitivity at 20% referral (%)", y = NULL, tag = tag) + common
}
common_benchmark <- subset(benchmark, outcome == "osteoporosis_fn_us")
local_benchmark <- subset(benchmark, outcome == "osteoporosis_fn_local" & grepl("KNHANES", target_period))
stopifnot(nrow(common_benchmark) == 15, nrow(local_benchmark) == 6)
s5 <- (benchmark_plot(common_benchmark, "a") | benchmark_plot(local_benchmark, "b")) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(figures, "FigureS5_Model_Capacity_Comparison.png"), s5,
       width = 6.35, height = 3.45, units = "in", dpi = 600, bg = "white")
ggsave(file.path(figures, "FigureS5_Model_Capacity_Comparison.svg"), s5,
       width = 6.35, height = 3.45, units = "in", device = svglite::svglite, bg = "white")
write.csv(data.frame(file = c("R1_paired_allocation_contrasts.csv", "R1_age_auc_summary.csv"),
                     bootstrap_replicates = B, seed = 20261028,
                     scope = "Post hoc exploratory; no multiplicity adjustment; paired survey replicates"),
          file.path(results, "R1_synthesis_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(results, "R1_synthesis_sessionInfo.txt"))
source(file.path(Sys.getenv("DXA_CODE_DIR", unset = getwd()), "09b_supplement_figure_layout.R"), encoding = "UTF-8")
cat("Paired contrasts, FRAX metadata audit, and Figure 2 completed.\n")
