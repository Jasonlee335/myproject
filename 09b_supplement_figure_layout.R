# Layout-only refresh from existing aggregate robustness results.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
project <- Sys.getenv("DXA_PROJECT_DIR", file.path(getwd(), "DXA_transportability_results"))
figures <- file.path(project, "reviewer1_revision", "figures")
img <- read.csv(file.path(project, "step7_results", "Table13_imaging_site_subgroup_robustness.csv"))
img <- subset(img, (cohort == "NHANES" & outcome_reference == "US common") |
                (cohort == "KNHANES" & outcome_reference == "Asian local"))
charls <- read.csv(file.path(project, "step7_results", "Table14_charls_weight_subgroup_sensitivity.csv"))
land <- read.csv(file.path(project, "step6_results", "TableS11_sequential_landmarks.csv"))
land <- subset(land, !is.na(landmark_year))
periods <- c("NHANES 2005-2010", "NHANES 2013-2014", "NHANES 2017-Mar 2020",
             "KNHANES 2008-2011", "KNHANES 2024")
short <- c("US 2005-2010", "US 2013-2014", "US 2017-2020", "Korea 2008-2011", "Korea 2024")
names(short) <- periods
cols <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
theme_set(theme_classic(base_size = 8.5, base_family = "Arial"))
style <- theme(axis.text = element_text(size = 8, colour = "black"),
               axis.title = element_text(size = 8.5),
               plot.title = element_text(size = 8.5, face = "bold"),
               plot.tag = element_text(size = 11, face = "bold"),
               plot.tag.position = "topleft", plot.margin = margin(5, 5, 5, 4),
               legend.position = "bottom", legend.title = element_blank(),
               legend.text = element_text(size = 8), legend.key.size = unit(.28, "cm"),
               panel.grid.major.x = element_line(colour = "grey90", linewidth = .25))
sites <- c("Femoral neck", "Total hip", "Lumbar spine", "Any complete site")
a <- subset(img, subgroup == "All eligible")
a$period_label <- factor(short[a$period], levels = rev(short))
a$site <- factor(a$site, levels = sites)
pd <- position_dodge(.60)
p1 <- ggplot(a, aes(x = auc, y = period_label, colour = site, shape = site, group = site)) +
  geom_errorbar(aes(xmin = auc_lower, xmax = auc_upper), orientation = "y",
                position = pd, width = .10, linewidth = .35) +
  geom_point(position = pd, size = 1.5) +
  scale_colour_manual(values = setNames(cols, sites)) +
  scale_shape_manual(values = setNames(c(16, 17, 15, 18), sites)) +
  scale_x_continuous(limits = c(min(a$auc_lower, na.rm = TRUE) - .01, 1),
                     breaks = c(.6, .7, .8, .9, 1)) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  labs(x = "OSTA AUC (95% CI)", y = NULL, title = "Skeletal-site definitions", tag = "a") + style
groups <- c("Age 50-64", "Age 65-79", "BMI <25 kg/m2", "BMI >=25 kg/m2", "Postmenopausal observed")
pick <- periods[c(1, 4, 5)]
b <- subset(img, site == "Femoral neck" & subgroup %in% groups & period %in% pick)
group_labels <- c("Age 50-64", "Age 65-79", "BMI <25", "BMI >=25", "Reported\nmenopause*")
b$group_label <- factor(group_labels[match(b$subgroup, groups)], levels = rev(group_labels))
b$period_label <- factor(short[b$period], levels = short[pick])
pc <- setNames(cols[c(1, 3, 2)], short[pick])
p2 <- ggplot(b, aes(x = auc, y = group_label, colour = period_label,
                    shape = period_label, group = period_label)) +
  geom_errorbar(aes(xmin = auc_lower, xmax = auc_upper), orientation = "y",
                position = position_dodge(.5), width = .1, linewidth = .35) +
  geom_point(position = position_dodge(.5), size = 1.5) +
  scale_colour_manual(values = pc) +
  scale_shape_manual(values = setNames(c(16, 15, 17), short[pick])) +
  scale_x_continuous(limits = c(.58, .92), breaks = c(.6, .7, .8, .9)) +
  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  labs(x = "OSTA AUC (95% CI)", y = NULL, title = "Imaging subgroups", tag = "b") + style
c_labels <- ifelse(charls$analysis == "Prespecified subgroup", charls$subgroup, charls$analysis)
c_labels <- sub("Complete-case primary model", "Primary model", c_labels, fixed = TRUE)
c_labels <- sub("Individual-weight complete-case sensitivity", "Individual weight", c_labels, fixed = TRUE)
c_labels <- sub(" sensitivity", "", c_labels, fixed = TRUE)
c_labels <- sub(" kg/m2", "", c_labels, fixed = TRUE)
charls$label <- factor(c_labels, levels = rev(c_labels))
p3 <- ggplot(charls, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .10, linewidth = .4, colour = cols[1]) +
  geom_point(size = 1.6, colour = cols[1]) +
  scale_x_continuous(limits = c(.72, 2.22), breaks = c(1, 1.5, 2)) +
  labs(x = "HR per SD (95% CI)", y = NULL, title = "CHARLS weighting and subgroups", tag = "c") + style
land$year <- factor(land$landmark_year, levels = rev(sort(land$landmark_year)))
p4 <- ggplot(land, aes(x = estimate, y = year)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .3) +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .10, linewidth = .4, colour = cols[3]) +
  geom_point(size = 1.6, shape = 15, colour = cols[3]) +
  scale_x_continuous(limits = c(.8, 2.55), breaks = c(1, 1.5, 2, 2.5)) +
  labs(x = "HR for top 20% (95% CI)", y = "Updated landmark year",
       title = "Sequential landmarks", tag = "d") + style
plot <- (p1 | p2) / (p3 | p4) + plot_layout(heights = c(1.12, 1))
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figures, "FigureS4_Robustness.png"), plot,
       width = 6.35, height = 6.2, units = "in", dpi = 600, bg = "white")
ggsave(file.path(figures, "FigureS4_Robustness.svg"), plot,
       width = 6.35, height = 6.2, units = "in", device = svglite::svglite, bg = "white")
ggsave(file.path(figures, "FigureS4_Robustness.pdf"), plot,
       width = 6.35, height = 6.2, units = "in", device = cairo_pdf, bg = "white")
cat("Figure S4 layout refreshed; no analysis estimates changed.\n")
