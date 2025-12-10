# phase 6: anova modeling with blocking and permutation tests
#
# fitting three models: raw scale no blocking, raw scale with era blocking,
# log scale with era blocking, plus effect sizes and permutation tests

# sourcing configuration and utilities
source("config.R")
source("utils.R")

# loading required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(car)  # for Type III SS (Marginal) on unbalanced designs
  library(patchwork)
})

print_section_header("Phase 6: ANOVA Modeling With Blocking And Permutation Tests")

# 1. loading analysis-ready data from phase 5
input_file = file.path(RESULTS_DIR, "analysis_ready_dataset.csv")
analysis_data = load_csv(input_file)

cat(paste("Analysis-ready sample: n =", nrow(analysis_data), "\n"))

# 2. verifying factor variables exist from phase 5
# Phase 5 already created: review_type_factor, therapeutic_area_factor, regulatory_era_factor
cat("Factor variables from Phase 5:\n")
cat("  therapeutic_area_factor:", paste(levels(analysis_data$therapeutic_area_factor), collapse = ", "), "\n")
cat("  review_type_factor:", paste(levels(analysis_data$review_type_factor), collapse = ", "), "\n")
cat("  regulatory_era_factor:", paste(levels(analysis_data$regulatory_era_factor), collapse = ", "), "\n")

# verifying we have all 860 drugs (not 853!)
cat(sprintf("Total observations: %d (expecting 860)\n", nrow(analysis_data)))

# 3. model 1: raw scale, no blocking
model1 = lm(
  review_time_days_response ~ therapeutic_area_factor * review_type_factor,
  data = analysis_data
)

anova1 = car::Anova(model1, type = 3)

cat("\nModel 1: Two-way ANOVA (raw scale, no blocking, Type III SS)\n")
print(anova1)

# calculating effect sizes (eta-squared)
# Type III SS: extract by row name (skip Intercept)
ss_area = anova1["therapeutic_area_factor", "Sum Sq"]
ss_review = anova1["review_type_factor", "Sum Sq"]
ss_interaction = anova1["therapeutic_area_factor:review_type_factor", "Sum Sq"]
ss_residual = anova1["Residuals", "Sum Sq"]
ss_total = ss_area + ss_review + ss_interaction + ss_residual

eta_sq_area = ss_area / ss_total
eta_sq_review = ss_review / ss_total
eta_sq_interaction = ss_interaction / ss_total

cat("\nEffect sizes (eta-squared):\n")
cat(sprintf("  Therapeutic area: %.4f\n", eta_sq_area))
cat(sprintf("  Review type: %.4f\n", eta_sq_review))
cat(sprintf("  Interaction: %.4f\n", eta_sq_interaction))

# extracting key statistics
f_area_m1 = anova1["therapeutic_area_factor", "F value"]
p_area_m1 = anova1["therapeutic_area_factor", "Pr(>F)"]
f_review_m1 = anova1["review_type_factor", "F value"]
p_review_m1 = anova1["review_type_factor", "Pr(>F)"]
f_interaction_m1 = anova1["therapeutic_area_factor:review_type_factor", "F value"]
p_interaction_m1 = anova1["therapeutic_area_factor:review_type_factor", "Pr(>F)"]

cat("\nKey statistics:\n")
cat(sprintf("  Therapeutic area: F=%.2f, p=%.2e\n", f_area_m1, p_area_m1))
cat(sprintf("  Review type: F=%.2f, p=%.2e\n", f_review_m1, p_review_m1))
cat(sprintf("  Interaction: F=%.2f, p=%.2e\n", f_interaction_m1, p_interaction_m1))

# 4. model 2: raw scale, with era blocking
model2 = lm(
  review_time_days_response ~ therapeutic_area_factor * review_type_factor + regulatory_era_factor,
  data = analysis_data
)

anova2 = car::Anova(model2, type = 3)

cat("\nModel 2: Two-way ANOVA with blocking (raw scale, Type III SS)\n")
print(anova2)

# extracting key statistics
f_area_m2 = anova2["therapeutic_area_factor", "F value"]
p_area_m2 = anova2["therapeutic_area_factor", "Pr(>F)"]
f_review_m2 = anova2["review_type_factor", "F value"]
p_review_m2 = anova2["review_type_factor", "Pr(>F)"]
f_interaction_m2 = anova2["therapeutic_area_factor:review_type_factor", "F value"]
p_interaction_m2 = anova2["therapeutic_area_factor:review_type_factor", "Pr(>F)"]
f_era_m2 = anova2["regulatory_era_factor", "F value"]
p_era_m2 = anova2["regulatory_era_factor", "Pr(>F)"]

cat("\nKey statistics:\n")
cat(sprintf("  Therapeutic area: F=%.2f, p=%.2e\n", f_area_m2, p_area_m2))
cat(sprintf("  Review type: F=%.2f, p=%.2e\n", f_review_m2, p_review_m2))
cat(sprintf("  Interaction: F=%.2f, p=%.2e\n", f_interaction_m2, p_interaction_m2))
cat(sprintf("  Regulatory era (block): F=%.2f, p=%.2e\n", f_era_m2, p_era_m2))

# 5. model 3: log scale, with era blocking (PREFERRED MODEL)
model3 = lm(
  log_review_time_days_response ~ therapeutic_area_factor * review_type_factor + regulatory_era_factor,
  data = analysis_data
)

anova3 = car::Anova(model3, type = 3)

cat("\nModel 3: Two-way ANOVA with blocking (log scale, Type III SS, PREFERRED)\n")
print(anova3)

# extracting key statistics
f_area_m3 = anova3["therapeutic_area_factor", "F value"]
p_area_m3 = anova3["therapeutic_area_factor", "Pr(>F)"]
f_review_m3 = anova3["review_type_factor", "F value"]
p_review_m3 = anova3["review_type_factor", "Pr(>F)"]
f_interaction_m3 = anova3["therapeutic_area_factor:review_type_factor", "F value"]
p_interaction_m3 = anova3["therapeutic_area_factor:review_type_factor", "Pr(>F)"]
f_era_m3 = anova3["regulatory_era_factor", "F value"]
p_era_m3 = anova3["regulatory_era_factor", "Pr(>F)"]

# effect sizes for model 3 (log scale, blocked)
ss_area_m3 = anova3["therapeutic_area_factor", "Sum Sq"]
ss_review_m3 = anova3["review_type_factor", "Sum Sq"]
ss_interaction_m3 = anova3["therapeutic_area_factor:review_type_factor", "Sum Sq"]
ss_era_m3 = anova3["regulatory_era_factor", "Sum Sq"]
ss_resid_m3 = anova3["Residuals", "Sum Sq"]
ss_total_m3 = ss_area_m3 + ss_review_m3 + ss_interaction_m3 + ss_era_m3 + ss_resid_m3

eta_sq_area_m3 = ss_area_m3 / ss_total_m3
eta_sq_review_m3 = ss_review_m3 / ss_total_m3
eta_sq_interaction_m3 = ss_interaction_m3 / ss_total_m3
eta_sq_era_m3 = ss_era_m3 / ss_total_m3

cat("\nKey statistics:\n")
cat(sprintf("  Therapeutic area: F=%.2f, p=%.2e\n", f_area_m3, p_area_m3))
cat(sprintf("  Review type: F=%.2f, p=%.2e\n", f_review_m3, p_review_m3))
cat(sprintf("  Interaction: F=%.2f, p=%.2e\n", f_interaction_m3, p_interaction_m3))
cat(sprintf("  Regulatory era (block): F=%.2f, p=%.2e\n", f_era_m3, p_era_m3))

# 6. calculating estimated marginal means (cell means)
cell_means = analysis_data %>%
  group_by(therapeutic_area, review_type_simplified) %>%
  summarise(
    mean_review_time = mean(review_time_days_response, na.rm = TRUE),
    median_review_time = median(review_time_days_response, na.rm = TRUE),
    sd_review_time = sd(review_time_days_response, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

cat("\nCell means (raw scale):\n")
print(cell_means)

# calculating main effects
main_effect_area = analysis_data %>%
  group_by(therapeutic_area) %>%
  summarise(mean_review_time = mean(review_time_days_response, na.rm = TRUE), .groups = "drop")

main_effect_review = analysis_data %>%
  group_by(review_type_simplified) %>%
  summarise(mean_review_time = mean(review_time_days_response, na.rm = TRUE), .groups = "drop")

cat("\nMain effect means:\n")
cat("Therapeutic area:\n")
print(main_effect_area)
cat("\nReview designation:\n")
print(main_effect_review)

# 6.5. calculating practical effect sizes with bootstrap confidence intervals
cat("\n========================================\n")
cat("PRACTICAL EFFECT SIZES (GEOMETRIC MEANS)\n")
cat("========================================\n\n")

# calculating geometric means for Priority Review groups
oncology_priority_data = analysis_data %>%
  filter(therapeutic_area_factor == "Oncology", review_type_factor == "Priority") %>%
  pull(review_time_days_response)

other_priority_data = analysis_data %>%
  filter(therapeutic_area_factor == "Other", review_type_factor == "Priority") %>%
  pull(review_time_days_response)

# geometric mean calculation
geom_mean = function(x) {
  exp(mean(log(x), na.rm = TRUE))
}

oncology_priority_geom = geom_mean(oncology_priority_data)
other_priority_geom = geom_mean(other_priority_data)
ratio_point = oncology_priority_geom / other_priority_geom

cat(sprintf("Geometric mean review time (Priority Review):\n"))
cat(sprintf("  Oncology: %.1f days\n", oncology_priority_geom))
cat(sprintf("  Other: %.1f days\n", other_priority_geom))
cat(sprintf("  Ratio (Oncology/Other): %.3f\n\n", ratio_point))

# bootstrap confidence interval for ratio
set.seed(RANDOM_SEED)
n_boot = 10000
boot_ratios = numeric(n_boot)

cat(sprintf("Calculating bootstrap 95%% CI (%d iterations)...\n", n_boot))

for (i in 1:n_boot) {
  # resampling with replacement
  boot_onc = sample(oncology_priority_data, replace = TRUE)
  boot_oth = sample(other_priority_data, replace = TRUE)

  # calculating geometric means for bootstrap sample
  boot_geom_onc = geom_mean(boot_onc)
  boot_geom_oth = geom_mean(boot_oth)
  boot_ratios[i] = boot_geom_onc / boot_geom_oth
}

# percentile method confidence interval
ratio_ci_lower = quantile(boot_ratios, 0.025)
ratio_ci_upper = quantile(boot_ratios, 0.975)

cat("\nBootstrap 95% Confidence Interval for Ratio:\n")
cat(sprintf("  Point estimate: %.3f\n", ratio_point))
cat(sprintf("  95%% CI: [%.3f, %.3f]\n\n", ratio_ci_lower, ratio_ci_upper))

# interpreting as percent faster/slower
pct_diff = (1 - ratio_point) * 100
pct_ci_lower = (1 - ratio_ci_upper) * 100
pct_ci_upper = (1 - ratio_ci_lower) * 100

cat("PRACTICAL INTERPRETATION:\n")
if (ratio_point < 1) {
  cat(sprintf("Oncology Priority drugs approved %.1f%% FASTER than Other Priority drugs\n",
              abs(pct_diff)))
  cat(sprintf("  [95%% CI: %.1f%% to %.1f%% faster]\n\n",
              pct_ci_lower, pct_ci_upper))
} else {
  cat(sprintf("Oncology Priority drugs approved %.1f%% SLOWER than Other Priority drugs\n",
              pct_diff))
  cat(sprintf("  [95%% CI: %.1f%% to %.1f%% slower]\n\n",
              abs(pct_ci_upper), abs(pct_ci_lower)))
}

# 7. permutation test for interaction effect
set.seed(RANDOM_SEED)

# using interaction F-stat from Model 2 (Type III SS)
observed_f_stat = f_interaction_m2

# creating permutation distribution
n_perms = N_PERMUTATIONS
perm_f_stats = numeric(n_perms)

cat(paste("\nRunning", n_perms, "permutations...\n"))

for (i in 1:n_perms) {
  # permuting therapeutic area labels
  perm_data = analysis_data
  perm_data$therapeutic_area_factor = sample(analysis_data$therapeutic_area_factor)

  # fitting permuted model
  perm_model = lm(
    review_time_days_response ~ therapeutic_area_factor * review_type_factor + regulatory_era_factor,
    data = perm_data
  )

  # extracting F-statistic for interaction using Type III SS
  perm_anova = car::Anova(perm_model, type = 3)
  perm_f_stats[i] = perm_anova["therapeutic_area_factor:review_type_factor", "F value"]

  # progress indicator every 2000 iterations
  if (i %% 2000 == 0) {
    cat(paste("  Completed", i, "permutations\n"))
  }
}

# calculating p-value
perm_p_value = mean(perm_f_stats >= observed_f_stat)

cat(paste("\nPermutation test results:\n"))
cat(sprintf("  Observed F-statistic: %.2f\n", observed_f_stat))
cat(sprintf("  Permutation p-value: %.4f\n", perm_p_value))
cat(sprintf("  Parametric p-value: %.2e\n", p_interaction_m2))

# saving permutation results
perm_results = data.frame(
  permutation = 1:n_perms,
  f_statistic = perm_f_stats
)

perm_output_file = file.path(RESULTS_DIR, "permutation_test_results.csv")
save_csv(perm_results, perm_output_file)

# 8. visualization 1: fig1_interaction_plot.png
interaction_plot_data = cell_means %>%
  rename(Review_Type = review_type_simplified)

interaction_plot = ggplot(
  interaction_plot_data,
  aes(x = Review_Type, y = mean_review_time,
      color = therapeutic_area, group = therapeutic_area)
) +
  geom_line(linewidth = 1.5, alpha = 0.8) +
  geom_point(size = 5, alpha = 0.9) +
  geom_errorbar(
    aes(ymin = mean_review_time - sd_review_time / sqrt(n),
        ymax = mean_review_time + sd_review_time / sqrt(n)),
    width = 0.15,
    linewidth = 1.2,
    alpha = 0.7
  ) +
  labs(
    title = "Interaction Effect: Therapeutic Area x Review Type",
    subtitle = sprintf("F=%.2f, p<%.2e", f_interaction_m2, p_interaction_m2),
    x = "Review Designation",
    y = "Mean Review Time (Days)",
    color = "Therapeutic Area"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text = element_text(size = 12)
  ) +
  scale_color_manual(values = c("Other" = "#1f77b4", "Oncology" = "#ff7f0e"))

print(interaction_plot)  # Display in R session

ggsave(
  file.path(FIGURES_DIR, "fig1_interaction_plot.png"),
  plot = interaction_plot,
  width = FIGURE_WIDTH,
  height = FIGURE_HEIGHT,
  dpi = DPI
)

cat("Saved: fig1_interaction_plot.png\n")

# 9. visualization 2: fig2_era_trends_plot.png
# CRITICAL FIX: Era ordering Post→Pre (modern first)
era_means = analysis_data %>%
  group_by(regulatory_era_factor, therapeutic_area, review_type_simplified) %>%
  summarise(mean_review_time = mean(review_time_days_response, na.rm = TRUE), .groups = "drop") %>%
  mutate(group = paste(therapeutic_area, review_type_simplified, sep = "-"))

# reorder factor levels for display: Post→Pre (modern first)
era_means$regulatory_era_factor = factor(
  era_means$regulatory_era_factor,
  levels = c("Post-FDASIA", "Mid-PDUFA", "Early-PDUFA", "Pre-PDUFA")
)

era_trends_plot = ggplot(
  era_means,
  aes(x = regulatory_era_factor, y = mean_review_time, color = group, group = group)
) +
  geom_line(linewidth = 1.5, alpha = 0.8) +
  geom_point(size = 4, alpha = 0.9) +
  labs(
    title = "Temporal Trends: Mean Review Time By Regulatory Era",
    x = "Regulatory Era",
    y = "Mean Review Time (Days)",
    color = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  scale_color_manual(
    values = c(
      "Oncology-Priority" = "#d95f02",
      "Oncology-Standard" = "#a05dbb",
      "Other-Priority" = "#1b9e77",
      "Other-Standard" = "#4b8ad1"
    )
  )

print(era_trends_plot)  # Display in R session

ggsave(
  file.path(FIGURES_DIR, "fig2_era_trends_plot.png"),
  plot = era_trends_plot,
  width = FIGURE_WIDTH,
  height = FIGURE_HEIGHT,
  dpi = DPI
)

cat("Saved: fig2_era_trends_plot.png\n")

# 10. visualization 3: fig3_permutation_histograms.png
perm_hist = ggplot(perm_results, aes(x = f_statistic)) +
  geom_histogram(bins = 50, fill = "#74a9cf", color = "white", alpha = 0.9, linewidth = 0.2) +
  geom_vline(xintercept = observed_f_stat, color = "#d62728", linewidth = 2, linetype = "dashed") +
  annotate(
    "text",
    x = observed_f_stat,
    y = Inf,
    label = sprintf("Observed F = %.2f", observed_f_stat),
    vjust = 1.5,
    hjust = -0.1,
    color = "#d62728",
    size = 5,
    fontface = "bold"
  ) +
  labs(
    title = "Permutation Test Distribution For Interaction Effect",
    subtitle = sprintf("p-value = %.4f (based on %d permutations)", perm_p_value, n_perms),
    x = "F-Statistic",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12)
  )

print(perm_hist)  # Display in R session

ggsave(
  file.path(FIGURES_DIR, "fig3_permutation_histograms.png"),
  plot = perm_hist,
  width = FIGURE_WIDTH,
  height = FIGURE_HEIGHT,
  dpi = DPI
)

cat("Saved: fig3_permutation_histograms.png\n")

# 11. visualization 4: fig4_cell_means_barchart.png
cell_means_plot = ggplot(
  interaction_plot_data,
  aes(x = therapeutic_area, y = mean_review_time, fill = Review_Type)
) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8,
           alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = mean_review_time - 1.96 * sd_review_time / sqrt(n),
        ymax = mean_review_time + 1.96 * sd_review_time / sqrt(n)),
    position = position_dodge(width = 0.9),
    width = 0.25,
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = sprintf("n=%d", n)),
    position = position_dodge(width = 0.9),
    vjust = -1.5,
    size = 4,
    fontface = "bold"
  ) +
  labs(
    title = "Cell Means With 95% Confidence Intervals",
    x = "Therapeutic Area",
    y = "Mean Review Time (Days)",
    fill = "Review Designation"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text = element_text(size = 12)
  ) +
  scale_fill_manual(values = c("Standard" = "#8da0cb", "Priority" = "#fc8d62"))

print(cell_means_plot)  # Display in R session

ggsave(
  file.path(FIGURES_DIR, "fig4_cell_means_barchart.png"),
  plot = cell_means_plot,
  width = FIGURE_WIDTH,
  height = FIGURE_HEIGHT,
  dpi = DPI
)

cat("Saved: fig4_cell_means_barchart.png\n")

# 12. visualization 5: fig5_effect_sizes_barchart.png
effect_sizes_data = data.frame(
  effect = c("Therapeutic Area", "Review Type", "Interaction"),
  eta_squared = c(eta_sq_area, eta_sq_review, eta_sq_interaction)
)

effect_sizes_plot = ggplot(
  effect_sizes_data,
  aes(x = reorder(effect, -eta_squared), y = eta_squared, fill = effect)
) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.4f", eta_squared)),
    vjust = -0.5,
    size = 5,
    fontface = "bold"
  ) +
  labs(
    title = "Effect Sizes (Eta-Squared)",
    subtitle = "Model 1: Raw Scale, No Blocking",
    x = "Effect",
    y = "Eta-Squared"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12)
  ) +
  scale_fill_manual(
    values = c(
      "Therapeutic Area" = "#4b8ad1",
      "Review Type" = "#f28e2b",
      "Interaction" = "#59a14f"
    )
  )

print(effect_sizes_plot)  # Display in R session

ggsave(
  file.path(FIGURES_DIR, "fig5_effect_sizes_barchart.png"),
  plot = effect_sizes_plot,
  width = FIGURE_WIDTH,
  height = FIGURE_HEIGHT,
  dpi = DPI
)

cat("Saved: fig5_effect_sizes_barchart.png\n")

# additional visualization 5 (alt): fig5_effect_sizes_barchart_alt.png (blocked log-scale η², horizontal)
effect_sizes_m3 = data.frame(
  factor = c("Regulatory Era (Block)", "Review Type", "Therapeutic Area", "Area × Review Interaction"),
  eta_sq = c(eta_sq_era_m3, eta_sq_review_m3, eta_sq_area_m3, eta_sq_interaction_m3)
) %>%
  mutate(
    eta_pct = eta_sq * 100,
    magnitude = case_when(
      eta_sq >= 0.14 ~ "large",
      eta_sq >= 0.06 ~ "medium",
      eta_sq >= 0.01 ~ "small",
      TRUE ~ "negligible"
    ),
    magnitude = factor(magnitude, levels = c("large", "medium", "small", "negligible"))
  )

magnitude_palette = c(
  "large" = "#3b5c7a",
  "medium" = "#587b9d",
  "small" = "#9ecae1",
  "negligible" = "#cbd5e8"
)

effect_sizes_plot_alt = ggplot(effect_sizes_m3, aes(x = eta_pct, y = factor, fill = magnitude)) +
  geom_col(width = 0.6, alpha = 0.95, color = "white") +
  geom_text(
    aes(label = sprintf("%.1f%% (%s)", eta_pct, magnitude)),
    hjust = -0.1,
    size = 4,
    fontface = "bold"
  ) +
  scale_fill_manual(values = magnitude_palette, name = "Effect size (η²)") +
  labs(
    title = "Effect Sizes: Blocked ANOVA (Log Scale)",
    subtitle = "Variance explained (η²) by factor",
    x = "Effect Size η² (%)",
    y = "Factor"
  ) +
  coord_cartesian(xlim = c(0, max(effect_sizes_m3$eta_pct) * 1.3)) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.text = element_text(size = 11)
  )

print(effect_sizes_plot_alt)

ggsave(
  file.path(FIGURES_DIR, "fig5_effect_sizes_barchart_alt.png"),
  plot = effect_sizes_plot_alt,
  width = FIGURE_WIDTH,
  height = FIGURE_HEIGHT,
  dpi = DPI
)

# 13. visualization 6: fig6_residual_diagnostics.png
residuals_data = data.frame(
  fitted = fitted(model3),
  residuals = residuals(model3),
  standardized_residuals = rstandard(model3)
)

png(
  file.path(FIGURES_DIR, "fig6_residual_diagnostics.png"),
  width = FIGURE_WIDTH * DPI,
  height = FIGURE_HEIGHT * DPI,
  res = DPI
)

par(mfrow = c(2, 2))

# plot 1: residuals vs fitted
plot(
  residuals_data$fitted,
  residuals_data$residuals,
  main = "Residuals vs Fitted",
  xlab = "Fitted Values",
  ylab = "Residuals",
  pch = 16,
  col = rgb(0, 0, 0, 0.3)
)
abline(h = 0, col = "red", lwd = 2, lty = 2)

# plot 2: QQ plot
qqnorm(residuals_data$standardized_residuals, main = "Normal Q-Q Plot", pch = 16, col = rgb(0, 0, 0, 0.3))
qqline(residuals_data$standardized_residuals, col = "red", lwd = 2)

# plot 3: scale-location
plot(
  residuals_data$fitted,
  sqrt(abs(residuals_data$standardized_residuals)),
  main = "Scale-Location",
  xlab = "Fitted Values",
  ylab = expression(sqrt("|Standardized Residuals|")),
  pch = 16,
  col = rgb(0, 0, 0, 0.3)
)

# plot 4: residuals histogram
hist(
  residuals_data$residuals,
  breaks = 30,
  main = "Histogram of Residuals",
  xlab = "Residuals",
  col = "skyblue",
  border = "black"
)

dev.off()

cat("Saved: fig6_residual_diagnostics.png\n")

# additional residual diagnostics (alt) with stats and density overlay (log scale)
residuals_log = residuals(model3)
shapiro_log = shapiro.test(residuals_log)
w_stat_log = shapiro_log$statistic
p_value_log = shapiro_log$p.value
mean_log = mean(residuals_log)
sd_log = sd(residuals_log)
n_log = length(residuals_log)

qq_log_alt = ggplot(data.frame(residuals_log = residuals_log), aes(sample = residuals_log)) +
  stat_qq(color = "#1f77b4", size = 1.4, alpha = 0.85) +
  stat_qq_line(color = "#d62728", linewidth = 1.2) +
  annotate(
    "label",
    x = -2.5,
    y = 2.5,
    label = sprintf("Shapiro-Wilk:\nW = %.4f\np = %.2e", w_stat_log, p_value_log),
    fill = "white",
    color = "black",
    label.size = 0.3
  ) +
  labs(
    title = "QQ-Plot: Log Scale Residuals",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

hist_log_alt = ggplot(data.frame(residuals_log = residuals_log), aes(x = residuals_log)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60,
                 fill = "#6c8ebf",
                 color = "white",
                 alpha = 0.9) +
  stat_function(fun = dnorm,
                args = list(mean = mean_log, sd = sd_log),
                color = "#d62728",
                linewidth = 1.2) +
  annotate(
    "label",
    x = min(residuals_log),
    y = max(density(residuals_log)$y),
    hjust = 0,
    label = sprintf("Mean = %.4f\nSD = %.4f\nn = %d", mean_log, sd_log, n_log),
    fill = "white",
    color = "black",
    label.size = 0.3
  ) +
  labs(
    title = "Histogram: Log Scale Residuals",
    x = "Residuals (log scale)",
    y = "Density"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

resid_alt = qq_log_alt + hist_log_alt + plot_layout(ncol = 2)

ggsave(
  file.path(FIGURES_DIR, "fig6_residual_diagnostics_alt.png"),
  plot = resid_alt,
  width = FIGURE_WIDTH,
  height = FIGURE_HEIGHT,
  dpi = DPI
)

# 14. saving model comparison results
model_comparison = data.frame(
  model = c("Model 1", "Model 1", "Model 1",
            "Model 2", "Model 2", "Model 2", "Model 2",
            "Model 3", "Model 3", "Model 3", "Model 3"),
  effect = c("Therapeutic Area", "Review Type", "Interaction",
             "Therapeutic Area", "Review Type", "Interaction", "Regulatory Era",
             "Therapeutic Area", "Review Type", "Interaction", "Regulatory Era"),
  f_statistic = c(f_area_m1, f_review_m1, f_interaction_m1,
                  f_area_m2, f_review_m2, f_interaction_m2, f_era_m2,
                  f_area_m3, f_review_m3, f_interaction_m3, f_era_m3),
  p_value = c(p_area_m1, p_review_m1, p_interaction_m1,
              p_area_m2, p_review_m2, p_interaction_m2, p_era_m2,
              p_area_m3, p_review_m3, p_interaction_m3, p_era_m3)
)

comparison_output_file = file.path(RESULTS_DIR, "model_comparison.csv")
save_csv(model_comparison, comparison_output_file)

# saving cell means
cell_means_output_file = file.path(RESULTS_DIR, "cell_means.csv")
save_csv(cell_means, cell_means_output_file)

cat("\nPhase 6 complete\n")
