#!/usr/bin/env Rscript

# HE Regression Results Visualization
# This script creates plots for Haseman-Elston regression results

library(ggplot2)
library(dplyr)
library(readr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript plot_he_regression.R <he_results_file> [output_prefix]\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  he_results_file - Tab-delimited file with HE regression results\n") 
  cat("  output_prefix   - Output file prefix (default: he_plots)\n")
  cat("\n")
  cat("Example:\n")
  cat("  Rscript plot_he_regression.R he_regression_results.txt my_plots\n")
  quit(status = 1)
}

results_file <- args[1]
output_prefix <- ifelse(length(args) >= 2, args[2], "he_plots")

# Read results
cat("Reading HE regression results from:", results_file, "\n")

tryCatch({
  results <- read_tsv(results_file, show_col_types = FALSE)
  cat("Loaded", nrow(results), "intervals\n")
}, error = function(e) {
  cat("Error reading results file:", e$message, "\n")
  quit(status = 1)
})

# Filter out intervals with no data
results_clean <- results %>%
  filter(!is.na(slope) & n_pairs >= 2) %>%
  mutate(
    heritability = slope,  # HE estimate is the slope directly
    h2_ci_lower = slope_ci_lower,
    h2_ci_upper = slope_ci_upper,
    interval_mid = as.numeric(sub("-.*", "", interval)) + 
                   (as.numeric(sub(".*-", "", interval)) - as.numeric(sub("-.*", "", interval))) / 2,
    significant = slope_ci_lower > 0 | slope_ci_upper < 0
  )

if (nrow(results_clean) == 0) {
  cat("Warning: No valid results to plot\n")
  quit(status = 0)
}

cat("Plotting", nrow(results_clean), "intervals with valid data\n")

# Color palette
colors <- c("TRUE" = "#d62728", "FALSE" = "#1f77b4")  # Red for significant, blue for non-significant

# Plot 1: Heritability estimates with confidence intervals
p1 <- ggplot(results_clean, aes(x = interval_mid, y = heritability, color = significant)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = h2_ci_lower, ymax = h2_ci_upper), width = 0.01) +
  scale_color_manual(values = colors, name = "Significant") +
  labs(
    title = "Heritability Estimates by GRM Interval",
    subtitle = "Points show heritability estimates (HE slope) with 95% bootstrap confidence intervals",
    x = "Genetic Relatedness (GRM interval midpoint)",
    y = "Heritability Estimate (h²)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray50"),
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5)

ggsave(paste0(output_prefix, "_heritability.png"), p1, width = 10, height = 6, dpi = 300)
cat("Saved heritability plot:", paste0(output_prefix, "_heritability.png"), "\n")

# Plot 2: R-squared values by interval
p2 <- ggplot(results_clean, aes(x = interval_mid, y = r_squared)) +
  geom_point(size = 3, color = "#2ca02c") +
  geom_line(color = "#2ca02c", alpha = 0.7) +
  labs(
    title = "Model Fit (R²) by GRM Interval", 
    subtitle = "Higher R² indicates better model fit within each interval",
    x = "Genetic Relatedness (GRM interval midpoint)",
    y = "R-squared"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray50")
  ) +
  ylim(0, 1)

ggsave(paste0(output_prefix, "_rsquared.png"), p2, width = 10, height = 6, dpi = 300)
cat("Saved R-squared plot:", paste0(output_prefix, "_rsquared.png"), "\n")

# Plot 3: Sample sizes by interval
p3 <- ggplot(results_clean, aes(x = interval_mid, y = n_pairs)) +
  geom_bar(stat = "identity", fill = "#ff7f0e", alpha = 0.7) +
  labs(
    title = "Sample Sizes by GRM Interval",
    subtitle = "Number of pairs used in each HE regression", 
    x = "Genetic Relatedness (GRM interval midpoint)",
    y = "Number of Pairs"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray50")
  )

ggsave(paste0(output_prefix, "_sample_sizes.png"), p3, width = 10, height = 6, dpi = 300)
cat("Saved sample sizes plot:", paste0(output_prefix, "_sample_sizes.png"), "\n")

# Create summary table
summary_table <- results_clean %>%
  select(interval, n_pairs, heritability, h2_ci_lower, h2_ci_upper, r_squared, significant) %>%
  arrange(interval)

cat("\n=== Heritability Summary ===\n")
print(summary_table)

# Save summary table
write_tsv(summary_table, paste0(output_prefix, "_summary.txt"))
cat("\nSaved summary table:", paste0(output_prefix, "_summary.txt"), "\n")

cat("\nHE regression analysis complete!\n")
cat("Generated files:\n")
cat("  -", paste0(output_prefix, "_heritability.png"), "\n") 
cat("  -", paste0(output_prefix, "_rsquared.png"), "\n")
cat("  -", paste0(output_prefix, "_sample_sizes.png"), "\n")
cat("  -", paste0(output_prefix, "_summary.txt"), "\n")
