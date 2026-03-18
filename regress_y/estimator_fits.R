#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(rsample)
  library(purrr)
})

# -------------------------
# Config / inputs
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_estimator_fits.R <crossprod_file> <grm_file>")
}
crossprod_file <- args[[1]]
grm_file        <- args[[2]]

# Column in crossprod file you want as y
y_col <- "yy"

# Bootstrap reps + CI level
B <- 100
conf_level <- 0.95
alpha <- (1 - conf_level) / 2

# Bin width for plotting
bin_width <- 0.005

# -------------------------
# Helpers
# -------------------------

# Fit a single estimator model
fit_estimator_model <- function(df, y_col, estimator_name) {
  y <- df[[y_col]]

  if (estimator_name == "h2_unrel_SE") {
    # Slope in model with intercept
    m <- lm(y ~ GRM, data = df)
  } else if (estimator_name == "h2_ped_W26") {
    # Slope with quadratic term, intercept forced to zero
    m <- lm(y ~ GRM + I(GRM^2) + 0, data = df)
  } else if (estimator_name %in% c("h2_ped_SE", "h2_HS_SE", "h2_FS_SE")) {
    # Slope in model with intercept
    m <- lm(y ~ GRM, data = df)
  } else if (estimator_name %in% c("h2_ped_SE_factors")) {
    # Slope in model with intercept
    m <- lm(y ~ GRM+ I(GRM>0.1 & GRM<=0.2) + I(GRM>0.2 & GRM<=0.4)+ I(GRM>0.4 & GRM<=0.6), data = df)
  } else {
    # h2_unrel, h2_ped, h2_HS, h2_FS: Slope in model without intercept
    m <- lm(y ~ GRM + 0, data = df)
  }

  m
}

# Bootstrap stats for estimator
bootstrap_stats_estimator <- function(df, y_col, estimator_name, B, alpha) {
  boots <- bootstraps(df, times = B, apparent = FALSE)

  boot_estimates <- boots %>%
    mutate(est_tbl = map(splits, ~ {
      model <- fit_estimator_model(analysis(.x), y_col = y_col, estimator_name = estimator_name)
      coefs <- coef(model)
      grm_slope <- coefs["GRM"]
      tibble(h2_estimate = as.numeric(grm_slope))
    })) %>%
    select(est_tbl) %>%
    unnest(est_tbl)

  boot_estimates %>%
    summarise(
      boot_se  = sd(h2_estimate, na.rm = TRUE),
      ci_low   = as.numeric(quantile(h2_estimate, probs = alpha, na.rm = TRUE, type = 7)),
      ci_high  = as.numeric(quantile(h2_estimate, probs = 1 - alpha, na.rm = TRUE, type = 7)),
      .groups  = "drop"
    )
}

# Run estimator with bootstrap
run_estimator <- function(df, y_col, estimator_name, B, alpha) {
  # Filter to appropriate GRM range
  range <- get_grm_filter(estimator_name)
  df_filtered <- df %>%
    filter(GRM >= range[1], GRM <= range[2])

  if (nrow(df_filtered) == 0) {
    return(
      tibble(
        estimator = estimator_name,
        n_pairs = 0,
        point_estimate = NA_real_,
        boot_se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      )
    )
  }

  # Get point estimate
  model <- fit_estimator_model(df_filtered, y_col = y_col, estimator_name = estimator_name)
  coefs <- coef(model)
  grm_slope <- coefs["GRM"]
  point <- tibble(
    estimator = estimator_name,
    h2_estimate = as.numeric(grm_slope)
  )

  # Get bootstrap statistics
  bs <- bootstrap_stats_estimator(df_filtered, y_col = y_col, estimator_name = estimator_name, B = B, alpha = alpha)

  point %>%
    bind_cols(bs) %>%
    rename(point_estimate = h2_estimate) %>%
    mutate(n_pairs = nrow(df_filtered)) %>%
    select(estimator, n_pairs, point_estimate, boot_se, ci_low, ci_high)
}

# Get the GRM range for each estimator
get_grm_filter <- function(estimator_name) {
  ranges <- list(
    h2_unrel = c(-Inf, 0.02),
    h2_ped = c(0.05, 0.7),
    h2_ped_W26 = c(0.05, 0.7),
    h2_HS = c(0.2, 0.3),
    h2_FS = c(0.4, 0.6),
    h2_HS_SE = c(0.2, 0.3),
    h2_FS_SE = c(0.4, 0.6),
    h2_ped_SE = c(0.05, 0.7),
    h2_ped_SE_factors = c(0.05, 0.7),
    h2_unrel_SE = c(-Inf, 0.02)  # All data
  )
  ranges[[estimator_name]]
}

# Create bins
make_bins <- function(grm, breaks) {
  cut(grm, breaks = breaks, include.lowest = TRUE, right = TRUE)
}

# Predict values for an estimator across its range
predict_estimator <- function(model, estimator_name, grm_range, bin_width) {
  # Handle infinite bounds
  range_min <- ifelse(is.infinite(grm_range[1]), -0.05, grm_range[1])
  range_max <- ifelse(is.infinite(grm_range[2]), 0.70, grm_range[2])

  # Create sequence of GRM values across the estimator's range
  grm_seq <- seq(range_min, range_max, by = bin_width)

  # Create prediction data frame
  pred_data <- data.frame(GRM = grm_seq)

  # For h2_ped_SE_factors, we need the factor variables
  if (estimator_name == "h2_ped_SE_factors") {
    pred_data$`I(GRM > 0.1 & GRM <= 0.2)` <- as.numeric(pred_data$GRM > 0.1 & pred_data$GRM <= 0.2)
    pred_data$`I(GRM > 0.2 & GRM <= 0.4)` <- as.numeric(pred_data$GRM > 0.2 & pred_data$GRM <= 0.4)
    pred_data$`I(GRM > 0.4 & GRM <= 0.6)` <- as.numeric(pred_data$GRM > 0.4 & pred_data$GRM <= 0.6)
  }

  # Predict
  pred <- predict(model, newdata = pred_data)

  data.frame(
    GRM = grm_seq,
    predicted_y = pred,
    estimator = estimator_name
  )
}

# pretty labels
estimator_labels <- list(
  h2_unrel = bquote({hat(h)^2}[unrel]),
  h2_unrel_SE = bquote({hat(h)^2}[unrel*"*"]),
  
  h2_ped = bquote({hat(h)^2}[ped]),
  h2_ped_SE = bquote({hat(h)^2}[ped*"*"]),
  h2_ped_SE_factors = bquote({hat(h)^2}[ped*"*f"]),
  
  h2_ped_W26 = bquote({hat(h)^2}[ped*"W26"]),
  
  h2_HS = bquote({hat(h)^2}[HS]),
  h2_HS_SE = bquote({hat(h)^2}[HS*"*"]),
  
  h2_FS = bquote({hat(h)^2}[FS]),
  h2_FS_SE = bquote({hat(h)^2}[FS*"*"])
)
# -------------------------
# Load + join data
# -------------------------
yy <- readr::read_delim(crossprod_file, delim = "\t", show_col_types = FALSE, col_names=TRUE)
colnames(yy) <- c("IID2", "IID1","yy")
yy %>% filter(as.numeric(IID1)>0, as.numeric(IID2)>0) -> yy
yy$IID1 <- as.character(yy$IID1)
yy$IID2 <- as.character(yy$IID2)
print("First few rows of crossprod data:")
print(head(yy))
yy <- yy %>%
  filter(!is.na(yy))
print("Number of pairs in crossprod data after filtering NA yy:")
print(nrow(yy))

x <- readr::read_delim(
  grm_file,
  delim = "\t",
  col_names = TRUE,
  col_types = readr::cols(.default = readr::col_guess()),
  show_col_types = FALSE
)

print("First few rows of GRM data:")
print(head(x))

if (ncol(x) < 4) stop("GRM file must have at least 4 columns: IID1 IID2 nsnps GRM")
colnames(x)[1:4] <- c("IID1", "IID2", "N_SNPs",  "GRM")
x$IID1 <- as.character(x$IID1)
x$IID2 <- as.character(x$IID2)

xyy <- right_join(x, yy, by = c("IID1", "IID2")) %>%
  mutate(GRM = as.numeric(GRM))

if (!y_col %in% names(xyy)) {
  stop(sprintf("Could not find y column '%s' in joined table.", y_col))
}

xyy <- xyy %>%
  filter(!is.na(GRM), !is.na(.data[[y_col]]))
print("Range of GRM values after join and filtering:")
print(range(xyy$GRM, na.rm = TRUE))

# -------------------------
# Create binned data
# -------------------------
bin_breaks <- seq(-0.05, 0.70, by = bin_width)

binned_data <- xyy %>%
  mutate(GRM_bin = make_bins(GRM, bin_breaks)) %>%
  filter(!is.na(GRM_bin), IID1>0,IID2>0) %>%
  group_by(GRM_bin) %>%
  summarise(
    x = mean(GRM),
    y = mean(.data[[y_col]]),
    n = n(),
    CI_low = quantile(.data[[y_col]], 0.025, na.rm = TRUE),
    CI_high = quantile(.data[[y_col]], 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------
# Run estimators with bootstrap
# -------------------------
estimators_list <- c("h2_unrel", "h2_ped", "h2_ped_W26", "h2_HS", "h2_FS", 
                     "h2_HS_SE", "h2_FS_SE", "h2_ped_SE", "h2_unrel_SE", "h2_ped_SE_factors")

results_table <- map_df(
  estimators_list, 
  ~ run_estimator(xyy, y_col = y_col, estimator_name = .x, B = B, alpha = alpha)
) %>%
  mutate(
    crossprod_file = crossprod_file,
    grm_file = grm_file,
    boot_reps = B,
    conf_level = conf_level
  ) %>%
  select(crossprod_file, grm_file, boot_reps, conf_level, everything())
estimator_predictions <- map_df(estimators_list, function(est_name) {
  # Filter data to estimator range
  grm_range <- get_grm_filter(est_name)
  df_filtered <- xyy %>%
    filter(GRM >= grm_range[1], GRM <= grm_range[2])

  if (nrow(df_filtered) == 0) {
    return(NULL)
  }

  # Fit model
  model <- fit_estimator_model(df_filtered, y_col, est_name)

  # Predict across range
  predict_estimator(model, est_name, grm_range, bin_width)
})

# -------------------------
# Create plots
# -------------------------

# Define colors for estimators
estimator_colors <- c(
  # Unrelated (blue family)
  "h2_unrel" = "#4E79A7",
  "h2_unrel_SE" = "#1F77B4",

  # Pedigree (red family)
  "h2_ped" = "#E15759",
  "h2_ped_SE" = "#D62728",
  "h2_ped_SE_factors" = "#A52A2A",

  # Ped + quadratic (purple family)
  "h2_ped_W26" = "#B279A2",

  # Half-sibs (green family)
  "h2_HS" = "#59A14F",
  "h2_HS_SE" = "#2E6B2E",

  # Full-sibs (gold / ochre family)
  "h2_FS" = "#F1A340",
  "h2_FS_SE" = "#B36A00"
)

# Define linetypes: unadjusted (without SE) get dashed lines
estimator_linetypes <- c(
  "h2_unrel" = "dashed",
  "h2_ped" = "dashed",
  "h2_ped_W26" = "dashed",
  "h2_HS" = "dashed",
  "h2_FS" = "dashed",
  "h2_HS_SE" = "solid",
  "h2_FS_SE" = "solid",
  "h2_ped_SE" = "solid",
  "h2_unrel_SE" = "solid",
  "h2_ped_SE_factors" = "solid"
)

# Plot 1: Binned data with estimator fits
p1 <- ggplot() +
  # Binned data points with error bars
  geom_point(data = binned_data %>% filter(n > 50),
                  aes(x = x, y = y),
                  alpha = 0.6, size = 0.9, col="black") +
  # Regression lines for each estimator
  geom_line(data = estimator_predictions,
            aes(x = GRM, y = predicted_y, color = estimator, linetype = estimator),
            linewidth = 0.6, alpha=0.8) +
  # Customize colors and linetypes
  scale_color_manual(
  values = estimator_colors,
  labels = estimator_labels
    ) +
    scale_linetype_manual(
    values = estimator_linetypes,
    labels = estimator_labels
    ) +
  # Theme and labels
  theme_bw() +
  labs(x = "Genetic Covariance", y = "yy (cross-product)",
       color = "Estimator",
       linetype = "Estimator",
       title = "Binned Data with Estimator Fits",
       subtitle = paste("Bin width:", bin_width)) +
  lims(x=c(-0.05,0.6)) +
  # Legend position
  theme(legend.position = "right")

# Plot 2: Estimator point estimates with confidence intervals
p2 <- results_table %>%
  ggplot(aes(x = estimator, y = point_estimate, color = estimator)) +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high, linetype = estimator)) +
  scale_color_manual(values = estimator_colors, labels = estimator_labels) +
  scale_linetype_manual(values = estimator_linetypes, labels = estimator_labels) +
  scale_x_discrete(labels = estimator_labels) +
  theme_bw() +
  geom_hline(aes(yintercept = 0), color = "gray") +
  labs(x = "Estimator", y = expression(hat(h)^2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  coord_flip()

# -------------------------
# Save plots
# -------------------------
out_plot1 <- sub("\\.[^.]*$", "", basename(crossprod_file))
out_plot1 <- paste0(out_plot1, "__estimator_fits_plot.pdf")
ggsave(out_plot1, p1, width = 7, height = 4)
message(sprintf("[INFO] Wrote: %s", out_plot1))

out_plot2 <- sub("\\.[^.]*$", "", basename(crossprod_file))
out_plot2 <- paste0(out_plot2, "__estimator_comparison_plot.pdf")
ggsave(out_plot2, p2, width = 3, height = 4)
message(sprintf("[INFO] Wrote: %s", out_plot2))

# Save results table
out_tsv <- sub("\\.[^.]*$", "", basename(crossprod_file))
out_tsv <- paste0(out_tsv, "__estimator_comparison.tsv")
readr::write_tsv(results_table, out_tsv)
message(sprintf("[INFO] Wrote: %s", out_tsv))

# Also save the binned data
# out_bin_tsv <- sub("\\.[^.]*$", "", basename(crossprod_file))
# out_bin_tsv <- paste0(out_bin_tsv, "__binned_data.tsv")
# readr::write_tsv(binned_data, out_bin_tsv)
# message(sprintf("[INFO] Wrote: %s", out_bin_tsv))

# And the predictions
# out_pred_tsv <- sub("\\.[^.]*$", "", basename(crossprod_file))
# out_pred_tsv <- paste0(out_pred_tsv, "__estimator_predictions.tsv")
# readr::write_tsv(estimator_predictions, out_pred_tsv)
# message(sprintf("[INFO] Wrote: %s", out_pred_tsv))