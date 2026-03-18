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
  stop("Usage: Rscript bootstrap_estimators.R <crossprod_file> <grm_file>")
}
crossprod_file <- args[[1]]
grm_file        <- args[[2]]

# Column in crossprod file you want as y
y_col <- "yy"

# Bootstrap reps + CI level
B <- 100
conf_level <- 0.95
alpha <- (1 - conf_level) / 2

# -------------------------
# Helpers
# -------------------------

# Fit a single estimator model and extract the GRM slope
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
    m <- lm(y ~ GRM+ I(GRM>0.2 & GRM<0.4)+ I(GRM>0.4 & GRM<0.6), data = df)
  } else {
    # h2_unrel, h2_ped, h2_HS, h2_FS: Slope in model without intercept
    m <- lm(y ~ GRM + 0, data = df)
  }
  
  coefs <- coef(m)
  grm_slope <- coefs["GRM"]
  
  tibble(
    estimator = estimator_name,
    h2_estimate = as.numeric(grm_slope)
  )
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

bootstrap_stats_estimator <- function(df, y_col, estimator_name, B, alpha) {
  boots <- bootstraps(df, times = B, apparent = FALSE)

  boot_estimates <- boots %>%
    mutate(est_tbl = map(splits, ~ fit_estimator_model(analysis(.x), y_col = y_col, estimator_name = estimator_name))) %>%
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
  point <- fit_estimator_model(df_filtered, y_col = y_col, estimator_name = estimator_name)
  
  # Get bootstrap statistics
  bs <- bootstrap_stats_estimator(df_filtered, y_col = y_col, estimator_name = estimator_name, B = B, alpha = alpha)

  point %>%
    bind_cols(bs) %>%
    rename(point_estimate = h2_estimate) %>%
    mutate(n_pairs = nrow(df_filtered)) %>%
    select(estimator, n_pairs, point_estimate, boot_se, ci_low, ci_high)
}


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
print("Range of yy values:")
print(range(yy$yy, na.rm = TRUE))


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
print("Range of GRM values:")
print(range(x$GRM, na.rm = TRUE))

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
# Run estimators
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

# -------------------------
# Output
# -------------------------
print(results_table, n = Inf)

out_tsv <- sub("\\.[^.]*$", "", basename(crossprod_file))
out_tsv <- paste0(out_tsv, "__boot_estimators.tsv")
readr::write_tsv(results_table, out_tsv)
message(sprintf("[INFO] Wrote: %s", out_tsv))

# make a quick plot of point estimates + CIs
results_table %>%
  ggplot(aes(x = estimator, y = point_estimate, color = estimator)) +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high)) +
  theme_bw() +
  labs(x = "Estimator", y = "h2 Estimate (slope)", color = "Estimator") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  -> p

out_plot <- sub("\\.[^.]*$", "", basename(crossprod_file))
out_plot <- paste0(out_plot, "__boot_estimators_plot.png")
ggsave(out_plot, p, width = 10, height = 6)
message(sprintf("[INFO] Wrote: %s", out_plot))

