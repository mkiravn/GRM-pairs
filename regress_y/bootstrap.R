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
  stop("Usage: Rscript grm_crossprod_boot.R <crossprod_file> <grm_file>")
}
crossprod_file <- args[[1]]
grm_file        <- args[[2]]

# Column in crossprod file you want as y
y_col <- "yy"

# Your requested bins
bin_breaks <- c(-0.05, 0.03, 0.05, 0.10, 0.20, 0.35, 0.70)

# Bootstrap reps + CI level
B <- 100
conf_level <- 0.95
alpha <- (1 - conf_level) / 2

# -------------------------
# Helpers
# -------------------------

make_bins <- function(grm, breaks) {
  cut(grm, breaks = breaks, include.lowest = TRUE, right = TRUE)
}

fit_models <- function(df, y_col) {
  y <- df[[y_col]]

  # with intercept
  m_int <- lm(y ~ GRM, data = df)

  # intercept forced to zero
  m0 <- lm(y ~ GRM + 0, data = df)

  tibble(
    model = c("with_intercept", "intercept_zero"),
    fit   = list(m_int, m0)
  ) %>%
    mutate(coefs = map(fit, coef)) %>%
    transmute(
      model,
      term = map(coefs, names),
      estimate = map(coefs, unname)
    ) %>%
    unnest(c(term, estimate))
}

bootstrap_stats <- function(df, y_col, B, alpha) {
  boots <- bootstraps(df, times = B, apparent = FALSE)

  boot_estimates <- boots %>%
    mutate(est_tbl = map(splits, ~ fit_models(analysis(.x), y_col = y_col))) %>%
    select(est_tbl) %>%
    unnest(est_tbl)

  boot_estimates %>%
    group_by(model, term) %>%
    summarise(
      boot_se  = sd(estimate, na.rm = TRUE),
      ci_low   = as.numeric(quantile(estimate, probs = alpha, na.rm = TRUE, type = 7)),
      ci_high  = as.numeric(quantile(estimate, probs = 1 - alpha, na.rm = TRUE, type = 7)),
      .groups  = "drop"
    )
}

run_group <- function(df, group_label, y_col, B, alpha) {
  point <- fit_models(df, y_col = y_col) %>%
    rename(point_estimate = estimate)

  bs <- bootstrap_stats(df, y_col = y_col, B = B, alpha = alpha)

  point %>%
    left_join(bs, by = c("model", "term")) %>%
    mutate(group = group_label, n_pairs = nrow(df)) %>%
    select(group, n_pairs, model, term, point_estimate, boot_se, ci_low, ci_high) %>%
    arrange(group, model, term)
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
# Build groups: All + bins
# -------------------------
xyy_bins <- xyy %>%
  mutate(GRM_bin = make_bins(GRM, bin_breaks)) %>%
  filter(!is.na(GRM_bin))

# -------------------------
# Run analyses
# -------------------------
out_all <- run_group(xyy, "All", y_col = y_col, B = B, alpha = alpha)

out_bins <- xyy_bins %>%
  group_by(GRM_bin) %>%
  group_modify(~ run_group(.x, as.character(.y$GRM_bin), y_col = y_col, B = B, alpha = alpha)) %>%
  ungroup()

results_table <- bind_rows(out_all, out_bins) %>%
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
out_tsv <- paste0(out_tsv, "__boot_table.tsv")
readr::write_tsv(results_table, out_tsv)
message(sprintf("[INFO] Wrote: %s", out_tsv))

# make a quick plot of point estimates + CIs
results_table %>%
  ggplot(aes(x = group, y = point_estimate, color = term)) +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.5)) +
  theme_bw() +
  facet_wrap(~ model, scales = "free_y") +
  labs(x = "Group", y = "Point Estimate (slope)", color = "Model") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  -> p

out_plot <- sub("\\.[^.]*$", "", basename(crossprod_file))
out_plot <- paste0(out_plot, "__boot_plot.png")
ggsave(out_plot, p, width = 8, height = 6)
message(sprintf("[INFO] Wrote: %s", out_plot))


# make a quick bin plot of GRM vs yy with the regression lines, bins should be 0.005
xyy %>%
    group_by(GRM_bin = make_bins(GRM, seq(-0.05, 0.70, by = 0.005))) %>%
    filter(!is.na(GRM_bin)) %>%
    reframe(x = mean(GRM), y = mean(.data[[y_col]]), n = n(), CI_low = quantile(.data[[y_col]], 0.025, na.rm = TRUE), CI_high = quantile(.data[[y_col]], 0.975, na.rm = TRUE)) -> binned_data

# write out the binned data for plotting
out_bin_tsv <- sub("\\.[^.]*$", "", basename(crossprod_file))
out_bin_tsv <- paste0(out_bin_tsv, "__binned_data.tsv")
readr::write_tsv(binned_data, out_bin_tsv)
message(sprintf("[INFO] Wrote: %s", out_bin_tsv))

binned_data %>%
  filter(n>50) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(alpha = 0.3) +
    theme_bw() +
    lims(y=c(-0.5, 1.5)) +
    labs(x = "GRM", y= "yy") -> p_bin

out_bin_plot <- sub("\\.[^.]*$", "", basename(crossprod_file))
out_bin_plot <- paste0(out_bin_plot, "__bin_plot.png")
ggsave(out_bin_plot, p_bin, width = 8, height = 6)
message(sprintf("[INFO] Wrote: %s", out_bin_plot))

