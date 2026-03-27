#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "fit_bin_regression_jackknife.R <full_bins.tsv> <jackknife_detail.tsv> <out.tsv>",
  "[--min-grm -Inf] [--max-grm Inf] [--weights n_pairs|uniform] [--no-intercept]",
  sep = "\n  "
)

if (length(args) < 3) {
  stop(usage)
}

full_file <- args[1]
detail_file <- args[2]
out_file <- args[3]
min_grm <- -Inf
max_grm <- Inf
weight_mode <- "n_pairs"
use_intercept <- TRUE

i <- 4L
while (i <= length(args)) {
  key <- args[i]
  if (key == "--no-intercept") {
    use_intercept <- FALSE
    i <- i + 1L
    next
  }
  if (i == length(args)) {
    stop("Missing value for ", key, "\n", usage)
  }
  value <- args[i + 1L]

  if (key == "--min-grm") {
    min_grm <- as.numeric(value)
  } else if (key == "--max-grm") {
    max_grm <- as.numeric(value)
  } else if (key == "--weights") {
    weight_mode <- value
  } else {
    stop("Unknown argument: ", key, "\n", usage)
  }
  i <- i + 2L
}

if (!(weight_mode %in% c("n_pairs", "uniform"))) {
  stop("--weights must be one of: n_pairs, uniform")
}

full_df <- read.table(full_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
detail_df <- read.table(detail_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

keep_full <- full_df$avg_grm >= min_grm & full_df$avg_grm <= max_grm & full_df$n_pairs > 0
keep_detail <- detail_df$avg_grm >= min_grm & detail_df$avg_grm <= max_grm & detail_df$n_pairs > 0

full_df <- full_df[keep_full, ]
detail_df <- detail_df[keep_detail, ]

if (!nrow(full_df) || !nrow(detail_df)) {
  stop("No bins remain after filtering.")
}

fit_block <- function(df) {
  weights <- if (weight_mode == "n_pairs") df$n_pairs else rep(1, nrow(df))
  formula <- if (use_intercept) avg_pheno_crossprod ~ avg_grm else avg_pheno_crossprod ~ 0 + avg_grm
  fit <- lm(formula, data = df, weights = weights)
  coefs <- coef(fit)
  out <- list(
    intercept = if (use_intercept) unname(coefs["(Intercept)"]) else 0,
    slope = unname(coefs["avg_grm"])
  )
  out
}

full_fit <- fit_block(full_df)
blocks <- sort(unique(detail_df$block))
loo <- lapply(blocks, function(b) fit_block(detail_df[detail_df$block == b, ]))

extract_coef <- function(name) {
  vals <- vapply(loo, function(x) x[[name]], numeric(1))
  mean_loo <- mean(vals)
  se <- sqrt((length(vals) - 1) / length(vals) * sum((vals - mean_loo)^2))
  c(full = full_fit[[name]], mean_loo = mean_loo, se = se)
}

intercept_stats <- extract_coef("intercept")
slope_stats <- extract_coef("slope")

out <- data.frame(
  parameter = c("intercept", "slope"),
  estimate = c(intercept_stats["full"], slope_stats["full"]),
  jackknife_mean = c(intercept_stats["mean_loo"], slope_stats["mean_loo"]),
  jackknife_se = c(intercept_stats["se"], slope_stats["se"]),
  n_blocks = length(blocks),
  min_grm = min_grm,
  max_grm = max_grm,
  weight_mode = weight_mode,
  use_intercept = use_intercept,
  stringsAsFactors = FALSE
)

write.table(out, file = out_file, quote = FALSE, sep = "\t", row.names = FALSE)
cat("Wrote regression coefficients to", out_file, "\n")
