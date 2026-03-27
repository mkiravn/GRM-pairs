#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: summarize_simulations.R <run_dir> <out_tsv>")
}

run_dir <- args[1]
out_tsv <- args[2]

list_rep_files <- function(pattern) {
  files <- Sys.glob(file.path(run_dir, "rep_*", pattern))
  files[file.info(files)$isdir %in% FALSE]
}

read_table_strict <- function(path, add_replicate = FALSE) {
  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (add_replicate && !("replicate" %in% names(df))) {
    rep_dir <- basename(dirname(path))
    rep_id <- as.integer(sub("^rep_", "", rep_dir))
    df$replicate <- rep_id
  }
  df
}

binned_files <- list_rep_files("binned.tsv")
jk_files <- list_rep_files("jackknife.tsv")
truth_files <- list_rep_files("truth.tsv")

if (!length(binned_files)) {
  stop("No replicate outputs found under ", run_dir)
}
if (length(jk_files) != length(binned_files)) {
  stop("Expected one jackknife.tsv per replicate in ", run_dir)
}
if (length(truth_files) != length(binned_files)) {
  stop("Expected one truth.tsv per replicate in ", run_dir)
}

binned <- do.call(rbind, lapply(binned_files, read_table_strict, add_replicate = TRUE))
jackknife <- do.call(rbind, lapply(jk_files, read_table_strict, add_replicate = TRUE))
truth <- do.call(rbind, lapply(truth_files, read_table_strict))

rename_with_prefix <- function(df, mapping) {
  for (nm in names(mapping)) {
    names(df)[names(df) == nm] <- mapping[[nm]]
  }
  df
}

binned <- rename_with_prefix(
  binned,
  c(avg_pheno_crossprod = "estimate", n_pairs = "n_pairs_obs")
)
jackknife <- rename_with_prefix(
  jackknife,
  c(avg_pheno_crossprod = "estimate_jk", n_pairs = "n_pairs_jk", se = "jackknife_se")
)

merged <- merge(
  binned[, c("bin", "lower", "upper", "avg_grm", "estimate", "n_pairs_obs", "replicate")],
  jackknife[, c("bin", "replicate", "estimate_jk", "n_pairs_jk", "jackknife_se")],
  by = c("bin", "replicate"),
  all = TRUE
)

merged <- merge(
  merged,
  truth[, c("bin", "replicate", "expected_avg_pheno_crossprod", "oracle_se_indep", "n_pairs")],
  by = c("bin", "replicate"),
  all = TRUE
)

summarize_bin <- function(df) {
  est <- df$estimate
  truth_vals <- df$expected_avg_pheno_crossprod
  jackknife_se <- df$jackknife_se
  avg_grm <- df$avg_grm

  data.frame(
    bin = df$bin[1],
    lower = df$lower[1],
    upper = df$upper[1],
    n_reps = sum(!is.na(est)),
    n_pairs = round(mean(df$n_pairs, na.rm = TRUE)),
    avg_grm = mean(avg_grm, na.rm = TRUE),
    expected_avg_pheno_crossprod = mean(truth_vals, na.rm = TRUE),
    mean_estimate = mean(est, na.rm = TRUE),
    empirical_sd = if (sum(!is.na(est)) > 1) sd(est, na.rm = TRUE) else NA_real_,
    mean_bias = mean(est - truth_vals, na.rm = TRUE),
    rmse = sqrt(mean((est - truth_vals)^2, na.rm = TRUE)),
    mean_jackknife_se = mean(jackknife_se, na.rm = TRUE),
    mean_oracle_se_indep = mean(df$oracle_se_indep, na.rm = TRUE),
    jackknife_to_empirical = mean(jackknife_se, na.rm = TRUE) /
      sd(est, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

split_bins <- split(merged, merged$bin)
summary_df <- do.call(rbind, lapply(split_bins, summarize_bin))
summary_df <- summary_df[order(summary_df$bin), ]

write.table(
  summary_df,
  file = out_tsv,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

cat("Wrote simulation summary to", out_tsv, "\n")
