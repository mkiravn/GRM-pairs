#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "fit_bin_regression_jackknife.R <full_bins.tsv> <jackknife_detail.tsv> <out.tsv>",
  "[--formula 'avg_pheno_crossprod ~ avg_grm']",
  "[--weights n_pairs|uniform]",
  "[--min-grm -Inf] [--max-grm Inf]",
  "[--range low:high] [--range low:high ...]",
  "[--cov-out cov.tsv]",
  sep = "\n  "
)

if (length(args) < 3) {
  stop(usage)
}

full_file <- args[1]
detail_file <- args[2]
out_file <- args[3]
formula_text <- "avg_pheno_crossprod ~ avg_grm"
weight_mode <- "n_pairs"
min_grm <- -Inf
max_grm <- Inf
range_specs <- character(0)
cov_out <- NA_character_

parse_range <- function(x) {
  parts <- strsplit(x, ":", fixed = TRUE)[[1]]
  if (length(parts) != 2) {
    stop("Range must have form low:high, got: ", x)
  }
  out <- as.numeric(parts)
  if (any(!is.finite(out))) {
    stop("Range bounds must be numeric, got: ", x)
  }
  if (out[1] > out[2]) {
    stop("Range lower bound exceeds upper bound: ", x)
  }
  out
}

i <- 4L
while (i <= length(args)) {
  key <- args[i]
  if (i == length(args)) {
    stop("Missing value for ", key, "\n", usage)
  }
  value <- args[i + 1L]

  if (key == "--formula") {
    formula_text <- value
  } else if (key == "--weights") {
    weight_mode <- value
  } else if (key == "--min-grm") {
    min_grm <- as.numeric(value)
  } else if (key == "--max-grm") {
    max_grm <- as.numeric(value)
  } else if (key == "--range") {
    range_specs <- c(range_specs, value)
  } else if (key == "--cov-out") {
    cov_out <- value
  } else {
    stop("Unknown argument: ", key, "\n", usage)
  }
  i <- i + 2L
}

if (!(weight_mode %in% c("n_pairs", "uniform"))) {
  stop("--weights must be one of: n_pairs, uniform")
}

formula_obj <- as.formula(formula_text)
vars_needed <- all.vars(formula_obj)
if (!("avg_pheno_crossprod" %in% vars_needed)) {
  stop("Formula must include avg_pheno_crossprod on the left-hand side.")
}

full_df <- read.table(full_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
detail_df <- read.table(detail_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

if (length(range_specs)) {
  parsed_ranges <- do.call(
    rbind,
    lapply(range_specs, parse_range)
  )
} else {
  parsed_ranges <- matrix(c(min_grm, max_grm), ncol = 2, byrow = TRUE)
  range_specs <- sprintf("%g:%g", min_grm, max_grm)
}

keep_by_range <- function(x) {
  keep <- rep(FALSE, length(x))
  for (r in seq_len(nrow(parsed_ranges))) {
    keep <- keep | (x >= parsed_ranges[r, 1] & x <= parsed_ranges[r, 2])
  }
  keep
}

keep_full <- keep_by_range(full_df$avg_grm) & full_df$n_pairs > 0
keep_detail <- keep_by_range(detail_df$avg_grm) & detail_df$n_pairs > 0

full_df <- full_df[keep_full, ]
detail_df <- detail_df[keep_detail, ]

if (!nrow(full_df) || !nrow(detail_df)) {
  stop("No bins remain after filtering.")
}

weight_vec <- function(df) {
  if (weight_mode == "n_pairs") {
    df$n_pairs
  } else {
    rep(1, nrow(df))
  }
}

fit_block <- function(df) {
  wts <- weight_vec(df)
  mf <- model.frame(formula_obj, data = df, na.action = na.fail)
  y <- model.response(mf)
  x <- model.matrix(formula_obj, data = mf)
  fit <- lm.wfit(x = x, y = y, w = wts)
  out <- fit$coefficients
  names(out) <- colnames(x)
  out
}

full_coef <- fit_block(full_df)
coef_names <- names(full_coef)

blocks <- sort(unique(detail_df$block))
loo_mat <- matrix(
  NA_real_,
  nrow = length(blocks),
  ncol = length(coef_names),
  dimnames = list(as.character(blocks), coef_names)
)

for (idx in seq_along(blocks)) {
  b <- blocks[idx]
  coef_b <- fit_block(detail_df[detail_df$block == b, ])
  if (!identical(names(coef_b), coef_names)) {
    stop("Coefficient mismatch in block ", b, ". Formula may be rank-deficient on a leave-one-block-out fit.")
  }
  loo_mat[idx, ] <- coef_b
}

loo_mean <- colMeans(loo_mat)
centered <- sweep(loo_mat, 2L, loo_mean, FUN = "-")
jk_scale <- (nrow(loo_mat) - 1) / nrow(loo_mat)
jk_cov <- jk_scale * crossprod(centered)
jk_se <- sqrt(diag(jk_cov))

out <- data.frame(
  parameter = coef_names,
  estimate = unname(full_coef),
  jackknife_mean = unname(loo_mean[coef_names]),
  jackknife_se = unname(jk_se[coef_names]),
  n_blocks = nrow(loo_mat),
  weight_mode = weight_mode,
  formula = formula_text,
  range_spec = paste(range_specs, collapse = ","),
  n_bins_used = nrow(full_df),
  stringsAsFactors = FALSE
)

write.table(out, file = out_file, quote = FALSE, sep = "\t", row.names = FALSE)

if (!is.na(cov_out)) {
  cov_df <- expand.grid(
    parameter_1 = coef_names,
    parameter_2 = coef_names,
    stringsAsFactors = FALSE
  )
  cov_df$jackknife_cov <- as.vector(jk_cov)
  cov_df$n_blocks <- nrow(loo_mat)
  cov_df$formula <- formula_text
  cov_df$range_spec <- paste(range_specs, collapse = ",")
  cov_df$weight_mode <- weight_mode
  write.table(cov_df, file = cov_out, quote = FALSE, sep = "\t", row.names = FALSE)
}

cat("Wrote regression coefficients to", out_file, "\n")
if (!is.na(cov_out)) {
  cat("Wrote jackknife covariance to", cov_out, "\n")
}
