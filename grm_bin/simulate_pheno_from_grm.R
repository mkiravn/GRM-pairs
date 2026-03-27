#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "simulate_pheno_from_grm.R <grm_prefix> <out_pheno>",
  "[--h2 0.4] [--seed 1] [--replicate 1]",
  "[--truth-out truth.tsv] [--meta-out meta.tsv]",
  sep = "\n  "
)

parse_args <- function(argv) {
  if (length(argv) < 2) {
    stop(usage)
  }

  opts <- list(
    grm_prefix = argv[1],
    out_pheno = argv[2],
    h2 = 0.4,
    seed = 1L,
    replicate = 1L,
    truth_out = NA_character_,
    meta_out = NA_character_
  )

  i <- 3L
  while (i <= length(argv)) {
    key <- argv[i]
    if (i == length(argv)) {
      stop("Missing value for argument: ", key, "\n", usage)
    }
    value <- argv[i + 1L]

    if (key == "--h2") {
      opts$h2 <- as.numeric(value)
    } else if (key == "--seed") {
      opts$seed <- as.integer(value)
    } else if (key == "--replicate") {
      opts$replicate <- as.integer(value)
    } else if (key == "--truth-out") {
      opts$truth_out <- value
    } else if (key == "--meta-out") {
      opts$meta_out <- value
    } else {
      stop("Unknown argument: ", key, "\n", usage)
    }
    i <- i + 2L
  }

  if (!is.finite(opts$h2) || opts$h2 < 0 || opts$h2 > 1) {
    stop("--h2 must be between 0 and 1.")
  }

  opts
}

make_bins <- function(min_edge = -0.3,
                      split1 = -0.02,
                      split2 = 0.02,
                      max_edge = 1.05,
                      small_width = 0.001,
                      large_width = 0.005) {
  edges <- c(min_edge, split1)

  x <- split1 + small_width
  while (x <= split2 + 1e-12) {
    edges <- c(edges, x)
    x <- x + small_width
  }

  x <- split2 + large_width
  while (x <= max_edge + 1e-12) {
    edges <- c(edges, x)
    x <- x + large_width
  }

  if (abs(tail(edges, 1) - max_edge) > 1e-12) {
    edges <- c(edges, max_edge)
  }
  edges
}

read_grm_matrix <- function(prefix) {
  id_path <- paste0(prefix, ".grm.id")
  bin_path <- paste0(prefix, ".grm.bin")

  ids <- read.table(id_path, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(ids) < 2) {
    stop("GRM ID file must have at least two columns: ", id_path)
  }
  colnames(ids)[1:2] <- c("FID", "IID")

  n <- nrow(ids)
  expected_values <- n * (n + 1) / 2

  con <- file(bin_path, "rb")
  on.exit(close(con), add = TRUE)
  vals <- readBin(con, what = numeric(), n = expected_values, size = 4L)
  if (length(vals) != expected_values) {
    stop("Expected ", expected_values, " float entries in ", bin_path,
         " but found ", length(vals), ".")
  }

  mat <- matrix(0, nrow = n, ncol = n)
  idx <- 1L
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      mat[i, j] <- vals[idx]
      mat[j, i] <- vals[idx]
      idx <- idx + 1L
    }
  }

  list(ids = ids, grm = mat)
}

compute_truth_by_bin <- function(grm, h2, replicate_id) {
  n <- nrow(grm)
  edges <- make_bins()
  lower <- edges[-length(edges)]
  upper <- edges[-1L]

  vals <- grm[lower.tri(grm)]
  bins <- findInterval(vals, edges, rightmost.closed = TRUE, all.inside = FALSE)
  valid <- bins >= 1L & bins <= length(lower)
  vals <- vals[valid]
  bins <- bins[valid]

  out <- data.frame(
    replicate = replicate_id,
    bin = seq_along(lower),
    lower = lower,
    upper = upper,
    avg_grm = NA_real_,
    expected_avg_pheno_crossprod = NA_real_,
    n_pairs = 0L,
    oracle_se_indep = NA_real_,
    stringsAsFactors = FALSE
  )

  if (!length(vals)) {
    return(out[out$n_pairs > 0, , drop = FALSE])
  }

  split_vals <- split(vals, bins)
  for (name in names(split_vals)) {
    k <- as.integer(name)
    grm_vals <- split_vals[[name]]
    rho_vals <- h2 * grm_vals
    m <- length(grm_vals)

    out$avg_grm[k] <- mean(grm_vals)
    out$expected_avg_pheno_crossprod[k] <- mean(rho_vals)
    out$n_pairs[k] <- m
    out$oracle_se_indep[k] <- sqrt(sum(1 + rho_vals^2)) / m
  }

  out[out$n_pairs > 0, , drop = FALSE]
}

opts <- parse_args(args)
set.seed(opts$seed)

grm_obj <- read_grm_matrix(opts$grm_prefix)
ids <- grm_obj$ids
grm <- (grm_obj$grm + t(grm_obj$grm)) / 2
n <- nrow(grm)

eig <- eigen(grm, symmetric = TRUE)
eig$values[eig$values < 0] <- 0

z_g <- rnorm(n)
g_component <- eig$vectors %*% (sqrt(eig$values) * z_g)
g_component <- as.numeric(g_component)
e_component <- rnorm(n)

diag_mean <- mean(diag(grm))
if (!is.finite(diag_mean) || diag_mean <= 0) {
  stop("Mean diagonal of GRM must be positive.")
}

g_component <- g_component / sqrt(diag_mean)
pheno <- sqrt(opts$h2) * g_component + sqrt(1 - opts$h2) * e_component

pheno_out <- data.frame(
  id = paste(ids$FID, ids$IID, sep = ":"),
  pheno = pheno,
  stringsAsFactors = FALSE
)

write.table(
  pheno_out,
  file = opts$out_pheno,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

if (!is.na(opts$truth_out)) {
  truth <- compute_truth_by_bin(grm, opts$h2, opts$replicate)
  write.table(
    truth,
    file = opts$truth_out,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}

if (!is.na(opts$meta_out)) {
  meta <- data.frame(
    replicate = opts$replicate,
    n = n,
    h2 = opts$h2,
    seed = opts$seed,
    mean_diag = diag_mean,
    min_eigenvalue = min(eig$values),
    max_eigenvalue = max(eig$values),
    stringsAsFactors = FALSE
  )
  write.table(
    meta,
    file = opts$meta_out,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
}

cat("Wrote phenotype to", opts$out_pheno, "\n")
cat("n =", n, " h2 =", opts$h2, " seed =", opts$seed, "\n")
