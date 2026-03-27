#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "sample_effect_sizes_from_bim.R <bim> <out_tsv>",
  "[--h2 0.4] [--pi 1.0] [--seed 1]",
  sep = "\n  "
)

if (length(args) < 2) {
  stop(usage)
}

bim_file <- args[1]
out_file <- args[2]
h2 <- 0.4
pi_causal <- 1.0
seed <- 1L

i <- 3L
while (i <= length(args)) {
  key <- args[i]
  if (i == length(args)) {
    stop("Missing value for ", key, "\n", usage)
  }
  value <- args[i + 1L]

  if (key == "--h2") {
    h2 <- as.numeric(value)
  } else if (key == "--pi") {
    pi_causal <- as.numeric(value)
  } else if (key == "--seed") {
    seed <- as.integer(value)
  } else {
    stop("Unknown argument: ", key, "\n", usage)
  }
  i <- i + 2L
}

if (!is.finite(h2) || h2 < 0 || h2 > 1) {
  stop("--h2 must be between 0 and 1.")
}
if (!is.finite(pi_causal) || pi_causal <= 0 || pi_causal > 1) {
  stop("--pi must be in (0, 1].")
}

set.seed(seed)

bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
if (ncol(bim) < 6) {
  stop("Expected PLINK .bim with 6 columns: ", bim_file)
}

colnames(bim)[1:6] <- c("CHR", "ID", "CM", "BP", "A1", "A2")
m <- nrow(bim)
causal <- runif(m) < pi_causal

if (!any(causal)) {
  causal[sample.int(m, 1L)] <- TRUE
}

beta <- numeric(m)
beta_sd <- sqrt(h2 / sum(causal))
beta[causal] <- rnorm(sum(causal), sd = beta_sd)

effects <- data.frame(
  ID = bim$ID,
  A1 = bim$A1,
  BETA = beta,
  CAUSAL = as.integer(causal),
  stringsAsFactors = FALSE
)

write.table(
  effects,
  file = out_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

cat("Wrote effect sizes to", out_file, "\n")
cat("Variants:", m, " Causal:", sum(causal), " Beta SD:", beta_sd, "\n")
