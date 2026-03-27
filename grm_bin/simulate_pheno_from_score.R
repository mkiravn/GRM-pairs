#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "simulate_pheno_from_score.R <score.sscore> <grm_prefix> <out_pheno>",
  "[--h2 0.4] [--seed 1] [--genetic-out genetic.txt]",
  sep = "\n  "
)

if (length(args) < 3) {
  stop(usage)
}

score_file <- args[1]
grm_prefix <- args[2]
out_pheno <- args[3]
h2 <- 0.4
seed <- 1L
genetic_out <- NA_character_

i <- 4L
while (i <= length(args)) {
  key <- args[i]
  if (i == length(args)) {
    stop("Missing value for ", key, "\n", usage)
  }
  value <- args[i + 1L]

  if (key == "--h2") {
    h2 <- as.numeric(value)
  } else if (key == "--seed") {
    seed <- as.integer(value)
  } else if (key == "--genetic-out") {
    genetic_out <- value
  } else {
    stop("Unknown argument: ", key, "\n", usage)
  }
  i <- i + 2L
}

if (!is.finite(h2) || h2 < 0 || h2 > 1) {
  stop("--h2 must be between 0 and 1.")
}

score <- read.table(
  score_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  comment.char = ""
)
names(score)[names(score) == "#FID"] <- "FID"
names(score)[names(score) == "X.FID"] <- "FID"
grm_ids <- read.table(
  paste0(grm_prefix, ".grm.id"),
  header = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("FID", "IID") %in% names(score))) {
  stop("Score file must contain FID and IID columns: ", score_file)
}

score_col <- names(score)[grepl("^SCORE", names(score))]
if (!length(score_col) && "BETA_SUM" %in% names(score)) {
  score_col <- "BETA_SUM"
}
if (!length(score_col)) {
  stop("Could not find SCORE column in ", score_file)
}

colnames(grm_ids)[1:2] <- c("FID", "IID")
match_idx <- match(
  paste(grm_ids$FID, grm_ids$IID, sep = ":"),
  paste(score$FID, score$IID, sep = ":")
)

if (anyNA(match_idx)) {
  stop("Some GRM IDs are missing from the score file.")
}

g_raw <- score[[score_col[1]]][match_idx]
g_centered <- g_raw - mean(g_raw)

set.seed(seed)
e_raw <- rnorm(length(g_centered), sd = sqrt(1 - h2))
y <- g_centered + e_raw
y <- y - mean(y)

aligned <- data.frame(
  id = paste(grm_ids$FID, grm_ids$IID, sep = ":"),
  pheno = y,
  stringsAsFactors = FALSE
)

write.table(
  aligned,
  file = out_pheno,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

if (!is.na(genetic_out)) {
  genetic_aligned <- data.frame(
    id = paste(grm_ids$FID, grm_ids$IID, sep = ":"),
    genetic_component = g_centered,
    stringsAsFactors = FALSE
  )
  write.table(
    genetic_aligned,
    file = genetic_out,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
}

cat("Wrote phenotype to", out_pheno, "\n")
cat("Score column:", score_col[1], "\n")
