#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: prepare_pheno_for_grm.R <pheno.txt> <grm_prefix.grm.id> <out.txt>")
}

pheno_file <- args[1]
grm_id_file <- args[2]
out_file <- args[3]

read_pheno <- function(path) {
  dat <- read.table(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "-999")
  )

  if (ncol(dat) == 2) {
    colnames(dat) <- c("id", "pheno")
    dat$key <- dat$id
  } else if (ncol(dat) >= 3) {
    colnames(dat)[1:3] <- c("fid", "iid", "pheno")
    dat$key <- paste(dat$fid, dat$iid, sep = ":")
  } else {
    stop("Phenotype file must have either 2 columns (ID pheno) or 3 columns (FID IID pheno).")
  }

  dat
}

pheno <- read_pheno(pheno_file)
grm <- read.table(grm_id_file, header = FALSE, stringsAsFactors = FALSE)
if (ncol(grm) < 2) stop("GRM ID file must have at least 2 columns")
colnames(grm)[1:2] <- c("fid", "iid")
grm$key <- paste(grm$fid, grm$iid, sep = ":")

# Prefer exact FID:IID match, but allow fallback to IID if phenotype file only has one ID column
if (!("fid" %in% names(pheno))) {
  match_idx <- match(grm$iid, pheno$key)
} else {
  match_idx <- match(grm$key, pheno$key)
}

aligned <- data.frame(
  id = grm$key,
  pheno = pheno$pheno[match_idx],
  stringsAsFactors = FALSE
)

write.table(
  aligned,
  file = out_file,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  na = "NA"
)

cat("Wrote aligned phenotype file:", out_file, "\n")
cat("Matched rows:", sum(!is.na(aligned$pheno)), "out of", nrow(aligned), "\n")