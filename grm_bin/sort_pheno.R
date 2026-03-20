#!/usr/bin/env Rscript

# Fast phenotype sorting using R (data.table for efficient joins)
# Usage: Rscript sort_pheno.R <grm_prefix> <pheno_input> <pheno_output>

library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript sort_pheno.R <grm_prefix> <pheno_input> <pheno_output>\n")
  cat("Sorts phenotype file to match GRM individual order\n")
  quit(status = 1)
}

grm_prefix <- args[1]
pheno_in <- args[2]
pheno_out <- args[3]

# Read GRM ID file
grm_id_file <- paste0(grm_prefix, ".grm.id")
cat("Reading GRM ID file:", grm_id_file, "\n")

grm_ids <- fread(grm_id_file, header = FALSE, col.names = c("FID", "IID"))
n_indiv <- nrow(grm_ids)
cat("Found", n_indiv, "individuals in GRM ID file\n")

# Create ID columns for matching (ensure character type)
grm_ids[, iid_only := as.character(IID)]
grm_ids[, has_valid_iid := IID >= 0]

# Read phenotype file
cat("Reading phenotype file:", pheno_in, "\n")
pheno <- fread(pheno_in, header = TRUE, select = 1:2)
colnames(pheno) <- c("id", "pheno")
pheno[, id := as.character(id)]

cat("Found", nrow(pheno), "samples in phenotype file\n")

# Attempting ID matching: try IID-only format
cat("Attempting ID matching...\n")

# Join on IID (phenotype id column matched to GRM iid_only)
# Add sequence to preserve original GRM order
grm_ids[, seq_grm := .I]

result <- merge(grm_ids[, .(FID, IID, iid_only, has_valid_iid, seq_grm)], 
                pheno[, .(id, pheno)], 
                by.x = "iid_only", by.y = "id", 
                all.x = TRUE)

# Order by original GRM order
result <- result[order(seq_grm)]

# Count matches
n_matched <- sum(!is.na(result$pheno))
n_missing_pheno <- sum(is.na(result$pheno))
n_negative_iid <- sum(!result$has_valid_iid)

cat("Matched", n_matched, "individuals with phenotype data\n")
cat("Found", n_missing_pheno - n_negative_iid, "missing values in phenotype\n")

# Replace missing with -999, and handle invalid IIDs
result[is.na(result$pheno), pheno := -999]
result[result$has_valid_iid == FALSE, pheno := -999]

# Create output ID column (FID_IID format to match C version)
result[, id_output := paste(FID, IID, sep = "_")]

# Write output
cat("Writing output to:", pheno_out, "\n")
fwrite(result[, .(id_output, pheno)], 
       file = pheno_out, 
       sep = " ", 
       col.names = FALSE, 
       quote = FALSE)

cat("Output written to", pheno_out, "\n")
cat("Individuals without phenotype data:", n_missing_pheno - n_negative_iid, "\n")
cat("Individuals with negative IID:", n_negative_iid, "\n")
