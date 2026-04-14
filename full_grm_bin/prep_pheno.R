library(readr)
library(dplyr)

standardize_within_sex <- function(x, sex) {
  ave(
    x, sex,
    FUN = function(z) {
      if (all(is.na(z))) return(rep(NA_real_, length(z)))
      s <- sd(z, na.rm = TRUE)
      m <- mean(z, na.rm = TRUE)
      if (is.na(s) || s == 0) return(rep(NA_real_, length(z)))
      (z - m) / s
    }
  )
}

write_grm_pheno <- function(df, y_col, file) {
  out <- df %>%
    transmute(
      FID = participant.eid,
      IID = participant.eid,
      Y = .data[[y_col]]
    )

  write.table(
    out,
    file = file,
    quote = FALSE,
    sep = ' ',
    row.names = FALSE,
    col.names = TRUE,
    na = 'NA'
  )
}

pheno <- read.delim('pheno_raw.tsv', sep = '\t', header = TRUE, check.names = FALSE)
covars <- read.delim('all_covars.tsv', sep = '\t', header = TRUE, check.names = FALSE)
pcs <- read.delim('ukb_subset_pca.eigenvec', header = TRUE, check.names = FALSE)

names(pheno) <- c('participant.eid', 'phenotype')
pheno$participant.eid <- as.character(pheno$participant.eid)
covars$participant.eid <- as.character(covars$participant.eid)

# Handle PLINK eigenvec format robustly
if ('IID' %in% names(pcs)) {
  pcs$participant.eid <- as.character(pcs$IID)
} else if ('participant.eid' %in% names(pcs)) {
  pcs$participant.eid <- as.character(pcs$participant.eid)
} else {
  stop('Could not identify participant ID column in ukb_subset_pca.eigenvec')
}

pc_names <- paste0('PC', 1:5)
if (!all(pc_names %in% names(pcs))) {
  stop(paste('Missing PCs in eigenvec file:', paste(setdiff(pc_names, names(pcs)), collapse = ', ')))
}

dat <- pheno %>%
  left_join(covars, by = 'participant.eid') %>%
  left_join(pcs %>% select(participant.eid, all_of(pc_names)), by = 'participant.eid')

dat <- dat %>%
  mutate(
    participant.eid = as.character(participant.eid),
    phenotype = suppressWarnings(as.numeric(phenotype))
  )

# 5 SD phenotype filter, Kemper-style
y_mean <- mean(dat$phenotype, na.rm = TRUE)
y_sd <- sd(dat$phenotype, na.rm = TRUE)

dat <- dat %>%
  mutate(
    phenotype_clean = ifelse(
      abs(phenotype - y_mean) > 5 * y_sd,
      NA_real_,
      phenotype
    ),
    sex = factor(sex),
    age = factor(age),
    yob = factor(yob),
    batch = factor(batch),
    gc = factor(gc)
  )

f1 <- as.formula('phenotype_clean ~ sex + age + yob + batch')
f2 <- as.formula('phenotype_clean ~ sex + age + yob + batch + gc')
f3 <- as.formula(
  paste('phenotype_clean ~ sex + age + yob + batch +', paste(pc_names, collapse = ' + '))
)
f4 <- as.formula(
  paste('phenotype_clean ~ sex + age + yob + batch + gc +', paste(pc_names, collapse = ' + '))
)

fit1 <- lm(f1, data = dat, na.action = na.exclude)
fit2 <- lm(f2, data = dat, na.action = na.exclude)
fit3 <- lm(f3, data = dat, na.action = na.exclude)
fit4 <- lm(f4, data = dat, na.action = na.exclude)

cat('==============================\n')
cat('Residualization summaries\n')
cat('==============================\n\n')

cat('Model m1: phenotype_clean ~ sex + age + yob + batch\n')
print(summary(fit1))
cat('\n')

cat('Model m2: phenotype_clean ~ sex + age + yob + batch + gc\n')
print(summary(fit2))
cat('\n')

cat('Model m3: phenotype_clean ~ sex + age + yob + batch + PC1-5\n')
print(summary(fit3))
cat('\n')

cat('Model m4: phenotype_clean ~ sex + age + yob + batch + gc + PC1-5\n')
print(summary(fit4))
cat('\n')

dat$resid_m1 <- standardize_within_sex(residuals(fit1), dat$sex)
dat$resid_m2 <- standardize_within_sex(residuals(fit2), dat$sex)
dat$resid_m3 <- standardize_within_sex(residuals(fit3), dat$sex)
dat$resid_m4 <- standardize_within_sex(residuals(fit4), dat$sex)

write_grm_pheno(dat, 'resid_m1', 'pheno_input_m1.txt')
write_grm_pheno(dat, 'resid_m2', 'pheno_input_m2.txt')
write_grm_pheno(dat, 'resid_m3', 'pheno_input_m3.txt')
write_grm_pheno(dat, 'resid_m4', 'pheno_input_m4.txt')