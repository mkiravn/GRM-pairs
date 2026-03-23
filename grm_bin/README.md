

# GRM binning + phenotype covariance

This directory contains a small, fast pipeline to:

1. Align a phenotype file to a GRM (.grm.id)
2. Standardize the phenotype
3. Compute phenotype cross-products for all pairs
4. Bin those pairs by genetic relatedness (GRM)
5. Output average GRM and phenotype covariance per bin

The implementation is intentionally simple, reproducible, and HPC-friendly.

---

## Files

- `grm_bin_cov.cpp`
  - C++ program that reads an **aligned phenotype file** and a PLINK GRM
  - Streams through `.grm.bin` and computes binned averages

- `prepare_pheno_for_grm.R`
  - Preprocessing script
  - Matches phenotype IDs to `.grm.id`
  - Reorders phenotype into GRM order
  - Writes aligned phenotype file

- `Makefile`
  - Compiles the C++ program

---

## Input formats

### GRM (PLINK format)

Expected files:

- `prefix.grm.bin`
- `prefix.grm.id`

`.grm.id` must contain:

```
FID IID
```

---

### Phenotype file

Two supported formats:

#### Option 1: IID only
```
IID phenotype
```

#### Option 2: FID IID
```
FID IID phenotype
```

Missing values can be encoded as:

- `NA`
- `NaN`
- `-999`

---

## Workflow

### Step 1: Align phenotype to GRM

```bash
Rscript prepare_pheno_for_grm.R \
    phenotype.txt \
    prefix.grm.id \
    aligned_pheno.txt
```

This:

- matches IDs
- reorders to GRM order
- inserts `NA` where missing

---

### Step 2: Compile

```bash
make
```

---

### Step 3: Run binning

```bash
./grm_bin_cov aligned_pheno.txt prefix output.tsv
```

---

## Output

Tab-separated file:

```
bin    lower    upper    avg_grm    avg_pheno_crossprod    n_pairs
```

Where:

- `avg_grm` = mean GRM value in bin
- `avg_pheno_crossprod` = mean of y_i * y_j
- `n_pairs` = number of pairs in bin

---

## Binning scheme

- One coarse bin: `[-0.3, -0.02)`
- Fine bins: width `0.001` between `-0.02` and `0.02`
- Coarser bins: width `0.005` from `0.02` to `1.05`

---

## Assumptions

- GRM is lower-triangular PLINK `.grm.bin`
- Phenotype file is aligned before C++ step
- No ID checking inside C++ (done in R step)

---

## Performance notes

- GRM is streamed sequentially (fast, no random access)
- Memory usage is minimal
- Runtime is dominated by number of pairs (~n²/2)

---

## Common pitfalls

- Misaligned phenotype and GRM → incorrect results
- Missing values not encoded properly
- `.grm.bin` not matching `.grm.id`

---

## Minimal example

```bash
Rscript prepare_pheno_for_grm.R pheno.txt grm.grm.id aligned.txt
make
./grm_bin_cov aligned.txt grm output.tsv
```

---

## Extensions (easy to add)

- weighted bins
- alternative binning schemes
- output full distribution instead of averages
- parallelization (OpenMP over rows)

---

## Author note

This is designed for fast exploratory analysis of GRM vs phenotype covariance relationships, especially for large datasets (e.g. UKB-scale pairwise analyses).