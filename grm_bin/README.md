

# GRM binning + phenotype covariance

This directory contains a small, fast pipeline to:

1. Align a phenotype file to a GRM (`.grm.id`)
2. Compute phenotype cross-products for all pairs
3. Bin those pairs by genetic relatedness
4. Output average GRM and phenotype covariance per bin
5. Optionally estimate bin-level noise with a delete-block jackknife
6. Run small phenotype simulations on top of a real genotype panel

The implementation is intentionally simple, reproducible, and HPC-friendly.

---

## Files

- `grm_bin_cov.cpp`
  - C++ program that reads an **aligned phenotype file** and a PLINK GRM
  - Streams through `.grm.bin` and computes binned averages

- `prepare_pheno.R`
  - Preprocessing script
  - Matches phenotype IDs to `.grm.id`
  - Reorders phenotype into GRM order
  - Writes aligned phenotype file

- `grm_bin_cov_jackknife.cpp`
  - Same binning pass, but also computes a delete-block jackknife SE per bin

- `simulate_pheno_from_grm.R`
  - Reads a realized GRM and simulates `y ~ N(0, h2 * GRM + (1-h2) * I)`
  - Writes phenotypes in GRM order plus per-bin truth under the simulation model

- `summarize_simulations.R`
  - Aggregates replicate outputs
  - Compares empirical replicate-to-replicate variability to jackknife SEs

- `run_simulations.sh`
  - PLINK-first end-to-end driver for building a small simulation panel from a `.bed/.bim/.fam` dataset
  - Uses GCTA only for the final `.grm.bin/.grm.id` build step because that is the format consumed here

- `Makefile`
  - Compiles the C++ program

---

## Input formats

### GRM (`.grm.bin` / `.grm.id` format)

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
Rscript prepare_pheno.R \
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

### Step 4: Optional jackknife SEs

```bash
./grm_bin_cov_jackknife aligned_pheno.txt prefix output_jk.tsv 50
```

The last argument is the number of contiguous individual blocks for the delete-block jackknife.

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

For the jackknife binary, the output has one additional column:

```
se
```

---

## Binning scheme

- One coarse bin: `[-0.3, -0.02)`
- Fine bins: width `0.001` between `-0.02` and `0.02`
- Coarser bins: width `0.005` from `0.02` to `1.05`

---

## Assumptions

- GRM is lower-triangular binary relatedness data in the `.grm.bin`/`.grm.id` layout expected by this code
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
Rscript prepare_pheno.R pheno.txt grm.grm.id aligned.txt
make
./grm_bin_cov aligned.txt grm output.tsv
```

## Simulation workflow

The simplest way to test whether the binned means behave correctly is:

1. Start from a real genotype panel such as 1000 Genomes (`.bed/.bim/.fam`)
2. Use PLINK to subset individuals and SNPs, and LD-prune to keep the run fast
3. Build a GRM once
4. Simulate phenotypes conditional on that realized GRM
5. Run `grm_bin_cov` and `grm_bin_cov_jackknife`
6. Compare the observed bin means to the known truth

Under the simulation model

```text
y ~ N(0, h2 * G + (1-h2) * I)
```

the conditional truth for any off-diagonal pair is

```text
E[y_i y_j | G] = h2 * G_ij
```

So for each bin, the expected mean cross-product is just `h2 * avg_grm`.

### End-to-end driver

```bash
./run_simulations.sh /path/to/1000g_prefix ./sim_1kg 800 20000 20 0.4 123
```

Arguments:

- `input_bfile`: input PLINK prefix
- `out_dir`: output directory
- `n_individuals`: number of sampled individuals
- `n_snps`: number of SNPs after pruning to keep
- `n_reps`: number of phenotype replicates
- `h2`: simulation heritability
- `seed`: top-level seed

Environment variables:

- `PLINK` (default `plink`)
- `GCTA` (default `gcta64`)
- `MAF` (default `0.05`)
- `GENO` (default `0.01`)
- `JK_BLOCKS` (default `50`)

The script does the following:

1. random individual subsample from the `.fam`
2. PLINK QC and LD pruning
3. PLINK extraction of a small simulation panel
4. GCTA `.grm.bin` construction
5. repeated phenotype simulation and binning
6. replicate summary

### Outputs

- `simulation_panel.*`: pruned genotype subset and GRM
- `rep_k/pheno.txt`: one simulated phenotype replicate
- `rep_k/truth.tsv`: bin-level truth for that replicate
- `rep_k/binned.tsv`: `grm_bin_cov` output
- `rep_k/jackknife.tsv`: `grm_bin_cov_jackknife` output
- `simulation_summary.tsv`: replicate-aggregated summary

The summary file includes:

- `expected_avg_pheno_crossprod`: truth implied by `h2 * avg_grm`
- `mean_estimate`: average binned estimate across replicates
- `empirical_sd`: replicate-to-replicate variability
- `mean_jackknife_se`: average jackknife SE reported by the fast C++ code
- `mean_oracle_se_indep`: independence-based oracle SE, useful only as a lower-dependence baseline

`empirical_sd` is the main calibration target. `mean_jackknife_se` should usually be in the same range, while `mean_oracle_se_indep` is expected to be optimistic because pairs within a bin are dependent.

---

## Extensions (easy to add)

- weighted bins
- alternative binning schemes
- output full distribution instead of averages
- parallelization (OpenMP over rows)

---

## Author note

This is designed for fast exploratory analysis of GRM vs phenotype covariance relationships, especially for large datasets (e.g. UKB-scale pairwise analyses).
