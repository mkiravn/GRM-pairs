# DNAnexus Example Commands

This directory contains sample DNAnexus commands for running the GRM-pairs
analysis pipeline on the UK Biobank Research Analysis Platform (RAP).

The pipeline covers the full workflow:

1. Uploading input data and compiled binaries
2. Residualising phenotypes (`regress_y`)
3. Computing pairwise GRM values from a list of IDs (`grm_pairs`)
4. Computing phenotype cross-products (`pheno_pairs`)
5. Joining and binning GRM + phenotype data (`grm_pheno_join`)
6. Haseman-Elston regression (`he_regression`)
7. Downloading results

All steps run as **Swiss Army Knife** (`swiss-army-knife`) jobs on a DNAnexus
Ubuntu 24.04 x86_64 node unless otherwise noted.

---

## Prerequisites

```bash
# Log in to DNAnexus
dx login

# Select your project
dx select "My UKB Project"
```

---

## 1. Upload Compiled Binaries and Input Files

Build each tool locally (or on any x86_64 Ubuntu 24.04 machine) and upload the
resulting binaries to your DNAnexus project.

```bash
# Build all tools
for tool in regress_y grm_pairs pheno_pairs grm_pheno_join he_regression; do
    (cd $tool && make)
done

# Upload binaries to a dedicated folder in your project
dx mkdir -p /tools
dx upload regress_y/regress_y         --path /tools/
dx upload grm_pairs/grm_pairs         --path /tools/
dx upload pheno_pairs/pheno_pairs     --path /tools/
dx upload grm_pheno_join/grm_pheno_join --path /tools/
dx upload he_regression/he_regression --path /tools/

# Upload PLINK binary files (bed/bim/fam)
dx mkdir -p /data/plink
dx upload my_data.bed --path /data/plink/
dx upload my_data.bim --path /data/plink/
dx upload my_data.fam --path /data/plink/

# Upload allele frequency file (produced by plink2 --freq or plink --freq)
dx upload my_data.afreq --path /data/plink/   # plink2
# or
dx upload my_data.frq   --path /data/plink/   # plink1.9

# Upload phenotype and covariate files
dx mkdir -p /data/pheno
dx upload phenotypes.txt  --path /data/pheno/
dx upload covariates.txt  --path /data/pheno/

# Upload pair list (one "IID1  IID2" per line, no header required)
dx upload pairs.txt --path /data/
```

---

## 2. Residualise Phenotypes (`regress_y`)

Regress out covariate effects, remove outliers (>5 SD), and apply
sex-stratified normalisation.

```bash
dx run swiss-army-knife \
    -iin="/tools/regress_y" \
    -iin="/data/pheno/phenotypes.txt" \
    -iin="/data/pheno/covariates.txt" \
    -icmd='
        chmod +x regress_y
        # Regress out sex (col 1) and age (col 2)
        ./regress_y phenotypes.txt covariates.txt "1,2" residualized.txt
    ' \
    --instance-type mem1_ssd1_v2_x4 \
    --destination /results/ \
    --name "regress_y" \
    -y
```

**Output:** `residualized.txt` — tab-delimited file with columns `IID` and
`normalized_residual` (plus additional phenotype columns if present).

---

## 3. Compute Pairwise GRM from a List of IDs (`grm_pairs`)

Compute genetic relatedness for a user-specified list of pairs without building
the full N×N matrix.

```bash
dx run swiss-army-knife \
    -iin="/tools/grm_pairs" \
    -iin="/data/plink/my_data.bed" \
    -iin="/data/plink/my_data.bim" \
    -iin="/data/plink/my_data.fam" \
    -iin="/data/plink/my_data.afreq" \
    -iin="/data/pairs.txt" \
    -icmd='
        chmod +x grm_pairs
        # bfile prefix must match the uploaded .bed/.bim/.fam names
        ./grm_pairs my_data pairs.txt grm_output.txt
    ' \
    --instance-type mem1_ssd1_v2_x4 \
    --destination /results/ \
    --name "grm_pairs" \
    -y
```

**Notes:**
- The program auto-detects whether a `.afreq` (plink2) or `.frq` (plink1.9)
  frequency file is present; upload whichever format you have.
- `pairs.txt` contains whitespace-delimited `IID1 IID2` pairs, one per line.

**Output:** `grm_output.txt` — tab-delimited file with columns
`IID1  IID2  N_SNPs  GRM`.

---

## 4. Compute Phenotype Cross-Products (`pheno_pairs`)

Calculate the product Yi × Yj for every pair and phenotype column.

```bash
dx run swiss-army-knife \
    -iin="/tools/pheno_pairs" \
    -iin="/data/pairs.txt" \
    -iin="/results/residualized.txt" \
    -icmd='
        chmod +x pheno_pairs
        ./pheno_pairs pairs.txt residualized.txt pheno_cross.txt
    ' \
    --instance-type mem1_ssd1_v2_x4 \
    --destination /results/ \
    --name "pheno_pairs" \
    -y
```

**Output:** `pheno_cross.txt` — tab-delimited file with columns
`IID1  IID2  crossproduct_<pheno1>  crossproduct_<pheno2>  ...`

---

## 5. Join and Bin GRM + Phenotype Data (`grm_pheno_join`)

Merge GRM and phenotype cross-product files, then aggregate statistics into
GRM-value bins.

```bash
dx run swiss-army-knife \
    -iin="/tools/grm_pheno_join" \
    -iin="/results/grm_output.txt" \
    -iin="/results/pheno_cross.txt" \
    -icmd='
        chmod +x grm_pheno_join
        # Default binning: bin_width=0.001, GRM range [-0.05, 1.3]
        ./grm_pheno_join grm_output.txt pheno_cross.txt joined.txt

        # Custom binning (bin_width=0.01, range -0.05 to 1.0):
        # ./grm_pheno_join grm_output.txt pheno_cross.txt joined.txt 0.01 -0.05 1.0
    ' \
    --instance-type mem1_ssd1_v2_x4 \
    --destination /results/ \
    --name "grm_pheno_join" \
    -y
```

**Arguments (all optional after the three required positional arguments):**

| Argument | Default | Description |
|---|---|---|
| `bin_width` | `0.001` | Width of each GRM bin |
| `min_grm` | `-0.05` | Lower bound of GRM range |
| `max_grm` | `1.3` | Upper bound of GRM range |

**Output:** `joined.txt` — per-bin statistics (count, sum, sum-of-squares per
GRM interval), ready for HE regression.

---

## 6. Haseman-Elston Regression (`he_regression`)

Fit E[Yi × Yj] = β₀ + β₁ × π_ij within specified GRM intervals and estimate
heritability (h² = slope / 2).

```bash
dx run swiss-army-knife \
    -iin="/tools/he_regression" \
    -iin="/results/joined.txt" \
    -icmd='
        chmod +x he_regression
        # Intervals: unrelated, distant relatives, close relatives
        ./he_regression \
            -f joined.txt \
            -i "0-0.05,0.05-0.35,0.35-0.75" \
            -b 1000 \
            -s 42 \
            -o he_results.txt
    ' \
    --instance-type mem1_ssd1_v2_x4 \
    --destination /results/ \
    --name "he_regression" \
    -y
```

**Options:**

| Flag | Default | Description |
|---|---|---|
| `-f <file>` | — | Input file (output of `grm_pheno_join`) |
| `-i <intervals>` | — | Comma-separated GRM intervals, e.g. `"0-0.05,0.05-0.35"` |
| `-b <n>` | `1000` | Number of bootstrap replicates |
| `-s <seed>` | random | Random seed for reproducible bootstraps |
| `-o <file>` | `he_regression_results.txt` | Output file |

**Output:** `he_results.txt` — tab-delimited file with one row per interval
containing: `interval`, `n_pairs`, `intercept`, `slope`, `r_squared`, `mse`,
standard errors, and 95 % bootstrap confidence intervals.

---

## 7. Download Results

```bash
# Download all result files to a local directory
mkdir -p ./grm_results
dx download /results/grm_output.txt   -o ./grm_results/
dx download /results/pheno_cross.txt  -o ./grm_results/
dx download /results/joined.txt       -o ./grm_results/
dx download /results/he_results.txt   -o ./grm_results/
```

---

## End-to-End Shell Script

The following script submits all pipeline steps as a sequential batch of
Swiss Army Knife jobs. Each step waits for the previous one to finish using
`--depends-on`.

```bash
#!/usr/bin/env bash
# run_pipeline.sh  --  submit the full GRM-pairs pipeline on DNAnexus
# Usage: bash run_pipeline.sh

set -euo pipefail

PROJECT="My UKB Project"
DEST="/results"

# ── 1. Residualise phenotypes ────────────────────────────────────────────────
JOB1=$(dx run swiss-army-knife \
    -iin="/tools/regress_y" \
    -iin="/data/pheno/phenotypes.txt" \
    -iin="/data/pheno/covariates.txt" \
    -icmd='
        chmod +x regress_y
        ./regress_y phenotypes.txt covariates.txt "1,2" residualized.txt
    ' \
    --instance-type mem1_ssd1_v2_x4 \
    --destination "${DEST}" \
    --name "1_regress_y" \
    -y --brief)

echo "regress_y job: ${JOB1}"

# ── 2. Compute pairwise GRM ──────────────────────────────────────────────────
JOB2=$(dx run swiss-army-knife \
    -iin="/tools/grm_pairs" \
    -iin="/data/plink/my_data.bed" \
    -iin="/data/plink/my_data.bim" \
    -iin="/data/plink/my_data.fam" \
    -iin="/data/plink/my_data.afreq" \
    -iin="/data/pairs.txt" \
    -icmd='
        chmod +x grm_pairs
        ./grm_pairs my_data pairs.txt grm_output.txt
    ' \
    --instance-type mem1_ssd1_v2_x4 \
    --destination "${DEST}" \
    --name "2_grm_pairs" \
    -y --brief)

echo "grm_pairs job: ${JOB2}"

# ── 3. Compute phenotype cross-products ─────────────────────────────────────
JOB3=$(dx run swiss-army-knife \
    -iin="/tools/pheno_pairs" \
    -iin="/data/pairs.txt" \
    -iin="${DEST}/residualized.txt" \
    -icmd='
        chmod +x pheno_pairs
        ./pheno_pairs pairs.txt residualized.txt pheno_cross.txt
    ' \
    --depends-on "${JOB1}" \
    --instance-type mem1_ssd1_v2_x4 \
    --destination "${DEST}" \
    --name "3_pheno_pairs" \
    -y --brief)

echo "pheno_pairs job: ${JOB3}"

# ── 4. Join and bin ──────────────────────────────────────────────────────────
JOB4=$(dx run swiss-army-knife \
    -iin="/tools/grm_pheno_join" \
    -iin="${DEST}/grm_output.txt" \
    -iin="${DEST}/pheno_cross.txt" \
    -icmd='
        chmod +x grm_pheno_join
        ./grm_pheno_join grm_output.txt pheno_cross.txt joined.txt 0.001
    ' \
    --depends-on "${JOB2}" \
    --depends-on "${JOB3}" \
    --instance-type mem1_ssd1_v2_x4 \
    --destination "${DEST}" \
    --name "4_grm_pheno_join" \
    -y --brief)

echo "grm_pheno_join job: ${JOB4}"

# ── 5. HE regression ────────────────────────────────────────────────────────
JOB5=$(dx run swiss-army-knife \
    -iin="/tools/he_regression" \
    -iin="${DEST}/joined.txt" \
    -icmd='
        chmod +x he_regression
        ./he_regression \
            -f joined.txt \
            -i "0-0.05,0.05-0.35,0.35-0.75" \
            -b 1000 \
            -s 42 \
            -o he_results.txt
    ' \
    --depends-on "${JOB4}" \
    --instance-type mem1_ssd1_v2_x4 \
    --destination "${DEST}" \
    --name "5_he_regression" \
    -y --brief)

echo "he_regression job: ${JOB5}"

echo ""
echo "All jobs submitted. Monitor progress with:"
echo "  dx watch ${JOB5}"
echo ""
echo "Download results when complete:"
echo "  dx download ${DEST}/he_results.txt"
```

---

## Monitoring Jobs

```bash
# Watch a running job's logs in real time
dx watch <job-id>

# List recent jobs and their status
dx find jobs --project "My UKB Project" --state "*" --brief

# Check the state of a specific job
dx describe <job-id> | grep state
```

---

## Tips

- **Instance types**: `mem1_ssd1_v2_x4` (4 cores, 7.5 GB RAM) is sufficient for
  most steps. For very large pair lists (>100 M pairs), consider
  `mem2_ssd1_v2_x8` or larger.
- **Frequency files**: `grm_pairs` auto-detects `.afreq` (plink2) vs `.frq`
  (plink1.9). Upload whichever format your pipeline produces; both work.
- **Allele-flip check**: `grm_pairs` compares allele coding in the frequency
  file against the `.bim` file and flips frequencies automatically when needed —
  no manual intervention required.
- **Custom binning**: pass `bin_width`, `min_grm`, and `max_grm` to
  `grm_pheno_join` to adjust resolution or GRM range for your dataset.
- **Reproducibility**: use `-s <seed>` in `he_regression` so that bootstrap
  confidence intervals are reproducible across runs.
