# HE Regression Program

This program performs Haseman-Elston (HE) regression analysis on GRM-phenotype pairs within specified genetic relatedness intervals, with bootstrapped confidence intervals.

## Features

- **Interval-based regression**: Define custom GRM intervals (e.g., 0-0.05, 0.05-0.35, 0.35-0.75)
- **Automatic format detection**: Auto-detects GRM and phenotype columns in input files
- **Bootstrap confidence intervals**: Calculates 95% confidence intervals for regression coefficients
- **Comprehensive statistics**: Reports intercept, slope, R², MSE, standard errors, and confidence intervals
- **Flexible input**: Works with output from grm_pheno_join or similar joined datasets

## Usage

```bash
./he_regression -f <input_file> -i <intervals> [options]
```

### Required Arguments
- `-f <file>`: Input file with GRM values and phenotype cross-products
- `-i <intervals>`: Comma-separated intervals (e.g., "0-0.05,0.05-0.35,0.35-0.75")

### Optional Arguments
- `-b <number>`: Bootstrap replicates (default: 1000)
- `-o <file>`: Output file (default: he_regression_results.txt)
- `-s <seed>`: Random seed for reproducible bootstrap results
- `-h`: Show help message

## Example

```bash
# Basic usage with three intervals
./he_regression -f my_results/p21001_i0_cov1234_w0.01_joined.txt -i "0-0.05,0.05-0.35,0.35-0.75"

# With custom bootstrap replicates and output file
./he_regression -f joined_data.txt -i "0-1,0.05-0.35,0.35-0.75" -b 5000 -o my_he_results.txt

# Using a specific random seed for reproducibility
./he_regression -f joined_data.txt -i "0-0.05,0.05-0.35,0.35-0.75" -s 12345
```

## Input File Format

The program expects a tab or space-delimited file with:
- Headers that contain "grm", "relatedness", or "kinship" for the GRM column
- Headers that contain "pheno", "product", or "cross" for the phenotype cross-product column
- If headers don't match, assumes GRM in column 3 and phenotype in column 4

Example input file:
```
iid1    iid2    grm_value    phenotype_crossproduct
1001    1002    0.023       1.45
1001    1003    0.156       2.33
1002    1003    0.089       0.87
...
```

## Output Format

Tab-delimited file with columns:
- `interval`: GRM interval name (e.g., "0.000-0.050")
- `n_pairs`: Number of pairs in the interval
- `intercept`: Regression intercept
- `slope`: Regression slope (heritability estimate)
- `r_squared`: R-squared value
- `mse`: Mean squared error
- `intercept_se`: Standard error of intercept
- `slope_se`: Standard error of slope
- `intercept_ci_lower/upper`: 95% CI for intercept (bootstrap)
- `slope_ci_lower/upper`: 95% CI for slope (bootstrap)

## HE Regression Theory

The Haseman-Elston regression models the relationship between genetic relatedness and phenotype similarity:

```
E[Yi * Yj] = β₀ + β₁ * π_ij
```

Where:
- Yi, Yj are phenotype values for individuals i and j
- π_ij is the genetic relatedness (GRM value)
- β₁ estimates the heritability (h²) directly when using cross-products of residualized phenotypes

## Compilation

```bash
make
```

Requirements: GCC compiler with C99 support and math library (-lm)

## Integration with Pipeline

This program works seamlessly with the existing pipeline:

1. **Residualize phenotypes**: Use `regress_y` to remove covariate effects
2. **Calculate cross-products**: Use `pheno_pairs` to compute phenotype cross-products
3. **Join with GRM**: Use `grm_pheno_join` to combine GRM and cross-products
4. **HE regression**: Use this program to perform interval-based regression analysis

Example workflow:
```bash
# Step 1-3: Run existing pipeline
./residualize_phenotypes.sh phenotypes.txt covariates.txt
./binned_analysis.sh residualized_phenotypes.txt

# Step 4: HE regression on results
./he_regression -f my_results/p21001_i0_cov1234_w0.01_joined.txt -i "0-0.05,0.05-0.35,0.35-0.75"

# Step 5: Visualize results
Rscript plot_he_regression.R he_regression_results.txt bmi_he_analysis
```

## Visualization

The `plot_he_regression.R` script creates comprehensive visualizations:

```bash
Rscript plot_he_regression.R <results_file> [output_prefix]
```

This generates:
- **Heritability plot**: Shows h² estimates with 95% confidence intervals by GRM interval
- **R-squared plot**: Model fit quality across intervals  
- **Sample sizes plot**: Number of pairs in each interval
- **Summary table**: Tab-delimited summary of all results

The heritability estimates are calculated as h² = slope/2 (appropriate for phenotype cross-products of residualized traits).
