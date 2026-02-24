# Regress Y

A C program for phenotype regression, outlier removal, and sex-stratified normalization.

## Compilation

```bash
make
```

or

```bash
gcc -Wall -Wextra -O2 -std=c99 -lm -o regress_y main.c
```

## Usage

```bash
./regress_y <pheno_file> <covar_file> <covar_indices> <out_file>
```

### Input Files

1. **Phenotype file**: Tab or space-delimited file with individual IDs and phenotype values
   ```
   IID     height
   sample1 175.2
   sample2 182.1
   sample3 168.5
   ```

2. **Covariate file**: Tab or space-delimited file with individual IDs and covariate values
   ```
   IID     sex  age  bmi
   sample1 0    25   23.5
   sample2 1    30   21.8
   sample3 0    28   24.1
   ```

### Arguments

- `<pheno_file>`: Phenotype file (IID pheno)
- `<covar_file>`: Covariate file (IID cov1 cov2 cov3 ...)
- `<covar_indices>`: Comma-separated covariate column indices (1-based, e.g., "1,3")
- `<out_file>`: Output file with normalized residuals

### Process

1. **Regression**: Fits `pheno = intercept + selected_covariates + residual`
2. **Outlier removal**: Removes observations >5 standard deviations from residual mean
3. **Sex-stratified normalization**: Normalizes residuals to zero mean and unit variance within each sex
   - Assumes first covariate is sex (0=male, 1=female)

### Output

Tab-delimited file with normalized residuals:
```
IID     normalized_residual
sample1 0.234
sample2 -1.156
sample3 0.922
```

## Example

```bash
# Regress out sex and age effects
./regress_y phenotypes.txt covariates.txt "1,2" normalized_residuals.txt

# Regress out sex, age, and BMI effects  
./regress_y phenotypes.txt covariates.txt "1,2,3" normalized_residuals.txt
```

## Features

- **Flexible covariate selection**: Select any combination of covariates by index
- **Robust regression**: Uses normal equations for multiple linear regression
- **Outlier detection**: Removes extreme values (>5 SD) before normalization
- **Sex-stratified normalization**: Separate normalization for males and females
- **Comprehensive logging**: Reports regression coefficients, outliers removed, and normalization statistics
- **Memory efficient**: Handles large datasets efficiently

## Mathematical Details

### Regression Model
```
y = Xβ + ε
```
where:
- `y` is the phenotype vector
- `X` is the design matrix (intercept + selected covariates)
- `β` is the coefficient vector (solved using normal equations: β = (X'X)^(-1)X'y)
- `ε` is the residual vector

### Normalization
For each sex group s:
```
normalized_residual = (residual - mean_s) / sd_s
```

This ensures that the final residuals have mean=0 and variance=1 within each sex group.
