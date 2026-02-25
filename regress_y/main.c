/*
 * regress_y.c -- Phenotype regression, outlier removal, and normalization
 * ROBUST VERSION: Multi-covariate support with complete case analysis
 *
 * Usage: regress_y <pheno_file> <covar_file> <covar_indices> <out_file>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/* Structure to store individual data */
typedef struct {
    char* iid;
    double pheno;
    double* covars;
    uint32_t n_covars;
    double residual;
    int sex; /* 0=male, 1=female, -1=unknown */
    int valid; /* 1 if valid complete case */
} individual_t;

/* Parse comma-separated indices string */
static uint32_t* parse_indices(const char* indices_str, uint32_t* n_indices_out)
{
    if (!indices_str || strlen(indices_str) == 0) {
        *n_indices_out = 0;
        return NULL;
    }
    
    /* Count commas to estimate number of indices */
    uint32_t n_commas = 0;
    for (const char* p = indices_str; *p; p++) {
        if (*p == ',') n_commas++;
    }
    uint32_t max_indices = n_commas + 1;
    
    uint32_t* indices = malloc(max_indices * sizeof(uint32_t));
    if (!indices) return NULL;
    
    uint32_t n_indices = 0;
    char* str_copy = malloc(strlen(indices_str) + 1);
    strcpy(str_copy, indices_str);
    
    char* token = strtok(str_copy, ",");
    while (token && n_indices < max_indices) {
        int idx = atoi(token);
        if (idx > 0) { /* Convert to 0-based indexing */
            indices[n_indices++] = idx - 1;
        }
        token = strtok(NULL, ",");
    }
    
    free(str_copy);
    *n_indices_out = n_indices;
    
    if (n_indices == 0) {
        free(indices);
        return NULL;
    }
    
    return indices;
}

/* Read covariate file */
static individual_t* read_covariate_file(const char* covar_fname, 
                                         uint32_t* count_out,
                                         char*** covar_names_out,
                                         uint32_t* n_covars_out)
{
    FILE* f = fopen(covar_fname, "r");
    if (!f) { perror(covar_fname); return NULL; }

    char line[2048];
    uint32_t n_covars = 0;
    char** covar_names = NULL;
    
    /* Read first line to get column headers */
    if (!fgets(line, sizeof(line), f)) {
        fprintf(stderr, "Error: Empty file %s\n", covar_fname);
        fclose(f);
        return NULL;
    }
    
    /* Count columns */
    char line_copy[2048];
    strcpy(line_copy, line);
    char* token = strtok(line_copy, " \t\n\r");
    while (token) {
        n_covars++;
        token = strtok(NULL, " \t\n\r");
    }
    
    if (n_covars < 2) {
        fprintf(stderr, "Error: Need at least 2 columns (IID + covariates)\n");
        fclose(f);
        return NULL;
    }
    
    /* Store column names */
    covar_names = malloc(n_covars * sizeof(char*));
    strcpy(line_copy, line);
    token = strtok(line_copy, " \t\n\r");
    for (uint32_t i = 0; i < n_covars && token; i++) {
        covar_names[i] = malloc(strlen(token) + 1);
        strcpy(covar_names[i], token);
        token = strtok(NULL, " \t\n\r");
    }
    
    /* Count data rows */
    uint32_t row_count = 0;
    while (fgets(line, sizeof(line), f)) {
        if (strlen(line) > 1) row_count++;
    }
    
    if (row_count == 0) {
        fprintf(stderr, "Error: No data rows in %s\n", covar_fname);
        fclose(f);
        return NULL;
    }
    
    /* Allocate individuals array */
    individual_t* individuals = calloc(row_count, sizeof(individual_t));
    if (!individuals) {
        perror("malloc");
        fclose(f);
        return NULL;
    }
    
    /* Reset file to data rows */
    rewind(f);
    fgets(line, sizeof(line), f); /* Skip header */
    
    /* Read data */
    uint32_t count = 0;
    while (fgets(line, sizeof(line), f) && count < row_count) {
        if (strlen(line) <= 1) continue;
        
        /* Allocate covariates array */
        individuals[count].covars = malloc((n_covars - 1) * sizeof(double));
        individuals[count].n_covars = n_covars - 1;
        individuals[count].valid = 1; /* Assume valid initially */
        
        /* Parse line */
        char* token = strtok(line, " \t\n\r");
        if (!token) continue;
        
        /* First column is IID */
        individuals[count].iid = malloc(strlen(token) + 1);
        strcpy(individuals[count].iid, token);
        
        /* Parse covariates */
        for (uint32_t j = 0; j < n_covars - 1; j++) {
            token = strtok(NULL, " \t\n\r");
            if (!token) {
                individuals[count].valid = 0;
                individuals[count].covars[j] = NAN;
                break;
            }
            
            char* endptr;
            double val = strtod(token, &endptr);
            if (*endptr != '\0' || !isfinite(val)) {
                individuals[count].valid = 0;
                individuals[count].covars[j] = NAN;
            } else {
                individuals[count].covars[j] = val;
            }
        }
        
        count++;
    }
    
    fclose(f);
    
    *count_out = count;
    *covar_names_out = covar_names;
    *n_covars_out = n_covars;
    
    return individuals;
}

/* Read phenotype file and match with individuals */
static uint32_t read_phenotype_file(const char* pheno_fname, individual_t* individuals, 
                                   uint32_t count, char** pheno_name_out)
{
    FILE* f = fopen(pheno_fname, "r");
    if (!f) { 
        perror(pheno_fname); 
        return 0; 
    }

    char line[1024];
    char* pheno_name = NULL;
    
    /* Read header line */
    if (fgets(line, sizeof(line), f)) {
        char* tab_pos = strchr(line, '\t');
        if (tab_pos) {
            pheno_name = malloc(strlen(tab_pos + 1) + 1);
            strcpy(pheno_name, tab_pos + 1);
            /* Remove trailing whitespace */
            char* end = pheno_name + strlen(pheno_name) - 1;
            while (end > pheno_name && (*end == '\n' || *end == '\r' || *end == ' ' || *end == '\t')) {
                *end = '\0';
                end--;
            }
        }
    }
    
    /* Initialize all phenotypes as missing */
    for (uint32_t i = 0; i < count; i++) {
        individuals[i].pheno = NAN;
    }
    
    uint32_t matched = 0;
    
    /* Read phenotype data */
    while (fgets(line, sizeof(line), f)) {
        if (strlen(line) <= 1) continue;
        
        char* iid_token = strtok(line, "\t");
        char* pheno_token = strtok(NULL, "\t\n\r");
        
        if (!iid_token || !pheno_token) continue;
        
        /* Find matching individual */
        for (uint32_t i = 0; i < count; i++) {
            if (strcmp(individuals[i].iid, iid_token) == 0) {
                char* endptr;
                double val = strtod(pheno_token, &endptr);
                if (*endptr == '\0' && isfinite(val)) {
                    individuals[i].pheno = val;
                    matched++;
                } else {
                    individuals[i].valid = 0; /* Invalid phenotype = exclude from analysis */
                }
                break;
            }
        }
    }
    
    fclose(f);
    *pheno_name_out = pheno_name;
    return matched;
}

/* Apply complete case analysis - exclude any individual with missing data */
static uint32_t apply_complete_case_analysis(individual_t* individuals, uint32_t count,
                                           const uint32_t* covar_indices, uint32_t n_covar_indices)
{
    uint32_t complete_cases = 0;
    
    for (uint32_t i = 0; i < count; i++) {
        /* Start with existing validity */
        int valid = individuals[i].valid;
        
        /* Check phenotype */
        if (!isfinite(individuals[i].pheno)) {
            valid = 0;
        }
        
        /* Check all required covariates */
        for (uint32_t c = 0; c < n_covar_indices && valid; c++) {
            if (covar_indices[c] >= individuals[i].n_covars || 
                !isfinite(individuals[i].covars[covar_indices[c]])) {
                valid = 0;
            }
        }
        
        individuals[i].valid = valid;
        if (valid) complete_cases++;
    }
    
    return complete_cases;
}

/* Robust linear regression with complete case analysis */
static int perform_regression(individual_t* individuals, uint32_t count,
                             const uint32_t* covar_indices, uint32_t n_covar_indices,
                             double* beta_out, uint32_t* n_complete_out)
{
    uint32_t n_complete = 0;
    
    /* Count complete cases */
    for (uint32_t i = 0; i < count; i++) {
        if (individuals[i].valid) n_complete++;
    }
    
    if (n_complete == 0) {
        fprintf(stderr, "Error: No complete cases for regression\n");
        return 0;
    }
    
    if (n_complete <= n_covar_indices + 1) {
        fprintf(stderr, "Error: Need more complete cases than parameters (%u cases, %u params)\n", 
                n_complete, n_covar_indices + 1);
        return 0;
    }
    
    uint32_t p = n_covar_indices + 1; /* +1 for intercept */
    
    /* Allocate design matrix and response vector for complete cases only */
    double* X = calloc(n_complete * p, sizeof(double));
    double* y = malloc(n_complete * sizeof(double));
    
    if (!X || !y) {
        perror("malloc");
        free(X);
        free(y);
        return 0;
    }
    
    /* Fill design matrix and response vector with complete cases only */
    uint32_t row = 0;
    for (uint32_t i = 0; i < count; i++) {
        if (!individuals[i].valid) continue;
        
        X[row * p + 0] = 1.0; /* intercept */
        for (uint32_t j = 0; j < n_covar_indices; j++) {
            X[row * p + j + 1] = individuals[i].covars[covar_indices[j]];
        }
        y[row] = individuals[i].pheno;
        row++;
    }
    
    /* Compute X'X and X'y using only complete cases */
    double* XtX = calloc(p * p, sizeof(double));
    double* Xty = calloc(p, sizeof(double));
    
    if (!XtX || !Xty) {
        perror("malloc");
        free(X);
        free(y);
        free(XtX);
        free(Xty);
        return 0;
    }
    
    /* X'X = sum over complete cases of X_i' * X_i */
    for (uint32_t i = 0; i < n_complete; i++) {
        for (uint32_t j = 0; j < p; j++) {
            for (uint32_t k = 0; k < p; k++) {
                XtX[j * p + k] += X[i * p + j] * X[i * p + k];
            }
        }
    }
    
    /* X'y = sum over complete cases of X_i' * y_i */
    for (uint32_t i = 0; i < n_complete; i++) {
        for (uint32_t j = 0; j < p; j++) {
            Xty[j] += X[i * p + j] * y[i];
        }
    }
    
    /* Solve normal equations: (X'X) * beta = X'y using Gaussian elimination */
    /* Create augmented matrix [X'X | X'y] */
    double* aug_matrix = malloc(p * (p + 1) * sizeof(double));
    if (!aug_matrix) {
        perror("malloc");
        free(X);
        free(y);
        free(XtX);
        free(Xty);
        return 0;
    }
    
    for (uint32_t i = 0; i < p; i++) {
        for (uint32_t j = 0; j < p; j++) {
            aug_matrix[i * (p + 1) + j] = XtX[i * p + j];
        }
        aug_matrix[i * (p + 1) + p] = Xty[i];
    }
    
    /* Forward elimination with partial pivoting */
    for (uint32_t i = 0; i < p; i++) {
        /* Find pivot */
        uint32_t max_row = i;
        double max_val = fabs(aug_matrix[i * (p + 1) + i]);
        
        for (uint32_t k = i + 1; k < p; k++) {
            double val = fabs(aug_matrix[k * (p + 1) + i]);
            if (val > max_val) {
                max_val = val;
                max_row = k;
            }
        }
        
        /* Check for near-singularity */
        if (max_val < 1e-12) {
            fprintf(stderr, "Error: Matrix is nearly singular (pivot=%.2e)\n", max_val);
            free(X);
            free(y);
            free(XtX);
            free(Xty);
            free(aug_matrix);
            return 0;
        }
        
        /* Swap rows if needed */
        if (max_row != i) {
            for (uint32_t j = 0; j <= p; j++) {
                double temp = aug_matrix[i * (p + 1) + j];
                aug_matrix[i * (p + 1) + j] = aug_matrix[max_row * (p + 1) + j];
                aug_matrix[max_row * (p + 1) + j] = temp;
            }
        }
        
        /* Eliminate below pivot */
        for (uint32_t k = i + 1; k < p; k++) {
            double factor = aug_matrix[k * (p + 1) + i] / aug_matrix[i * (p + 1) + i];
            for (uint32_t j = i; j <= p; j++) {
                aug_matrix[k * (p + 1) + j] -= factor * aug_matrix[i * (p + 1) + j];
            }
        }
    }
    
    /* Back substitution */
    for (int i = p - 1; i >= 0; i--) {
        beta_out[i] = aug_matrix[i * (p + 1) + p];
        for (uint32_t j = i + 1; j < p; j++) {
            beta_out[i] -= aug_matrix[i * (p + 1) + j] * beta_out[j];
        }
        beta_out[i] /= aug_matrix[i * (p + 1) + i];
    }
    
    /* Verify coefficients are finite */
    for (uint32_t i = 0; i < p; i++) {
        if (!isfinite(beta_out[i])) {
            fprintf(stderr, "Error: Non-finite regression coefficient %u: %.6f\n", i, beta_out[i]);
            free(X);
            free(y);
            free(XtX);
            free(Xty);
            free(aug_matrix);
            return 0;
        }
    }
    
    free(X);
    free(y);
    free(XtX);
    free(Xty);
    free(aug_matrix);
    
    *n_complete_out = n_complete;
    return 1;
}

/* Calculate residuals for valid individuals only */
static void calculate_residuals(individual_t* individuals, uint32_t count,
                               const uint32_t* covar_indices, uint32_t n_covar_indices,
                               const double* beta)
{
    for (uint32_t i = 0; i < count; i++) {
        if (!individuals[i].valid) {
            individuals[i].residual = NAN;
            continue;
        }
        
        /* Calculate predicted value */
        double predicted = beta[0]; /* intercept */
        for (uint32_t j = 0; j < n_covar_indices; j++) {
            predicted += beta[j + 1] * individuals[i].covars[covar_indices[j]];
        }
        
        /* Residual = observed - predicted */
        individuals[i].residual = individuals[i].pheno - predicted;
        
        if (!isfinite(individuals[i].residual)) {
            fprintf(stderr, "Warning: Non-finite residual for individual %s\n", individuals[i].iid);
            individuals[i].valid = 0;
            individuals[i].residual = NAN;
        }
    }
}

/* Remove outliers (>5 SD from mean) from valid individuals */
static uint32_t remove_outliers(individual_t* individuals, uint32_t count)
{
    /* Calculate mean and SD of residuals for valid individuals */
    double sum = 0.0, sum_sq = 0.0;
    uint32_t n_valid = 0;
    
    for (uint32_t i = 0; i < count; i++) {
        if (individuals[i].valid && isfinite(individuals[i].residual)) {
            sum += individuals[i].residual;
            sum_sq += individuals[i].residual * individuals[i].residual;
            n_valid++;
        }
    }
    
    if (n_valid == 0) {
        fprintf(stderr, "Error: No valid residuals for outlier detection\n");
        return 0;
    }
    
    double mean = sum / n_valid;
    double variance = (sum_sq - sum * sum / n_valid) / (n_valid - 1);
    double sd = sqrt(variance);
    
    if (!isfinite(mean) || !isfinite(sd) || sd <= 0.0) {
        fprintf(stderr, "Error: Invalid residual statistics (mean=%.6f, sd=%.6f)\n", mean, sd);
        return 0;
    }
    
    printf("Residual mean: %.6f, SD: %.6f (n=%u)\n", mean, sd, n_valid);
    
    /* Remove outliers */
    uint32_t removed = 0;
    for (uint32_t i = 0; i < count; i++) {
        if (individuals[i].valid && isfinite(individuals[i].residual)) {
            if (fabs(individuals[i].residual - mean) > 5.0 * sd) {
                individuals[i].valid = 0;
                removed++;
            }
        }
    }
    
    return removed;
}

/* Normalize residuals by sex */
static void normalize_by_sex(individual_t* individuals, uint32_t count,
                           const uint32_t* covar_indices, uint32_t n_covar_indices)
{
    /* Find sex covariate (assume it's covariate index 0 if present) */
    int sex_covar_idx = -1;
    if (n_covar_indices > 0) {
        sex_covar_idx = covar_indices[0];
    }
    
    if (sex_covar_idx < 0) {
        /* No sex covariate, normalize all together */
        double sum = 0.0, sum_sq = 0.0;
        uint32_t n_valid = 0;
        
        for (uint32_t i = 0; i < count; i++) {
            if (individuals[i].valid && isfinite(individuals[i].residual)) {
                sum += individuals[i].residual;
                sum_sq += individuals[i].residual * individuals[i].residual;
                n_valid++;
            }
        }
        
        if (n_valid == 0) return;
        
        double mean = sum / n_valid;
        double variance = (sum_sq - sum * sum / n_valid) / (n_valid - 1);
        double sd = sqrt(variance);
        
        if (!isfinite(sd) || sd <= 1e-10) {
            fprintf(stderr, "Error: Invalid overall SD: %.6f\n", sd);
            return;
        }
        
        /* Normalize */
        for (uint32_t i = 0; i < count; i++) {
            if (individuals[i].valid && isfinite(individuals[i].residual)) {
                individuals[i].residual = (individuals[i].residual - mean) / sd;
            }
        }
        
        printf("Overall normalization: mean=%.6f, sd=%.6f (n=%u)\n", mean, sd, n_valid);
        return;
    }
    
    /* Separate by sex (assume 0=male, 1=female) */
    double sum_male = 0.0, sum_sq_male = 0.0, sum_female = 0.0, sum_sq_female = 0.0;
    uint32_t n_male = 0, n_female = 0;
    
    for (uint32_t i = 0; i < count; i++) {
        if (individuals[i].valid && isfinite(individuals[i].residual)) {
            double sex_val = individuals[i].covars[sex_covar_idx];
            if (fabs(sex_val - 0.0) < 0.5) { /* Male */
                sum_male += individuals[i].residual;
                sum_sq_male += individuals[i].residual * individuals[i].residual;
                n_male++;
            } else if (fabs(sex_val - 1.0) < 0.5) { /* Female */
                sum_female += individuals[i].residual;
                sum_sq_female += individuals[i].residual * individuals[i].residual;
                n_female++;
            }
        }
    }
    
    /* Calculate sex-specific statistics */
    double mean_male = (n_male > 0) ? sum_male / n_male : 0.0;
    double mean_female = (n_female > 0) ? sum_female / n_female : 0.0;
    
    double sd_male = 1.0, sd_female = 1.0;
    if (n_male > 1) {
        double var_male = (sum_sq_male - sum_male * sum_male / n_male) / (n_male - 1);
        sd_male = sqrt(var_male);
    }
    if (n_female > 1) {
        double var_female = (sum_sq_female - sum_female * sum_female / n_female) / (n_female - 1);
        sd_female = sqrt(var_female);
    }
    
    if (!isfinite(sd_male) || sd_male <= 1e-10) {
        fprintf(stderr, "Warning: Invalid male SD, using 1.0\n");
        sd_male = 1.0;
    }
    if (!isfinite(sd_female) || sd_female <= 1e-10) {
        fprintf(stderr, "Warning: Invalid female SD, using 1.0\n");  
        sd_female = 1.0;
    }
    
    printf("Male normalization: mean=%.6f, sd=%.6f (n=%u)\n", mean_male, sd_male, n_male);
    printf("Female normalization: mean=%.6f, sd=%.6f (n=%u)\n", mean_female, sd_female, n_female);
    
    /* Apply sex-specific normalization */
    uint32_t normalized = 0;
    for (uint32_t i = 0; i < count; i++) {
        if (individuals[i].valid && isfinite(individuals[i].residual)) {
            double sex_val = individuals[i].covars[sex_covar_idx];
            if (fabs(sex_val - 0.0) < 0.5) { /* Male */
                individuals[i].residual = (individuals[i].residual - mean_male) / sd_male;
                normalized++;
            } else if (fabs(sex_val - 1.0) < 0.5) { /* Female */
                individuals[i].residual = (individuals[i].residual - mean_female) / sd_female;
                normalized++;
            }
        }
    }
    
    printf("Normalized %u individuals by sex\n", normalized);
}

/* Write output file with header inheritance */
static int write_output(const char* output_fname, individual_t* individuals, uint32_t count,
                       const char* pheno_name)
{
    FILE* f = fopen(output_fname, "w");
    if (!f) {
        perror(output_fname);
        return 0;
    }
    
    /* Write header with .resid suffix */
    if (pheno_name) {
        fprintf(f, "IID\t%s.resid\n", pheno_name);
    } else {
        fprintf(f, "IID\tresidual\n");
    }
    
    /* Write data for valid individuals only */
    for (uint32_t i = 0; i < count; i++) {
        if (individuals[i].valid && isfinite(individuals[i].residual)) {
            fprintf(f, "%s\t%.6f\n", individuals[i].iid, individuals[i].residual);
        }
    }
    
    fclose(f);
    return 1;
}

int main(int argc, char* argv[])
{
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <pheno_file> <covar_file> <covar_indices> <out_file>\n", argv[0]);
        fprintf(stderr, "  covar_indices: comma-separated 1-based indices (e.g., '1,2,3,4')\n");
        return 1;
    }
    
    const char* pheno_fname = argv[1];
    const char* covar_fname = argv[2];
    const char* indices_str = argv[3];
    const char* output_fname = argv[4];
    
    /* Parse covariate indices */
    uint32_t n_covar_indices;
    uint32_t* covar_indices = parse_indices(indices_str, &n_covar_indices);
    
    if (!covar_indices || n_covar_indices == 0) {
        fprintf(stderr, "Error: No valid covariate indices provided\n");
        return 1;
    }
    
    printf("Using covariate columns:");
    for (uint32_t i = 0; i < n_covar_indices; i++) {
        printf(" %u", covar_indices[i] + 1);
    }
    printf("\n");
    
    /* Read covariate file */
    uint32_t count;
    char** covar_names;
    uint32_t n_covars;
    individual_t* individuals = read_covariate_file(covar_fname, &count, &covar_names, &n_covars);
    
    if (!individuals) {
        free(covar_indices);
        return 1;
    }
    
    printf("Read %u individuals with %u covariates from %s\n", count, n_covars - 1, covar_fname);
    
    /* Validate covariate indices */
    for (uint32_t i = 0; i < n_covar_indices; i++) {
        if (covar_indices[i] >= n_covars - 1) {
            fprintf(stderr, "Error: Covariate index %u out of range (max: %u)\n", 
                    covar_indices[i] + 1, n_covars - 1);
            free(covar_indices);
            return 1;
        }
    }
    
    /* Read phenotype file */
    char* pheno_name;
    uint32_t matched = read_phenotype_file(pheno_fname, individuals, count, &pheno_name);
    printf("Matched %u phenotype values from %s\n", matched, pheno_fname);
    
    if (matched == 0) {
        fprintf(stderr, "Error: No phenotype matches found\n");
        free(covar_indices);
        return 1;
    }
    
    /* Apply complete case analysis */
    uint32_t complete_cases = apply_complete_case_analysis(individuals, count, covar_indices, n_covar_indices);
    printf("Using %u complete cases for regression\n", complete_cases);
    
    if (complete_cases <= n_covar_indices + 1) {
        fprintf(stderr, "Error: Not enough complete cases for regression (%u cases, need > %u)\n",
                complete_cases, n_covar_indices + 1);
        free(covar_indices);
        return 1;
    }
    
    /* Perform regression */
    double* beta = malloc((n_covar_indices + 1) * sizeof(double));
    uint32_t n_complete;
    
    if (!perform_regression(individuals, count, covar_indices, n_covar_indices, beta, &n_complete)) {
        fprintf(stderr, "Error: Regression failed\n");
        free(covar_indices);
        free(beta);
        return 1;
    }
    
    printf("Regression coefficients: intercept=%.6f", beta[0]);
    for (uint32_t i = 0; i < n_covar_indices; i++) {
        printf(", coeff%u=%.6f", i + 1, beta[i + 1]);
    }
    printf("\n");
    
    /* Calculate residuals */
    calculate_residuals(individuals, count, covar_indices, n_covar_indices, beta);
    
    /* Remove outliers */
    uint32_t outliers_removed = remove_outliers(individuals, count);
    printf("Removed %u outliers (>5 SD from mean)\n", outliers_removed);
    
    /* Normalize by sex */
    normalize_by_sex(individuals, count, covar_indices, n_covar_indices);
    
    /* Write output */
    if (!write_output(output_fname, individuals, count, pheno_name)) {
        fprintf(stderr, "Error: Failed to write output\n");
        free(covar_indices);
        free(beta);
        return 1;
    }
    
    printf("Output written to %s\n", output_fname);
    
    /* Cleanup */
    for (uint32_t i = 0; i < count; i++) {
        free(individuals[i].iid);
        free(individuals[i].covars);
    }
    free(individuals);
    free(covar_indices);
    free(beta);
    if (pheno_name) free(pheno_name);
    for (uint32_t i = 0; i < n_covars; i++) {
        free(covar_names[i]);
    }
    free(covar_names);
    
    return 0;
}
