/*
 * regress_y.c -- Simple phenotype regression with complete case analysis
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
    int valid; /* 1 if valid after complete case analysis */
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
    
    /* Allocate space for covariate names (skip first column which is IID) */
    n_covars--; /* Subtract 1 for IID column */
    covar_names = malloc(n_covars * sizeof(char*));
    if (!covar_names) { fclose(f); return NULL; }
    
    /* Extract column names */
    strcpy(line_copy, line);
    token = strtok(line_copy, " \t\n\r"); /* Skip IID column */
    for (uint32_t i = 0; i < n_covars; i++) {
        token = strtok(NULL, " \t\n\r");
        if (token) {
            covar_names[i] = malloc(strlen(token) + 1);
            if (!covar_names[i]) {
                for (uint32_t j = 0; j < i; j++) free(covar_names[j]);
                free(covar_names); fclose(f); return NULL;
            }
            strcpy(covar_names[i], token);
        } else {
            char default_name[32];
            snprintf(default_name, sizeof(default_name), "cov%u", i+1);
            covar_names[i] = malloc(strlen(default_name) + 1);
            strcpy(covar_names[i], default_name);
        }
    }

    /* Count data lines */
    uint32_t n = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '\0' && line[0] != '\n') n++;
    }
    rewind(f);
    if (!fgets(line, sizeof(line), f)) {
        fprintf(stderr, "Error: Unable to skip header in %s\n", covar_fname);
        fclose(f);
        for (uint32_t i = 0; i < n_covars; i++) free(covar_names[i]);
        free(covar_names);
        return NULL;
    } /* Skip header */

    individual_t* individuals = malloc(n * sizeof(individual_t));
    if (!individuals) { 
        for (uint32_t i = 0; i < n_covars; i++) free(covar_names[i]);
        free(covar_names); 
        fclose(f); 
        return NULL; 
    }

    uint32_t i = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '\n') continue;
        
        /* Parse IID */
        token = strtok(line, " \t\n\r");
        if (!token) continue;
        
        individuals[i].iid = malloc(strlen(token) + 1);
        if (!individuals[i].iid) {
            for (uint32_t j = 0; j < i; j++) {
                free(individuals[j].iid);
                free(individuals[j].covars);
            }
            free(individuals);
            for (uint32_t j = 0; j < n_covars; j++) free(covar_names[j]);
            free(covar_names);
            fclose(f);
            return NULL;
        }
        strcpy(individuals[i].iid, token);
        
        /* Allocate and parse covariates */
        individuals[i].covars = malloc(n_covars * sizeof(double));
        individuals[i].n_covars = n_covars;
        individuals[i].pheno = 0.0/0.0; /* NaN - will be filled later */
        individuals[i].residual = 0.0;
        individuals[i].sex = -1; /* Unknown initially */
        individuals[i].valid = 1;
        
        if (!individuals[i].covars) {
            free(individuals[i].iid);
            for (uint32_t j = 0; j < i; j++) {
                free(individuals[j].iid);
                free(individuals[j].covars);
            }
            free(individuals);
            for (uint32_t j = 0; j < n_covars; j++) free(covar_names[j]);
            free(covar_names);
            fclose(f);
            return NULL;
        }
        
        /* Parse all covariate values */
        for (uint32_t c = 0; c < n_covars; c++) {
            token = strtok(NULL, " \t\n\r");
            if (token) {
                char* endptr;
                double value = strtod(token, &endptr);
                if (*endptr == '\0' && isfinite(value)) {
                    individuals[i].covars[c] = value;
                    /* Set sex from first covariate */
                    if (c == 0) {
                        if (value == 0.0) individuals[i].sex = 0; /* Male */
                        else if (value == 1.0) individuals[i].sex = 1; /* Female */
                    }
                } else {
                    individuals[i].covars[c] = 0.0/0.0; /* NaN for non-numeric */
                }
            } else {
                individuals[i].covars[c] = 0.0/0.0; /* NaN for missing */
            }
        }
        
        i++;
    }
    fclose(f);

    *count_out = i;
    *covar_names_out = covar_names;
    *n_covars_out = n_covars;
    fprintf(stderr, "Read %u individuals with %u covariates from %s\n", 
            i, n_covars, covar_fname);
    return individuals;
}

/* Read phenotype values and match to individuals */
static int read_phenotypes(const char* pheno_fname, individual_t* individuals, uint32_t n,
                          char** pheno_header_out)
{
    FILE* f = fopen(pheno_fname, "r");
    if (!f) { perror(pheno_fname); return 1; }

    char line[2048];
    uint32_t count = 0;
    char* pheno_header = NULL;
    
    /* Read first line to check for header */
    if (fgets(line, sizeof(line), f)) {
        /* Check if first line looks like a header */
        char line_copy[2048];
        strcpy(line_copy, line);
        char* token = strtok(line_copy, " \t\n\r");
        if (token) {
            char* endptr;
            double test_val = strtod(token, &endptr);
            /* If first token is non-numeric, assume it's a header */
            if (*endptr != '\0' || isnan(test_val)) {
                /* Extract phenotype column name (second column) */
                token = strtok(NULL, " \t\n\r");
                if (token) {
                    pheno_header = malloc(strlen(token) + 1);
                    if (pheno_header) {
                        strcpy(pheno_header, token);
                    }
                }
                /* Skip this header line */
            } else {
                /* This is data, rewind to read it */
                rewind(f);
            }
        }
    }
    
    /* If no header found, use default */
    if (!pheno_header) {
        pheno_header = malloc(6);
        if (pheno_header) {
            strcpy(pheno_header, "pheno");
        }
    }

    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '\n') continue;
        
        char* token = strtok(line, " \t\n\r");
        if (!token) continue;
        
        char* iid = token;
        token = strtok(NULL, " \t\n\r");
        if (!token) continue;
        
        char* endptr;
        double pheno = strtod(token, &endptr);
        if (*endptr != '\0') continue; /* Skip non-numeric */
        
        /* Find matching individual */
        for (uint32_t i = 0; i < n; i++) {
            if (strcmp(individuals[i].iid, iid) == 0) {
                individuals[i].pheno = pheno;
                count++;
                break;
            }
        }
    }
    
    fclose(f);
    fprintf(stderr, "Matched %u phenotype values from %s\n", count, pheno_fname);
    *pheno_header_out = pheno_header;
    return 0;
}

/* Simple multiple linear regression with complete case analysis */
static int fit_regression(individual_t* individuals, uint32_t n,
                          const uint32_t* covar_indices, uint32_t n_selected_covars)
{
    /* Mark individuals with any missing data as invalid */
    uint32_t n_complete = 0;
    for (uint32_t i = 0; i < n; i++) {
        individuals[i].valid = 1; /* Start as valid */
        
        /* Check phenotype */
        if (isnan(individuals[i].pheno) || !isfinite(individuals[i].pheno)) {
            individuals[i].valid = 0;
        }
        
        /* Check covariates */
        if (individuals[i].valid) {
            for (uint32_t c = 0; c < n_selected_covars; c++) {
                if (isnan(individuals[i].covars[covar_indices[c]]) || 
                    !isfinite(individuals[i].covars[covar_indices[c]])) {
                    individuals[i].valid = 0;
                    break;
                }
            }
        }
        
        if (individuals[i].valid) {
            n_complete++;
        }
    }
    
    if (n_complete == 0) {
        fprintf(stderr, "Error: No individuals with complete data\n");
        return 1;
    }
    
    if (n_complete <= n_selected_covars + 1) {
        fprintf(stderr, "Error: Too few complete cases (%u) for %u parameters\n", 
                n_complete, n_selected_covars + 1);
        return 1;
    }
    
    fprintf(stderr, "Using %u complete cases for regression\n", n_complete);
    
    /* Simple matrix setup for complete cases only */
    uint32_t p = n_selected_covars + 1; /* Parameters including intercept */
    
    /* Use direct solution for small problems */
    if (p == 2) {
        /* Simple linear regression: y = a + b*x */
        double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0;
        
        for (uint32_t i = 0; i < n; i++) {
            if (individuals[i].valid) {
                double x = individuals[i].covars[covar_indices[0]];
                double y = individuals[i].pheno;
                sum_x += x;
                sum_y += y;
                sum_xy += x * y;
                sum_xx += x * x;
            }
        }
        
        double mean_x = sum_x / n_complete;
        double mean_y = sum_y / n_complete;
        double beta1 = (sum_xy - n_complete * mean_x * mean_y) / (sum_xx - n_complete * mean_x * mean_x);
        double beta0 = mean_y - beta1 * mean_x;
        
        fprintf(stderr, "Regression coefficients: intercept=%.6f, coeff1=%.6f\n", beta0, beta1);
        
        /* Calculate residuals */
        for (uint32_t i = 0; i < n; i++) {
            if (individuals[i].valid) {
                double fitted = beta0 + beta1 * individuals[i].covars[covar_indices[0]];
                individuals[i].residual = individuals[i].pheno - fitted;
            } else {
                individuals[i].residual = 0.0/0.0;
            }
        }
        
    } else {
        fprintf(stderr, "Error: Only simple linear regression (1 covariate) is supported\n");
        return 1;
    }
    
    return 0;
}

/* Remove outliers (> 5 SD from mean) */
static void remove_outliers(individual_t* individuals, uint32_t n)
{
    /* Calculate mean and SD of residuals for valid individuals only */
    double sum = 0.0, sum_sq = 0.0;
    uint32_t count = 0;
    
    for (uint32_t i = 0; i < n; i++) {
        if (individuals[i].valid && !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            sum += individuals[i].residual;
            count++;
        }
    }
    
    if (count < 3) {
        fprintf(stderr, "Warning: Too few valid residuals (%u) for outlier removal\n", count);
        return;
    }
    
    double mean = sum / count;
    
    for (uint32_t i = 0; i < n; i++) {
        if (individuals[i].valid && !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            double diff = individuals[i].residual - mean;
            sum_sq += diff * diff;
        }
    }
    
    double variance = sum_sq / (count - 1);
    double sd = sqrt(variance);
    
    if (!isfinite(sd) || sd <= 0.0) {
        fprintf(stderr, "Warning: Invalid SD (%.6f), skipping outlier removal\n", sd);
        return;
    }
    
    fprintf(stderr, "Residual mean: %.6f, SD: %.6f (n=%u)\n", mean, sd, count);
    
    /* Mark outliers as invalid */
    uint32_t n_outliers = 0;
    for (uint32_t i = 0; i < n; i++) {
        if (individuals[i].valid && !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            double z_score = fabs(individuals[i].residual - mean) / sd;
            if (z_score > 5.0) {
                individuals[i].valid = 0;
                n_outliers++;
            }
        }
    }
    
    fprintf(stderr, "Removed %u outliers (>5 SD from mean)\n", n_outliers);
}

/* Normalize residuals within each sex */
static void normalize_by_sex(individual_t* individuals, uint32_t n)
{
    /* Process males and females separately */
    double sum_male = 0.0, sum_female = 0.0;
    uint32_t count_male = 0, count_female = 0;
    
    for (uint32_t i = 0; i < n; i++) {
        if (individuals[i].valid && !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            if (individuals[i].sex == 0) {
                sum_male += individuals[i].residual;
                count_male++;
            } else if (individuals[i].sex == 1) {
                sum_female += individuals[i].residual;
                count_female++;
            }
        }
    }
    
    double mean_male = (count_male > 0) ? sum_male / count_male : 0.0;
    double mean_female = (count_female > 0) ? sum_female / count_female : 0.0;
    
    /* Calculate SD for each sex */
    double sum_sq_male = 0.0, sum_sq_female = 0.0;
    
    for (uint32_t i = 0; i < n; i++) {
        if (individuals[i].valid && !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            if (individuals[i].sex == 0) {
                double diff = individuals[i].residual - mean_male;
                sum_sq_male += diff * diff;
            } else if (individuals[i].sex == 1) {
                double diff = individuals[i].residual - mean_female;
                sum_sq_female += diff * diff;
            }
        }
    }
    
    double sd_male = (count_male > 1) ? sqrt(sum_sq_male / (count_male - 1)) : 1.0;
    double sd_female = (count_female > 1) ? sqrt(sum_sq_female / (count_female - 1)) : 1.0;
    
    if (sd_male <= 0.0) sd_male = 1.0;
    if (sd_female <= 0.0) sd_female = 1.0;
    
    fprintf(stderr, "Male normalization: mean=%.6f, sd=%.6f (n=%u)\n", mean_male, sd_male, count_male);
    fprintf(stderr, "Female normalization: mean=%.6f, sd=%.6f (n=%u)\n", mean_female, sd_female, count_female);
    
    /* Apply normalization */
    uint32_t normalized = 0;
    for (uint32_t i = 0; i < n; i++) {
        if (individuals[i].valid && !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            if (individuals[i].sex == 0) {
                individuals[i].residual = (individuals[i].residual - mean_male) / sd_male;
                normalized++;
            } else if (individuals[i].sex == 1) {
                individuals[i].residual = (individuals[i].residual - mean_female) / sd_female;
                normalized++;
            } else {
                /* Unknown sex - leave as is but mark invalid */
                individuals[i].valid = 0;
            }
        }
    }
    
    fprintf(stderr, "Normalized %u individuals by sex\n", normalized);
}

/* Write output */
static int write_output(const char* out_fname, const individual_t* individuals, uint32_t n,
                        const char* pheno_header)
{
    FILE* f = fopen(out_fname, "w");
    if (!f) { perror(out_fname); return 1; }
    
    /* Write header with inherited phenotype name + .resid suffix */
    fprintf(f, "IID\t%s.resid\n", pheno_header ? pheno_header : "pheno");
    
    for (uint32_t i = 0; i < n; i++) {
        if (individuals[i].valid && !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            fprintf(f, "%s\t%.6f\n", individuals[i].iid, individuals[i].residual);
        } else {
            fprintf(f, "%s\tNA\n", individuals[i].iid);
        }
    }
    
    fclose(f);
    return 0;
}

int main(int argc, char* argv[])
{
    if (argc != 5) {
        fprintf(stderr, "Usage: %s <pheno_file> <covar_file> <covar_indices> <out_file>\n", argv[0]);
        fprintf(stderr, "  covar_indices: comma-separated 1-based column indices (e.g., '1')\n");
        fprintf(stderr, "  NOTE: This version only supports simple linear regression (1 covariate)\n");
        return 1;
    }
    
    const char* pheno_fname = argv[1];
    const char* covar_fname = argv[2];
    const char* indices_str = argv[3];
    const char* out_fname = argv[4];
    
    /* Parse covariate indices */
    uint32_t n_selected_covars;
    uint32_t* covar_indices = parse_indices(indices_str, &n_selected_covars);
    if (!covar_indices || n_selected_covars != 1) {
        fprintf(stderr, "Error: This version only supports 1 covariate (simple linear regression)\n");
        if (covar_indices) free(covar_indices);
        return 1;
    }
    
    fprintf(stderr, "Using covariate column: %u\n", covar_indices[0] + 1);
    
    /* Read covariate data */
    uint32_t n_individuals;
    char** covar_names;
    uint32_t n_covars;
    individual_t* individuals = read_covariate_file(covar_fname, &n_individuals,
                                                   &covar_names, &n_covars);
    if (!individuals) {
        free(covar_indices);
        return 1;
    }
    
    /* Validate covariate indices */
    if (covar_indices[0] >= n_covars) {
        fprintf(stderr, "Error: Covariate index %u exceeds available columns (%u)\n",
                covar_indices[0] + 1, n_covars);
        free(covar_indices);
        return 1;
    }
    
    /* Read phenotype data */
    char* pheno_header = NULL;
    if (read_phenotypes(pheno_fname, individuals, n_individuals, &pheno_header)) {
        free(covar_indices);
        return 1;
    }
    
    /* Fit regression with complete case analysis */
    if (fit_regression(individuals, n_individuals, covar_indices, n_selected_covars)) {
        free(covar_indices);
        if (pheno_header) free(pheno_header);
        return 1;
    }
    
    /* Remove outliers */
    remove_outliers(individuals, n_individuals);
    
    /* Normalize by sex */
    normalize_by_sex(individuals, n_individuals);
    
    /* Write output */
    if (write_output(out_fname, individuals, n_individuals, pheno_header)) {
        free(covar_indices);
        if (pheno_header) free(pheno_header);
        return 1;
    }
    
    fprintf(stderr, "Output written to %s\n", out_fname);
    
    /* Cleanup */
    for (uint32_t i = 0; i < n_individuals; i++) {
        free(individuals[i].iid);
        free(individuals[i].covars);
    }
    free(individuals);
    for (uint32_t i = 0; i < n_covars; i++) {
        free(covar_names[i]);
    }
    free(covar_names);
    free(covar_indices);
    if (pheno_header) free(pheno_header);
    
    return 0;
}
