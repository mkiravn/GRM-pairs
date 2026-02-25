/*
 * regress_y.c -- Phenotype regression, outlier removal, and normalization
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
    int valid; /* 1 if valid after outlier removal */
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
                if (*endptr == '\0') {
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
        
        double pheno = strtod(token, NULL);
        
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

/* Multiple linear regression using normal equations */
static int fit_regression(individual_t* individuals, uint32_t n,
                          const uint32_t* covar_indices, uint32_t n_selected_covars)
{
    /* Count valid observations */
    uint32_t n_valid = 0;
    for (uint32_t i = 0; i < n; i++) {
        if (!isnan(individuals[i].pheno)) {
            int valid = 1;
            for (uint32_t c = 0; c < n_selected_covars; c++) {
                if (isnan(individuals[i].covars[covar_indices[c]])) {
                    valid = 0;
                    break;
                }
            }
            if (valid) n_valid++;
        }
    }
    
    if (n_valid == 0) {
        fprintf(stderr, "Error: No valid observations for regression\n");
        return 1;
    }
    
    if (n_valid <= n_selected_covars + 1) {
        fprintf(stderr, "Error: Too few observations (%u) for %u parameters\n", 
                n_valid, n_selected_covars + 1);
        return 1;
    }
    
    fprintf(stderr, "Fitting regression with %u valid observations, %u covariates\n",
            n_valid, n_selected_covars);
    
    /* Build design matrix X and response vector y */
    double* X = malloc(n_valid * (n_selected_covars + 1) * sizeof(double));
    double* y = malloc(n_valid * sizeof(double));
    
    if (!X || !y) {
        free(X); free(y);
        fprintf(stderr, "Error: Memory allocation failed\n");
        return 1;
    }
    
    uint32_t row = 0;
    for (uint32_t i = 0; i < n; i++) {
        if (!isnan(individuals[i].pheno)) {
            int valid = 1;
            for (uint32_t c = 0; c < n_selected_covars; c++) {
                if (isnan(individuals[i].covars[covar_indices[c]])) {
                    valid = 0;
                    break;
                }
            }
            if (valid) {
                X[row * (n_selected_covars + 1)] = 1.0; /* Intercept */
                for (uint32_t c = 0; c < n_selected_covars; c++) {
                    X[row * (n_selected_covars + 1) + c + 1] = 
                        individuals[i].covars[covar_indices[c]];
                }
                y[row] = individuals[i].pheno;
                row++;
            }
        }
    }
    
    /* Solve normal equations: (X'X)β = X'y */
    uint32_t p = n_selected_covars + 1; /* Number of parameters */
    double* XtX = calloc(p * p, sizeof(double));
    double* Xty = calloc(p, sizeof(double));
    double* beta = calloc(p, sizeof(double));
    
    if (!XtX || !Xty || !beta) {
        free(X); free(y); free(XtX); free(Xty); free(beta);
        fprintf(stderr, "Error: Memory allocation failed\n");
        return 1;
    }
    
    /* Compute X'X */
    for (uint32_t i = 0; i < p; i++) {
        for (uint32_t j = 0; j < p; j++) {
            for (uint32_t k = 0; k < n_valid; k++) {
                XtX[i * p + j] += X[k * p + i] * X[k * p + j];
            }
        }
    }
    
    /* Compute X'y */
    for (uint32_t i = 0; i < p; i++) {
        for (uint32_t k = 0; k < n_valid; k++) {
            Xty[i] += X[k * p + i] * y[k];
        }
    }
    
    /* Solve normal equations using Gaussian elimination */
    /* Create augmented matrix [XtX | Xty] */
    double* aug_matrix = malloc(p * (p + 1) * sizeof(double));
    if (!aug_matrix) {
        free(X); free(y); free(XtX); free(Xty); free(beta);
        fprintf(stderr, "Error: Memory allocation failed\n");
        return 1;
    }
    
    /* Fill augmented matrix */
    for (uint32_t i = 0; i < p; i++) {
        for (uint32_t j = 0; j < p; j++) {
            aug_matrix[i * (p + 1) + j] = XtX[i * p + j];
        }
        aug_matrix[i * (p + 1) + p] = Xty[i]; /* Right-hand side */
    }
    
    /* Forward elimination */
    for (uint32_t i = 0; i < p; i++) {
        /* Find pivot */
        uint32_t max_row = i;
        for (uint32_t k = i + 1; k < p; k++) {
            if (fabs(aug_matrix[k * (p + 1) + i]) > fabs(aug_matrix[max_row * (p + 1) + i])) {
                max_row = k;
            }
        }
        
        /* Check for singular matrix */
        if (fabs(aug_matrix[max_row * (p + 1) + i]) < 1e-12) {
            fprintf(stderr, "Error: Singular matrix in regression\n");
            free(X); free(y); free(XtX); free(Xty); free(beta); free(aug_matrix);
            return 1;
        }
        
        /* Swap rows */
        if (max_row != i) {
            for (uint32_t j = 0; j <= p; j++) {
                double temp = aug_matrix[i * (p + 1) + j];
                aug_matrix[i * (p + 1) + j] = aug_matrix[max_row * (p + 1) + j];
                aug_matrix[max_row * (p + 1) + j] = temp;
            }
        }
        
        /* Make diagonal element 1 */
        double pivot = aug_matrix[i * (p + 1) + i];
        for (uint32_t j = 0; j <= p; j++) {
            aug_matrix[i * (p + 1) + j] /= pivot;
        }
        
        /* Eliminate column */
        for (uint32_t k = i + 1; k < p; k++) {
            double factor = aug_matrix[k * (p + 1) + i];
            for (uint32_t j = 0; j <= p; j++) {
                aug_matrix[k * (p + 1) + j] -= factor * aug_matrix[i * (p + 1) + j];
            }
        }
    }
    
    /* Back substitution */
    for (int i = (int)p - 1; i >= 0; i--) {
        beta[i] = aug_matrix[i * (p + 1) + p];
        for (uint32_t j = i + 1; j < p; j++) {
            beta[i] -= aug_matrix[i * (p + 1) + j] * beta[j];
        }
    }
    
    free(aug_matrix);
    
    /* Validate beta coefficients */
    for (uint32_t i = 0; i < p; i++) {
        if (!isfinite(beta[i])) {
            fprintf(stderr, "Error: Non-finite regression coefficient beta[%u] = %.6f\n", i, beta[i]);
            free(X); free(y); free(XtX); free(Xty); free(beta);
            return 1;
        }
    }
    
    fprintf(stderr, "Regression coefficients: intercept=%.6f", beta[0]);
    for (uint32_t i = 1; i < p; i++) {
        fprintf(stderr, ", coeff%u=%.6f", i, beta[i]);
    }
    fprintf(stderr, "\n");
    
    /* Calculate residuals */
    row = 0;
    for (uint32_t i = 0; i < n; i++) {
        if (!isnan(individuals[i].pheno)) {
            int valid = 1;
            for (uint32_t c = 0; c < n_selected_covars; c++) {
                if (isnan(individuals[i].covars[covar_indices[c]])) {
                    valid = 0;
                    break;
                }
            }
            if (valid) {
                double predicted = beta[0]; /* Intercept */
                for (uint32_t c = 0; c < n_selected_covars; c++) {
                    predicted += beta[c + 1] * individuals[i].covars[covar_indices[c]];
                }
                individuals[i].residual = individuals[i].pheno - predicted;
                row++;
            } else {
                individuals[i].residual = 0.0/0.0; /* NaN */
            }
        } else {
            individuals[i].residual = 0.0/0.0; /* NaN */
        }
    }
    
    free(X); free(y); free(XtX); free(Xty); free(beta);
    return 0;
}

/* Remove outliers (> 5 SD from mean) */
static void remove_outliers(individual_t* individuals, uint32_t n)
{
    /* Calculate mean and SD of residuals */
    double sum = 0.0, sum_sq = 0.0;
    uint32_t count = 0;
    
    /* First pass: calculate mean and count valid residuals */
    for (uint32_t i = 0; i < n; i++) {
        if (!isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            sum += individuals[i].residual;
            count++;
        }
    }
    
    if (count == 0) {
        fprintf(stderr, "Warning: No valid residuals for outlier removal\n");
        return;
    }
    
    if (count < 3) {
        fprintf(stderr, "Warning: Too few residuals (%u) for reliable outlier detection\n", count);
        return;
    }
    
    double mean = sum / count;
    
    /* Second pass: calculate sum of squares */
    for (uint32_t i = 0; i < n; i++) {
        if (!isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            double diff = individuals[i].residual - mean;
            sum_sq += diff * diff;
        }
    }
    
    double variance = sum_sq / (count - 1);
    double sd = sqrt(variance);
    
    /* Check for degenerate cases */
    if (!isfinite(mean) || !isfinite(sd) || sd <= 0.0) {
        fprintf(stderr, "Warning: Invalid statistics (mean=%.6f, sd=%.6f), skipping outlier removal\n", mean, sd);
        return;
    }
    
    fprintf(stderr, "Residual mean: %.6f, SD: %.6f (n=%u)\n", mean, sd, count);
    
    /* Mark outliers */
    uint32_t n_outliers = 0;
    for (uint32_t i = 0; i < n; i++) {
        if (!isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
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
    /* Process males (sex = 0) */
    double sum_male = 0.0, sum_sq_male = 0.0;
    uint32_t count_male = 0;
    
    for (uint32_t i = 0; i < n; i++) {
        if (individuals[i].sex == 0 && individuals[i].valid && 
            !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            sum_male += individuals[i].residual;
            count_male++;
        }
    }
    
    double mean_male = 0.0, sd_male = 1.0;
    if (count_male > 1) {
        mean_male = sum_male / count_male;
        
        /* Second pass for variance calculation */
        for (uint32_t i = 0; i < n; i++) {
            if (individuals[i].sex == 0 && individuals[i].valid && 
                !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
                double diff = individuals[i].residual - mean_male;
                sum_sq_male += diff * diff;
            }
        }
        
        double var_male = sum_sq_male / (count_male - 1);
        sd_male = sqrt(var_male);
        
        /* Check for degenerate cases */
        if (!isfinite(sd_male) || sd_male <= 1e-10) {
            fprintf(stderr, "Warning: Invalid male SD (%.6f), using SD=1.0\n", sd_male);
            sd_male = 1.0;
        }
    } else if (count_male == 1) {
        fprintf(stderr, "Warning: Only 1 male observation, using mean=0.0, SD=1.0\n");
        mean_male = 0.0;
        sd_male = 1.0;
    }
    
    /* Process females (sex = 1) */
    double sum_female = 0.0, sum_sq_female = 0.0;
    uint32_t count_female = 0;
    
    for (uint32_t i = 0; i < n; i++) {
        if (individuals[i].sex == 1 && individuals[i].valid && 
            !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            sum_female += individuals[i].residual;
            count_female++;
        }
    }
    
    double mean_female = 0.0, sd_female = 1.0;
    if (count_female > 1) {
        mean_female = sum_female / count_female;
        
        /* Second pass for variance calculation */
        for (uint32_t i = 0; i < n; i++) {
            if (individuals[i].sex == 1 && individuals[i].valid && 
                !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
                double diff = individuals[i].residual - mean_female;
                sum_sq_female += diff * diff;
            }
        }
        
        double var_female = sum_sq_female / (count_female - 1);
        sd_female = sqrt(var_female);
        
        /* Check for degenerate cases */
        if (!isfinite(sd_female) || sd_female <= 1e-10) {
            fprintf(stderr, "Warning: Invalid female SD (%.6f), using SD=1.0\n", sd_female);
            sd_female = 1.0;
        }
    } else if (count_female == 1) {
        fprintf(stderr, "Warning: Only 1 female observation, using mean=0.0, SD=1.0\n");
        mean_female = 0.0;
        sd_female = 1.0;
    }
    
    fprintf(stderr, "Male normalization: mean=%.6f, sd=%.6f (n=%u)\n", 
            mean_male, sd_male, count_male);
    fprintf(stderr, "Female normalization: mean=%.6f, sd=%.6f (n=%u)\n", 
            mean_female, sd_female, count_female);
    
    /* Apply normalization */
    uint32_t normalized_count = 0;
    for (uint32_t i = 0; i < n; i++) {
        if (individuals[i].valid && !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            if (individuals[i].sex == 0) {
                individuals[i].residual = (individuals[i].residual - mean_male) / sd_male;
                normalized_count++;
            } else if (individuals[i].sex == 1) {
                individuals[i].residual = (individuals[i].residual - mean_female) / sd_female;
                normalized_count++;
            } else {
                /* Unknown sex - mark as invalid */
                fprintf(stderr, "Warning: Unknown sex for individual %s, marking as invalid\n", individuals[i].iid);
                individuals[i].valid = 0;
                individuals[i].residual = 0.0/0.0;
            }
        }
    }
    
    fprintf(stderr, "Normalized %u individuals by sex\n", normalized_count);
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
        if (individuals[i].valid && !isnan(individuals[i].residual)) {
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
        fprintf(stderr, "  covar_indices: comma-separated 1-based column indices (e.g., '1,2,3')\n");
        return 1;
    }
    
    const char* pheno_fname = argv[1];
    const char* covar_fname = argv[2];
    const char* indices_str = argv[3];
    const char* out_fname = argv[4];
    
    /* Parse covariate indices */
    uint32_t n_selected_covars;
    uint32_t* covar_indices = parse_indices(indices_str, &n_selected_covars);
    if (!covar_indices || n_selected_covars == 0) {
        fprintf(stderr, "Error: No valid covariate indices specified\n");
        return 1;
    }
    
    fprintf(stderr, "Using %u covariates: ", n_selected_covars);
    for (uint32_t i = 0; i < n_selected_covars; i++) {
        fprintf(stderr, "%u", covar_indices[i] + 1);
        if (i < n_selected_covars - 1) fprintf(stderr, ",");
    }
    fprintf(stderr, "\n");
    
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
    for (uint32_t i = 0; i < n_selected_covars; i++) {
        if (covar_indices[i] >= n_covars) {
            fprintf(stderr, "Error: Covariate index %u exceeds available columns (%u)\n",
                    covar_indices[i] + 1, n_covars);
            free(covar_indices);
            return 1;
        }
    }
    
    /* Read phenotype data */
    char* pheno_header = NULL;
    if (read_phenotypes(pheno_fname, individuals, n_individuals, &pheno_header)) {
        free(covar_indices);
        return 1;
    }
    
    /* Data quality diagnostics */
    uint32_t n_pheno_missing = 0, n_covar_missing = 0, n_complete = 0;
    uint32_t n_male = 0, n_female = 0, n_sex_unknown = 0;
    
    for (uint32_t i = 0; i < n_individuals; i++) {
        if (isnan(individuals[i].pheno)) n_pheno_missing++;
        
        int covar_missing = 0;
        for (uint32_t c = 0; c < n_selected_covars; c++) {
            if (isnan(individuals[i].covars[covar_indices[c]])) {
                covar_missing = 1;
                break;
            }
        }
        if (covar_missing) n_covar_missing++;
        
        if (!isnan(individuals[i].pheno) && !covar_missing) n_complete++;
        
        if (individuals[i].sex == 0) n_male++;
        else if (individuals[i].sex == 1) n_female++;
        else n_sex_unknown++;
    }
    
    fprintf(stderr, "Data quality summary:\n");
    fprintf(stderr, "  Total individuals: %u\n", n_individuals);
    fprintf(stderr, "  Missing phenotypes: %u (%.1f%%)\n", n_pheno_missing, 100.0 * n_pheno_missing / n_individuals);
    fprintf(stderr, "  Missing covariates: %u (%.1f%%)\n", n_covar_missing, 100.0 * n_covar_missing / n_individuals);
    fprintf(stderr, "  Complete cases: %u (%.1f%%)\n", n_complete, 100.0 * n_complete / n_individuals);
    fprintf(stderr, "  Sex distribution: %u male, %u female, %u unknown\n", n_male, n_female, n_sex_unknown);
    
    if (n_complete < 10) {
        fprintf(stderr, "Error: Too few complete cases (%u) for reliable analysis\n", n_complete);
        free(covar_indices);
        return 1;
    }
    
    /* Fit regression */
    if (fit_regression(individuals, n_individuals, covar_indices, n_selected_covars)) {
        free(covar_indices);
        return 1;
    }
    
    /* Remove outliers */
    remove_outliers(individuals, n_individuals);
    
    /* Normalize by sex */
    normalize_by_sex(individuals, n_individuals);
    
    /* Final verification: calculate overall mean and SD after normalization */
    double final_sum = 0.0, final_sum_sq = 0.0;
    uint32_t final_count = 0;
    
    for (uint32_t i = 0; i < n_individuals; i++) {
        if (individuals[i].valid && !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
            final_sum += individuals[i].residual;
            final_count++;
        }
    }
    
    if (final_count > 1) {
        double final_mean = final_sum / final_count;
        
        for (uint32_t i = 0; i < n_individuals; i++) {
            if (individuals[i].valid && !isnan(individuals[i].residual) && isfinite(individuals[i].residual)) {
                double diff = individuals[i].residual - final_mean;
                final_sum_sq += diff * diff;
            }
        }
        
        double final_var = final_sum_sq / (final_count - 1);
        double final_sd = sqrt(final_var);
        
        fprintf(stderr, "Final normalized residuals: mean=%.6f, sd=%.6f (n=%u)\n", 
                final_mean, final_sd, final_count);
        
        if (fabs(final_mean) > 0.01 || fabs(final_sd - 1.0) > 0.01) {
            fprintf(stderr, "Warning: Normalization may not be perfect (expected mean≈0, sd≈1)\n");
        }
    } else {
        fprintf(stderr, "Warning: Too few final observations (%u) to verify normalization\n", final_count);
    }
    
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
