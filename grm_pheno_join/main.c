/*
 * grm_pheno_join.c -- Join GRM and phenotype cross-product results and perform binned analysis
 *
 * Usage: grm_pheno_join <grm_file> <pheno_cross_file> <out_file>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define MAX_LINE 4096
#define BIN_WIDTH 0.001
#define MAX_BINS 2000  /* For GRM values from -1.0 to 1.0 */

/* Structure to store paired data */
typedef struct {
    char* iid1;
    char* iid2;
    double grm_value;
    double* pheno_values;
    uint32_t n_phenos;
} pair_data_t;

/* Structure for binned statistics */
typedef struct {
    double bin_center;
    double* sum_pheno;
    double* sum_pheno_sq;
    uint32_t* count;
    uint32_t n_phenos;
} bin_stats_t;

/* Parse GRM file */
static pair_data_t* read_grm_file(const char* grm_fname, uint32_t* n_pairs_out)
{
    FILE* f = fopen(grm_fname, "r");
    if (!f) { 
        perror(grm_fname); 
        return NULL; 
    }

    char line[MAX_LINE];
    
    /* Count data lines (skip header if present) */
    uint32_t n_lines = 0;
    int has_header = 0;
    
    if (fgets(line, sizeof(line), f)) {
        /* Check if first line looks like a header */
        char line_copy[MAX_LINE];
        strcpy(line_copy, line);
        char* token = strtok(line_copy, " \t\n\r");
        if (token) {
            char* endptr;
            strtod(token, &endptr);
            if (*endptr != '\0') {
                has_header = 1;
            } else {
                n_lines = 1; /* First line is data */
            }
        }
    }
    
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '\0' && line[0] != '\n') n_lines++;
    }
    
    if (n_lines == 0) {
        fprintf(stderr, "Error: No data found in %s\n", grm_fname);
        fclose(f);
        return NULL;
    }
    
    /* Allocate array */
    pair_data_t* pairs = malloc(n_lines * sizeof(pair_data_t));
    if (!pairs) {
        fclose(f);
        return NULL;
    }
    
    /* Read data */
    rewind(f);
    if (has_header) {
        if (!fgets(line, sizeof(line), f)) {
            free(pairs);
            fclose(f);
            return NULL;
        }
    }
    
    uint32_t i = 0;
    while (fgets(line, sizeof(line), f) && i < n_lines) {
        if (line[0] == '\0' || line[0] == '\n') continue;
        
        /* Parse: IID1 IID2 GRM_value */
        char* token = strtok(line, " \t\n\r");
        if (!token) continue;
        
        pairs[i].iid1 = malloc(strlen(token) + 1);
        if (!pairs[i].iid1) {
            for (uint32_t j = 0; j < i; j++) {
                free(pairs[j].iid1);
                free(pairs[j].iid2);
            }
            free(pairs);
            fclose(f);
            return NULL;
        }
        strcpy(pairs[i].iid1, token);
        
        token = strtok(NULL, " \t\n\r");
        if (!token) {
            free(pairs[i].iid1);
            continue;
        }
        
        pairs[i].iid2 = malloc(strlen(token) + 1);
        if (!pairs[i].iid2) {
            free(pairs[i].iid1);
            for (uint32_t j = 0; j < i; j++) {
                free(pairs[j].iid1);
                free(pairs[j].iid2);
            }
            free(pairs);
            fclose(f);
            return NULL;
        }
        strcpy(pairs[i].iid2, token);
        
        token = strtok(NULL, " \t\n\r");
        if (!token) {
            free(pairs[i].iid1);
            free(pairs[i].iid2);
            continue;
        }
        
        pairs[i].grm_value = strtod(token, NULL);
        pairs[i].pheno_values = NULL;
        pairs[i].n_phenos = 0;
        
        i++;
    }
    
    fclose(f);
    *n_pairs_out = i;
    fprintf(stderr, "Read %u GRM pairs from %s\n", i, grm_fname);
    return pairs;
}

/* Join with phenotype cross-products */
static int join_pheno_data(pair_data_t* pairs, uint32_t n_pairs, 
                          const char* pheno_fname, uint32_t* n_phenos_out)
{
    FILE* f = fopen(pheno_fname, "r");
    if (!f) { 
        perror(pheno_fname); 
        return 1; 
    }

    char line[MAX_LINE];
    uint32_t n_phenos = 0;
    
    /* Read header to count phenotype columns */
    if (!fgets(line, sizeof(line), f)) {
        fprintf(stderr, "Error: Empty file %s\n", pheno_fname);
        fclose(f);
        return 1;
    }
    
    char line_copy[MAX_LINE];
    strcpy(line_copy, line);
    char* token = strtok(line_copy, " \t\n\r");
    while (token) {
        n_phenos++;
        token = strtok(NULL, " \t\n\r");
    }
    
    if (n_phenos < 3) {
        fprintf(stderr, "Error: Expected at least 3 columns (IID1, IID2, pheno), got %u\n", n_phenos);
        fclose(f);
        return 1;
    }
    
    n_phenos -= 2; /* Subtract IID1 and IID2 columns */
    *n_phenos_out = n_phenos;
    
    /* Allocate phenotype arrays for all pairs */
    for (uint32_t i = 0; i < n_pairs; i++) {
        pairs[i].pheno_values = calloc(n_phenos, sizeof(double));
        pairs[i].n_phenos = n_phenos;
        if (!pairs[i].pheno_values) {
            for (uint32_t j = 0; j < i; j++) {
                free(pairs[j].pheno_values);
            }
            fclose(f);
            return 1;
        }
        /* Initialize to NaN */
        for (uint32_t p = 0; p < n_phenos; p++) {
            pairs[i].pheno_values[p] = 0.0/0.0;
        }
    }
    
    /* Read phenotype data and match to pairs */
    uint32_t matched = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '\n') continue;
        
        /* Parse: IID1 IID2 pheno1 [pheno2 ...] */
        token = strtok(line, " \t\n\r");
        if (!token) continue;
        char* iid1 = token;
        
        token = strtok(NULL, " \t\n\r");
        if (!token) continue;
        char* iid2 = token;
        
        /* Find matching pair */
        for (uint32_t i = 0; i < n_pairs; i++) {
            if ((strcmp(pairs[i].iid1, iid1) == 0 && strcmp(pairs[i].iid2, iid2) == 0) ||
                (strcmp(pairs[i].iid1, iid2) == 0 && strcmp(pairs[i].iid2, iid1) == 0)) {
                
                /* Read phenotype values */
                for (uint32_t p = 0; p < n_phenos; p++) {
                    token = strtok(NULL, " \t\n\r");
                    if (token) {
                        char* endptr;
                        double val = strtod(token, &endptr);
                        if (*endptr == '\0' && !isnan(val)) {
                            pairs[i].pheno_values[p] = val;
                        }
                    }
                }
                matched++;
                break;
            }
        }
    }
    
    fclose(f);
    fprintf(stderr, "Matched %u pairs with phenotype data from %s\n", matched, pheno_fname);
    fprintf(stderr, "Found %u phenotype columns\n", n_phenos);
    return 0;
}

/* Perform binned analysis */
static int binned_analysis(const pair_data_t* pairs, uint32_t n_pairs, uint32_t n_phenos,
                          const char* out_fname)
{
    /* Initialize bins */
    int n_bins = (int)(2.0 / BIN_WIDTH) + 1; /* From -1.0 to 1.0 */
    bin_stats_t* bins = calloc(n_bins, sizeof(bin_stats_t));
    if (!bins) return 1;
    
    for (int b = 0; b < n_bins; b++) {
        bins[b].bin_center = -1.0 + b * BIN_WIDTH;
        bins[b].n_phenos = n_phenos;
        bins[b].sum_pheno = calloc(n_phenos, sizeof(double));
        bins[b].sum_pheno_sq = calloc(n_phenos, sizeof(double));
        bins[b].count = calloc(n_phenos, sizeof(uint32_t));
        
        if (!bins[b].sum_pheno || !bins[b].sum_pheno_sq || !bins[b].count) {
            for (int i = 0; i <= b; i++) {
                free(bins[i].sum_pheno);
                free(bins[i].sum_pheno_sq);
                free(bins[i].count);
            }
            free(bins);
            return 1;
        }
    }
    
    /* Accumulate data into bins */
    uint32_t total_processed = 0;
    for (uint32_t i = 0; i < n_pairs; i++) {
        if (isnan(pairs[i].grm_value)) continue;
        
        /* Find appropriate bin */
        int bin_idx = (int)((pairs[i].grm_value + 1.0) / BIN_WIDTH);
        if (bin_idx < 0) bin_idx = 0;
        if (bin_idx >= n_bins) bin_idx = n_bins - 1;
        
        /* Accumulate phenotype values */
        for (uint32_t p = 0; p < n_phenos; p++) {
            if (!isnan(pairs[i].pheno_values[p])) {
                double val = pairs[i].pheno_values[p];
                bins[bin_idx].sum_pheno[p] += val;
                bins[bin_idx].sum_pheno_sq[p] += val * val;
                bins[bin_idx].count[p]++;
            }
        }
        total_processed++;
    }
    
    fprintf(stderr, "Processed %u pairs into %d bins\n", total_processed, n_bins);
    
    /* Write results */
    FILE* f = fopen(out_fname, "w");
    if (!f) { 
        perror(out_fname);
        for (int b = 0; b < n_bins; b++) {
            free(bins[b].sum_pheno);
            free(bins[b].sum_pheno_sq);
            free(bins[b].count);
        }
        free(bins);
        return 1; 
    }
    
    /* Write header */
    fprintf(f, "grm_bin");
    for (uint32_t p = 0; p < n_phenos; p++) {
        fprintf(f, "\tpheno%u_mean\tpheno%u_n\tpheno%u_se", p+1, p+1, p+1);
    }
    fprintf(f, "\n");
    
    /* Write bin statistics */
    for (int b = 0; b < n_bins; b++) {
        /* Check if bin has any data */
        int has_data = 0;
        for (uint32_t p = 0; p < n_phenos; p++) {
            if (bins[b].count[p] > 0) {
                has_data = 1;
                break;
            }
        }
        
        if (!has_data) continue;
        
        fprintf(f, "%.3f", bins[b].bin_center);
        
        for (uint32_t p = 0; p < n_phenos; p++) {
            if (bins[b].count[p] > 0) {
                double mean = bins[b].sum_pheno[p] / bins[b].count[p];
                double se = 0.0;
                
                if (bins[b].count[p] > 1) {
                    double var = (bins[b].sum_pheno_sq[p] - 
                                 bins[b].sum_pheno[p] * bins[b].sum_pheno[p] / bins[b].count[p]) 
                                 / (bins[b].count[p] - 1);
                    se = sqrt(var / bins[b].count[p]);
                }
                
                fprintf(f, "\t%.6f\t%u\t%.6f", mean, bins[b].count[p], se);
            } else {
                fprintf(f, "\tNA\t0\tNA");
            }
        }
        fprintf(f, "\n");
    }
    
    fclose(f);
    
    /* Cleanup */
    for (int b = 0; b < n_bins; b++) {
        free(bins[b].sum_pheno);
        free(bins[b].sum_pheno_sq);
        free(bins[b].count);
    }
    free(bins);
    
    return 0;
}

int main(int argc, char* argv[])
{
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <grm_file> <pheno_cross_file> <out_file>\n", argv[0]);
        fprintf(stderr, "  Bins GRM values in %.3f increments and calculates phenotype statistics\n", BIN_WIDTH);
        return 1;
    }
    
    const char* grm_fname = argv[1];
    const char* pheno_fname = argv[2];
    const char* out_fname = argv[3];
    
    /* Read GRM data */
    uint32_t n_pairs;
    pair_data_t* pairs = read_grm_file(grm_fname, &n_pairs);
    if (!pairs) return 1;
    
    /* Join with phenotype data */
    uint32_t n_phenos;
    if (join_pheno_data(pairs, n_pairs, pheno_fname, &n_phenos)) {
        for (uint32_t i = 0; i < n_pairs; i++) {
            free(pairs[i].iid1);
            free(pairs[i].iid2);
            free(pairs[i].pheno_values);
        }
        free(pairs);
        return 1;
    }
    
    /* Perform binned analysis */
    if (binned_analysis(pairs, n_pairs, n_phenos, out_fname)) {
        for (uint32_t i = 0; i < n_pairs; i++) {
            free(pairs[i].iid1);
            free(pairs[i].iid2);
            free(pairs[i].pheno_values);
        }
        free(pairs);
        return 1;
    }
    
    fprintf(stderr, "Binned analysis complete. Output written to %s\n", out_fname);
    
    /* Cleanup */
    for (uint32_t i = 0; i < n_pairs; i++) {
        free(pairs[i].iid1);
        free(pairs[i].iid2);
        free(pairs[i].pheno_values);
    }
    free(pairs);
    
    return 0;
}
