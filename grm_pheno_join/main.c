/*
 * grm_pheno_join.c -- Join GRM and phenotype cross-product results and perform binned analysis
 *
 * Usage: grm_pheno_join <grm_file> <pheno_cross_file> <out_file> [bin_width] [min_grm] [max_grm]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define MAX_LINE 4096
#define DEFAULT_BIN_WIDTH 0.001
#define DEFAULT_MIN_GRM -0.05
#define DEFAULT_MAX_GRM 1.3
#define MAX_BINS 5000  /* Enough for wider range */

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
    double grm_sum;      /* Sum of GRM values in this bin */
    double grm_sum_sq;   /* Sum of squared GRM values */
    double* sum_pheno;
    double* sum_pheno_sq;
    uint32_t* count;
    uint32_t grm_count;  /* Count of GRM values */
    uint32_t n_phenos;
} bin_stats_t;

/* Parse GRM file with robust column detection */
static pair_data_t* read_grm_file(const char* grm_fname, uint32_t* n_pairs_out)
{
    FILE* f = fopen(grm_fname, "r");
    if (!f) { 
        perror(grm_fname); 
        return NULL; 
    }

    char line[MAX_LINE];
    
    /* Detect file format by examining header */
    int has_header = 0;
    int grm_column = -1;  /* 0-based index of GRM column */
    
    if (fgets(line, sizeof(line), f)) {
        /* Check if first line looks like a header */
        char line_copy[MAX_LINE];
        strcpy(line_copy, line);
        
        char* token = strtok(line_copy, " \t\n\r");
        int col = 0;
        while (token) {
            /* Look for column headers that indicate GRM column */
            if (strcasecmp(token, "GRM") == 0 || strcasecmp(token, "KINSHIP") == 0) {
                grm_column = col;
                has_header = 1;
                break;
            }
            col++;
            token = strtok(NULL, " \t\n\r");
        }
        
        /* If no header found, assume format by number of columns */
        if (grm_column == -1) {
            strcpy(line_copy, line);
            token = strtok(line_copy, " \t\n\r");
            col = 0;
            while (token) {
                col++;
                token = strtok(NULL, " \t\n\r");
            }
            
            if (col == 3) {
                grm_column = 2;  /* IID1 IID2 GRM */
            } else if (col == 4) {
                grm_column = 3;  /* IID1 IID2 N_SNPs GRM */
            } else {
                fprintf(stderr, "Error: Unexpected number of columns (%d) in %s\n", col, grm_fname);
                fprintf(stderr, "Expected: IID1 IID2 [N_SNPs] GRM\n");
                fclose(f);
                return NULL;
            }
            has_header = 0; /* First line is data */
        }
    }
    
    if (grm_column == -1) {
        fprintf(stderr, "Error: Could not determine GRM column in %s\n", grm_fname);
        fclose(f);
        return NULL;
    }
    
    fprintf(stderr, "Detected GRM in column %d (0-based) of %s\n", grm_column, grm_fname);
    
    /* Count data lines */
    uint32_t n_lines = 0;
    if (!has_header) n_lines = 1; /* First line was data */
    
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
        
        /* Parse tokens to reach GRM column */
        char line_copy[MAX_LINE];
        strcpy(line_copy, line);
        
        char* tokens[10];  /* Assume max 10 columns */
        int n_tokens = 0;
        char* token = strtok(line_copy, " \t\n\r");
        while (token && n_tokens < 10) {
            tokens[n_tokens++] = token;
            token = strtok(NULL, " \t\n\r");
        }
        
        if (n_tokens < grm_column + 1) {
            fprintf(stderr, "Warning: Not enough columns in line %u, skipping\n", i + 1);
            continue;
        }
        
        /* Extract IID1, IID2, and GRM */
        pairs[i].iid1 = malloc(strlen(tokens[0]) + 1);
        pairs[i].iid2 = malloc(strlen(tokens[1]) + 1);
        if (!pairs[i].iid1 || !pairs[i].iid2) {
            for (uint32_t j = 0; j < i; j++) {
                free(pairs[j].iid1);
                free(pairs[j].iid2);
            }
            free(pairs);
            fclose(f);
            return NULL;
        }
        
        strcpy(pairs[i].iid1, tokens[0]);
        strcpy(pairs[i].iid2, tokens[1]);
        pairs[i].grm_value = strtod(tokens[grm_column], NULL);
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
                          const char* pheno_fname, uint32_t* n_phenos_out, 
                          char*** pheno_names_out)
{
    FILE* f = fopen(pheno_fname, "r");
    if (!f) { 
        perror(pheno_fname); 
        return 1; 
    }

    char line[MAX_LINE];
    uint32_t n_phenos = 0;
    char** pheno_names = NULL;
    
    /* Read header */
    if (fgets(line, sizeof(line), f)) {
        char* token = strtok(line, " \t\n\r");
        while (token) {
            n_phenos++;
            token = strtok(NULL, " \t\n\r");
        }
        
        if (n_phenos < 3) {
            fprintf(stderr, "Error: Need at least 3 columns in %s (IID1 IID2 pheno...)\n", pheno_fname);
            fclose(f);
            return 1;
        }
        
        n_phenos -= 2; /* Subtract IID columns */
        pheno_names = malloc(n_phenos * sizeof(char*));
        
        /* Re-read header to extract phenotype names */
        rewind(f);
        if (!fgets(line, sizeof(line), f)) {
            free(pheno_names);
            fclose(f);
            return 1;
        }
        
        token = strtok(line, " \t\n\r"); /* Skip IID1 */
        token = strtok(NULL, " \t\n\r"); /* Skip IID2 */
        
        for (uint32_t p = 0; p < n_phenos; p++) {
            token = strtok(NULL, " \t\n\r");
            if (token) {
                /* Remove newline if present */
                char* newline = strchr(token, '\n');
                if (newline) *newline = '\0';
                char* carriage = strchr(token, '\r');
                if (carriage) *carriage = '\0';
                
                pheno_names[p] = malloc(strlen(token) + 1);
                strcpy(pheno_names[p], token);
            } else {
                char default_name[32];
                snprintf(default_name, sizeof(default_name), "pheno%u", p + 1);
                pheno_names[p] = malloc(strlen(default_name) + 1);
                strcpy(pheno_names[p], default_name);
            }
        }
    }
    
    /* Allocate phenotype arrays for all pairs */
    for (uint32_t i = 0; i < n_pairs; i++) {
        pairs[i].pheno_values = malloc(n_phenos * sizeof(double));
        pairs[i].n_phenos = n_phenos;
        for (uint32_t p = 0; p < n_phenos; p++) {
            pairs[i].pheno_values[p] = NAN;
        }
    }
    
    /* Read phenotype data */
    uint32_t matched = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '\n') continue;
        
        char* token = strtok(line, " \t\n\r");
        if (!token) continue;
        char* iid1 = token;
        
        token = strtok(NULL, " \t\n\r");
        if (!token) continue;
        char* iid2 = token;
        
        /* Find matching pair - bidirectional search */
        for (uint32_t i = 0; i < n_pairs; i++) {
            if ((strcmp(pairs[i].iid1, iid1) == 0 && strcmp(pairs[i].iid2, iid2) == 0) ||
                (strcmp(pairs[i].iid1, iid2) == 0 && strcmp(pairs[i].iid2, iid1) == 0)) {
                
                /* Read phenotype values */
                for (uint32_t p = 0; p < n_phenos; p++) {
                    token = strtok(NULL, " \t\n\r");
                    if (token) {
                        char* endptr;
                        double val = strtod(token, &endptr);
                        if (*endptr == '\0' && isfinite(val)) {
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
    
    *n_phenos_out = n_phenos;
    *pheno_names_out = pheno_names;
    return 0;
}

/* Enhanced binned analysis with configurable parameters */
static int binned_analysis(const pair_data_t* pairs, uint32_t n_pairs, uint32_t n_phenos,
                          char** pheno_names, const char* out_fname,
                          double bin_width, double min_grm, double max_grm)
{
    /* Calculate number of bins */
    int n_bins = (int)((max_grm - min_grm) / bin_width) + 1;
    if (n_bins > MAX_BINS) {
        fprintf(stderr, "Error: Too many bins (%d), increase bin_width or reduce range\n", n_bins);
        return 1;
    }
    
    /* Initialize bins */
    bin_stats_t* bins = calloc(n_bins, sizeof(bin_stats_t));
    if (!bins) return 1;
    
    for (int b = 0; b < n_bins; b++) {
        bins[b].bin_center = min_grm + b * bin_width;
        bins[b].sum_pheno = calloc(n_phenos, sizeof(double));
        bins[b].sum_pheno_sq = calloc(n_phenos, sizeof(double));
        bins[b].count = calloc(n_phenos, sizeof(uint32_t));
        bins[b].n_phenos = n_phenos;
        bins[b].grm_sum = 0.0;
        bins[b].grm_sum_sq = 0.0;
        bins[b].grm_count = 0;
    }
    
    /* Accumulate statistics */
    uint32_t total_processed = 0;
    for (uint32_t i = 0; i < n_pairs; i++) {
        if (pairs[i].grm_value < min_grm || pairs[i].grm_value > max_grm) continue;
        
        int bin_idx = (int)((pairs[i].grm_value - min_grm) / bin_width);
        if (bin_idx >= n_bins) bin_idx = n_bins - 1;
        if (bin_idx < 0) bin_idx = 0;
        
        /* Add GRM statistics */
        bins[bin_idx].grm_sum += pairs[i].grm_value;
        bins[bin_idx].grm_sum_sq += pairs[i].grm_value * pairs[i].grm_value;
        bins[bin_idx].grm_count++;
        
        /* Add phenotype statistics */
        if (pairs[i].pheno_values) {
            for (uint32_t p = 0; p < n_phenos; p++) {
                if (isfinite(pairs[i].pheno_values[p])) {
                    bins[bin_idx].sum_pheno[p] += pairs[i].pheno_values[p];
                    bins[bin_idx].sum_pheno_sq[p] += pairs[i].pheno_values[p] * pairs[i].pheno_values[p];
                    bins[bin_idx].count[p]++;
                }
            }
        }
        total_processed++;
    }
    
    /* Write results */
    FILE* out = fopen(out_fname, "w");
    if (!out) {
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
    fprintf(out, "grm_bin\tgrm_mean\tgrm_n");
    for (uint32_t p = 0; p < n_phenos; p++) {
        fprintf(out, "\t%s_mean\t%s_n\t%s_se", pheno_names[p], pheno_names[p], pheno_names[p]);
    }
    fprintf(out, "\n");
    
    /* Write data for non-empty bins */
    for (int b = 0; b < n_bins; b++) {
        if (bins[b].grm_count > 0) {
            double grm_mean = bins[b].grm_sum / bins[b].grm_count;
            fprintf(out, "%.3f\t%.6f\t%u", bins[b].bin_center, grm_mean, bins[b].grm_count);
            
            for (uint32_t p = 0; p < n_phenos; p++) {
                if (bins[b].count[p] > 0) {
                    double mean = bins[b].sum_pheno[p] / bins[b].count[p];
                    double se = NAN;
                    if (bins[b].count[p] > 1) {
                        double var = (bins[b].sum_pheno_sq[p] - bins[b].sum_pheno[p] * bins[b].sum_pheno[p] / bins[b].count[p]) / (bins[b].count[p] - 1);
                        se = sqrt(var / bins[b].count[p]);
                    }
                    fprintf(out, "\t%.6f\t%u\t%.6f", mean, bins[b].count[p], se);
                } else {
                    fprintf(out, "\tNA\t0\tNA");
                }
            }
            fprintf(out, "\n");
        }
    }
    
    fclose(out);
    
    /* Cleanup */
    for (int b = 0; b < n_bins; b++) {
        free(bins[b].sum_pheno);
        free(bins[b].sum_pheno_sq);
        free(bins[b].count);
    }
    free(bins);
    
    fprintf(stderr, "Processed %u pairs into %d bins (%.3f to %.3f, width %.4f)\n", 
            total_processed, n_bins, min_grm, max_grm, bin_width);
    fprintf(stderr, "Binned analysis complete. Output written to %s\n", out_fname);
    
    return 0;
}

int main(int argc, char* argv[])
{
    if (argc < 4 || argc > 7) {
        fprintf(stderr, "Usage: %s <grm_file> <pheno_cross_file> <out_file> [bin_width] [min_grm] [max_grm]\n", argv[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "Arguments:\n");
        fprintf(stderr, "  grm_file        -- GRM results file (IID1 IID2 [N_SNPs] GRM)\n");
        fprintf(stderr, "  pheno_cross_file -- Phenotype cross-product file (IID1 IID2 pheno...)\n");
        fprintf(stderr, "  out_file        -- Output file for binned analysis\n");
        fprintf(stderr, "  bin_width       -- Bin width (default: %.4f)\n", DEFAULT_BIN_WIDTH);
        fprintf(stderr, "  min_grm         -- Minimum GRM value (default: %.3f)\n", DEFAULT_MIN_GRM);
        fprintf(stderr, "  max_grm         -- Maximum GRM value (default: %.3f)\n", DEFAULT_MAX_GRM);
        fprintf(stderr, "\n");
        fprintf(stderr, "Features:\n");
        fprintf(stderr, "  - Auto-detects GRM column (handles both 3 and 4 column formats)\n");
        fprintf(stderr, "  - Bidirectional pair matching (IID1,IID2 matches IID2,IID1)\n");
        fprintf(stderr, "  - Configurable binning parameters\n");
        fprintf(stderr, "  - Outputs mean GRM values per bin\n");
        return 1;
    }
    
    const char* grm_fname = argv[1];
    const char* pheno_fname = argv[2];
    const char* out_fname = argv[3];
    
    double bin_width = DEFAULT_BIN_WIDTH;
    double min_grm = DEFAULT_MIN_GRM;
    double max_grm = DEFAULT_MAX_GRM;
    
    if (argc > 4) bin_width = atof(argv[4]);
    if (argc > 5) min_grm = atof(argv[5]);
    if (argc > 6) max_grm = atof(argv[6]);
    
    if (bin_width <= 0 || min_grm >= max_grm) {
        fprintf(stderr, "Error: Invalid parameters (bin_width=%.4f, min_grm=%.3f, max_grm=%.3f)\n",
                bin_width, min_grm, max_grm);
        return 1;
    }
    
    fprintf(stderr, "Binning parameters: width=%.4f, range=[%.3f, %.3f]\n", bin_width, min_grm, max_grm);
    
    /* Read GRM file */
    uint32_t n_pairs;
    pair_data_t* pairs = read_grm_file(grm_fname, &n_pairs);
    if (!pairs) return 1;
    
    /* Join with phenotype data */
    uint32_t n_phenos;
    char** pheno_names;
    if (join_pheno_data(pairs, n_pairs, pheno_fname, &n_phenos, &pheno_names)) {
        for (uint32_t i = 0; i < n_pairs; i++) {
            free(pairs[i].iid1);
            free(pairs[i].iid2);
        }
        free(pairs);
        return 1;
    }
    
    /* Perform binned analysis */
    int result = binned_analysis(pairs, n_pairs, n_phenos, pheno_names, out_fname,
                                bin_width, min_grm, max_grm);
    
    /* Cleanup */
    for (uint32_t i = 0; i < n_pairs; i++) {
        free(pairs[i].iid1);
        free(pairs[i].iid2);
        if (pairs[i].pheno_values) free(pairs[i].pheno_values);
    }
    free(pairs);
    
    for (uint32_t p = 0; p < n_phenos; p++) {
        free(pheno_names[p]);
    }
    free(pheno_names);
    
    return result;
}
