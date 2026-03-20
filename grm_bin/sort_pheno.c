#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char *id;
    char *iid_only;
    float pheno;
    int has_pheno;
    int has_valid_iid;  // Track if IID is >= 0
} Individual;

int compare_strings(const void *a, const void *b) {
    return strcmp(*(const char **)a, *(const char **)b);
}

int compare_individuals_by_id(const void *a, const void *b) {
    char *id_a = ((Individual *)a)->id;
    char *id_b = ((Individual *)b)->id;
    return strcmp(id_a, id_b);
}

void print_usage(char *program_name) {
    fprintf(stderr, "Usage: %s fileNameGrmPrefix fileNamePhenoIn fileNamePhenoOut\n", program_name);
    fprintf(stderr, "Sorts phenotype file to match GRM individual order\n");
    fprintf(stderr, "Missing values in input are coded as NA, empty, or any non-numeric value\n");
    fprintf(stderr, "Output missing values are coded as -999\n");
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }
    
    char *grm_prefix = argv[1];
    char *pheno_in = argv[2];
    char *pheno_out = argv[3];
    
    // Read GRM ID file to get individual order
    char grm_id_file[500];
    sprintf(grm_id_file, "%s.grm.id", grm_prefix);
    
    FILE *f = fopen(grm_id_file, "r");
    if (!f) {
        fprintf(stderr, "Error: cannot open file %s\n", grm_id_file);
        return 1;
    }
    
    // Count individuals
    long long n_indiv = 0;
    char line[500];
    while (fgets(line, sizeof(line), f)) {
        n_indiv++;
    }
    
    printf("Found %lld individuals in GRM ID file\n", n_indiv);
    
    // Read GRM IDs in order
    Individual *grm_order = (Individual *)malloc(n_indiv * sizeof(Individual));
    rewind(f);
    
    for (long long i = 0; i < n_indiv; i++) {
        if (fgets(line, sizeof(line), f) == NULL) {
            fprintf(stderr, "Error reading GRM ID file at line %lld\n", i+1);
            return 1;
        }
        
        // Parse line: typically "FID IID" format
        char fid[256], iid_str[256];
        if (sscanf(line, "%s %s", fid, iid_str) != 2) {
            fprintf(stderr, "Error parsing GRM ID file at line %lld\n", i+1);
            return 1;
        }
        
        // Parse IID as integer to check if it's >= 0
        long long iid_val = atoll(iid_str);
        
        // Store individual ID (use IID, or FID_IID if needed)
        grm_order[i].id = (char *)malloc(512 * sizeof(char));
        sprintf(grm_order[i].id, "%s_%s", fid, iid_str);
        
        // Also store just the IID for matching phenotype files that use only IID
        grm_order[i].iid_only = (char *)malloc(256 * sizeof(char));
        strcpy(grm_order[i].iid_only, iid_str);
        
        // Mark as missing if IID < 0
        if (iid_val < 0) {
            grm_order[i].has_pheno = 0;
            grm_order[i].pheno = -999.0f;
            grm_order[i].has_valid_iid = 0;
        } else {
            grm_order[i].has_pheno = 0;
            grm_order[i].pheno = -999.0f;
            grm_order[i].has_valid_iid = 1;
        }
    }
    
    fclose(f);
    
    // Read phenotype input file
    f = fopen(pheno_in, "r");
    if (!f) {
        fprintf(stderr, "Error: cannot open file %s\n", pheno_in);
        return 1;
    }
    
    long long n_matched = 0;
    long long n_missing = 0;
    
    // Skip header line if it exists
    char header_line[500];
    if (fgets(header_line, sizeof(header_line), f) != NULL) {
        // Check if it looks like a header (contains non-numeric characters)
        char test_id[256], test_pheno[256];
        if (sscanf(header_line, "%s %s", test_id, test_pheno) == 2) {
            char *endptr;
            strtol(test_id, &endptr, 10);
            if (*endptr != '\0') {
                // First line looks like header, skip it
                printf("Skipping header line in phenotype file\n");
            } else {
                // First line looks like data, rewind to process it
                rewind(f);
            }
        }
    }
    
    while (fgets(line, sizeof(line), f)) {
        char id[256];
        char pheno_str[256];
        
        // Parse line
        if (sscanf(line, "%s %s", id, pheno_str) != 2) {
            continue;
        }
        
        // Check if pheno value is numeric
        char *endptr;
        float pheno_val = strtof(pheno_str, &endptr);
        
        // Search for this individual in GRM order
        // Try matching both FID_IID format and IID-only format
        int matched = 0;
        for (long long i = 0; i < n_indiv; i++) {
            if (strcmp(grm_order[i].id, id) == 0 || strcmp(grm_order[i].iid_only, id) == 0) {
                // Only assign phenotype if IID is valid (>= 0)
                if (grm_order[i].has_valid_iid) {
                    // Check if value is actually numeric
                    if (*endptr == '\0' || *endptr == '\n') {
                        grm_order[i].pheno = pheno_val;
                        grm_order[i].has_pheno = 1;
                        n_matched++;
                    } else {
                        // Not numeric, treat as missing
                        grm_order[i].pheno = -999.0f;
                        grm_order[i].has_pheno = 0;
                        n_missing++;
                    }
                } else {
                    // IID < 0, mark as missing regardless
                    grm_order[i].pheno = -999.0f;
                    grm_order[i].has_pheno = 0;
                }
                matched = 1;
                break;
            }
        }
        
        if (!matched) {
            // Individual not found in GRM - this shouldn't happen but log it
            printf("Warning: Individual %s not found in GRM ID file\n", id);
        }
    }
    
    fclose(f);
    
    printf("Matched %lld individuals with phenotype data\n", n_matched);
    printf("Found %lld missing values in phenotype\n", n_missing);
    
    // Write output in GRM order
    FILE *out = fopen(pheno_out, "w");
    if (!out) {
        fprintf(stderr, "Error: cannot create output file %s\n", pheno_out);
        return 1;
    }
    
    long long n_unmatched = 0;
    long long n_negative_iid = 0;
    for (long long i = 0; i < n_indiv; i++) {
        if (grm_order[i].has_pheno) {
            fprintf(out, "%s %.4f\n", grm_order[i].id, grm_order[i].pheno);
        } else {
            fprintf(out, "%s -999\n", grm_order[i].id);
            if (!grm_order[i].has_valid_iid) {
                n_negative_iid++;
            } else {
                n_unmatched++;
            }
        }
    }
    
    fclose(out);
    
    printf("Output written to %s\n", pheno_out);
    printf("Individuals without phenotype data: %lld\n", n_unmatched);
    printf("Individuals with negative IID: %lld\n", n_negative_iid);
    
    // Free memory
    for (long long i = 0; i < n_indiv; i++) {
        free(grm_order[i].id);
        free(grm_order[i].iid_only);
    }
    free(grm_order);
    
    return 0;
}
