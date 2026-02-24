/*
 * main.c -- Calculate cross products of phenotypic values for pairs of individuals
 *
 * Usage: pheno_pairs <pair_file> <pheno_file> <out_file>
 *
 * Expects:
 *   <pair_file> -- tab-delimited file with pairs (IID1 IID2 ...)
 *   <pheno_file> -- tab-delimited file with phenotypes (IID pheno)
 *   <out_file>  -- output file for cross products
 *
 * For each pair (IID1, IID2), calculates pheno1 * pheno2 and outputs:
 *   IID1 IID2 pheno1 pheno2 cross_product
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* Structure to store phenotype data */
typedef struct {
    char* iid;
    double* phenos;  /* Array of phenotype values */
    uint32_t n_phenos;
} pheno_entry_t;

/* ------------------------------------------------------------------ */
/* Read phenotype file -- returns array of pheno_entry_t, sets count  */
/* Also reads column headers and sets pheno_names and n_phenos        */
/* ------------------------------------------------------------------ */
static pheno_entry_t* read_pheno_file(const char* pheno_fname, 
                                       uint32_t* count_out,
                                       char*** pheno_names_out,
                                       uint32_t* n_phenos_out)
{
    FILE* f = fopen(pheno_fname, "r");
    if (!f) { perror(pheno_fname); return NULL; }

    char line[1024];
    uint32_t n_phenos = 0;
    char** pheno_names = NULL;
    
    /* Read first line to get column headers */
    if (fgets(line, sizeof(line), f)) {
        /* Count columns */
        char* token = strtok(line, " \t\n\r");
        while (token) {
            n_phenos++;
            token = strtok(NULL, " \t\n\r");
        }
        
        if (n_phenos < 2) {
            fprintf(stderr, "Error: Need at least 2 columns (IID + phenotypes)\n");
            fclose(f);
            return NULL;
        }
        
        /* Allocate space for phenotype names (skip first column which is IID) */
        n_phenos--; /* Subtract 1 for IID column */
        pheno_names = malloc(n_phenos * sizeof(char*));
        if (!pheno_names) { fclose(f); return NULL; }
        
        /* Re-read first line to extract column names */
        rewind(f);
        if (!fgets(line, sizeof(line), f)) {
            fprintf(stderr, "Error: Unable to read header line from %s\n", pheno_fname);
            fclose(f);
            free(pheno_names);
            return NULL;
        }
        
        token = strtok(line, " \t\n\r"); /* Skip IID column */
        for (uint32_t i = 0; i < n_phenos; i++) {
            token = strtok(NULL, " \t\n\r");
            if (token) {
                pheno_names[i] = malloc(strlen(token) + 1);
                if (!pheno_names[i]) {
                    for (uint32_t j = 0; j < i; j++) free(pheno_names[j]);
                    free(pheno_names); fclose(f); return NULL;
                }
                strcpy(pheno_names[i], token);
            } else {
                /* Fallback name if header parsing fails */
                char default_name[32];
                snprintf(default_name, sizeof(default_name), "pheno%u", i+1);
                pheno_names[i] = malloc(strlen(default_name) + 1);
                strcpy(pheno_names[i], default_name);
            }
        }
    } else {
        fprintf(stderr, "Error: Empty file %s\n", pheno_fname);
        fclose(f);
        return NULL;
    }

    /* count data lines */
    uint32_t n = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '\0' && line[0] != '\n') n++;
    }
    rewind(f);
    if (!fgets(line, sizeof(line), f)) {
        fprintf(stderr, "Error: Unable to skip header in %s\n", pheno_fname);
        fclose(f);
        for (uint32_t i = 0; i < n_phenos; i++) free(pheno_names[i]);
        free(pheno_names);
        return NULL;
    } /* Skip header */

    pheno_entry_t* entries = malloc(n * sizeof(pheno_entry_t));
    if (!entries) { 
        for (uint32_t i = 0; i < n_phenos; i++) free(pheno_names[i]);
        free(pheno_names); 
        fclose(f); 
        return NULL; 
    }

    uint32_t i = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '\n') continue;
        
        /* Parse IID */
        char* token = strtok(line, " \t\n\r");
        if (!token) continue;
        
        entries[i].iid = malloc(strlen(token) + 1);
        if (!entries[i].iid) {
            for (uint32_t j = 0; j < i; j++) {
                free(entries[j].iid);
                free(entries[j].phenos);
            }
            free(entries);
            for (uint32_t j = 0; j < n_phenos; j++) free(pheno_names[j]);
            free(pheno_names);
            fclose(f);
            return NULL;
        }
        strcpy(entries[i].iid, token);
        
        /* Allocate and parse phenotypes */
        entries[i].phenos = malloc(n_phenos * sizeof(double));
        entries[i].n_phenos = n_phenos;
        if (!entries[i].phenos) {
            free(entries[i].iid);
            for (uint32_t j = 0; j < i; j++) {
                free(entries[j].iid);
                free(entries[j].phenos);
            }
            free(entries);
            for (uint32_t j = 0; j < n_phenos; j++) free(pheno_names[j]);
            free(pheno_names);
            fclose(f);
            return NULL;
        }
        
        /* Parse all phenotype values */
        uint32_t parsed_phenos = 0;
        for (uint32_t p = 0; p < n_phenos; p++) {
            token = strtok(NULL, " \t\n\r");
            if (token) {
                char* endptr;
                double value = strtod(token, &endptr);
                if (*endptr == '\0') {
                    entries[i].phenos[p] = value;
                    parsed_phenos++;
                } else {
                    entries[i].phenos[p] = 0.0/0.0; /* NaN for non-numeric */
                }
            } else {
                entries[i].phenos[p] = 0.0/0.0; /* NaN for missing */
            }
        }
        
        /* Only include entry if at least one phenotype was successfully parsed */
        if (parsed_phenos > 0) {
            i++;
        } else {
            free(entries[i].iid);
            free(entries[i].phenos);
        }
    }
    fclose(f);

    *count_out = i;
    *pheno_names_out = pheno_names;
    *n_phenos_out = n_phenos;
    fprintf(stderr, "Read %u phenotype entries with %u phenotype columns from %s\n", 
            i, n_phenos, pheno_fname);
    return entries;
}

/* ------------------------------------------------------------------ */
/* Find phenotypes for a given IID -- returns pointer to pheno array  */
/* ------------------------------------------------------------------ */
static const double* find_phenos(const pheno_entry_t* entries, uint32_t count, const char* iid)
{
    for (uint32_t i = 0; i < count; i++) {
        if (strcmp(entries[i].iid, iid) == 0) {
            return entries[i].phenos;
        }
    }
    return NULL; /* Not found */
}

/* ------------------------------------------------------------------ */
/* Process pairs file and calculate cross products                    */
/* ------------------------------------------------------------------ */
static int process_pairs(const char* pair_fname,
                        const pheno_entry_t* pheno_entries,
                        uint32_t pheno_count,
                        char** pheno_names,
                        uint32_t n_phenos,
                        const char* out_fname)
{
    FILE* pair_f = fopen(pair_fname, "r");
    if (!pair_f) { perror(pair_fname); return 1; }

    FILE* out_f = fopen(out_fname, "w");
    if (!out_f) { perror(out_fname); fclose(pair_f); return 1; }

    /* Write header to output file */
    fprintf(out_f, "IID1\tIID2");
    for (uint32_t p = 0; p < n_phenos; p++) {
        fprintf(out_f, "\t%s.crossprod", pheno_names[p]);
    }
    fprintf(out_f, "\n");

    char line[1024];
    uint32_t processed = 0, found_both = 0, missing = 0;
    
    while (fgets(line, sizeof(line), pair_f)) {
        if (line[0] == '\0' || line[0] == '\n') continue;
        
        char iid1[256], iid2[256];
        /* Parse first two columns as IID1 and IID2 */
        if (sscanf(line, "%255s %255s", iid1, iid2) == 2) {
            processed++;
            
            const double* phenos1 = find_phenos(pheno_entries, pheno_count, iid1);
            const double* phenos2 = find_phenos(pheno_entries, pheno_count, iid2);
            
            fprintf(out_f, "%s\t%s", iid1, iid2);
            
            if (phenos1 && phenos2) {
                /* Calculate cross products for each phenotype */
                int any_valid = 0;
                for (uint32_t p = 0; p < n_phenos; p++) {
                    double p1 = phenos1[p];
                    double p2 = phenos2[p];
                    
                    /* Check if both phenotypes are valid (not NaN) */
                    if (p1 == p1 && p2 == p2) { /* NaN != NaN */
                        double cross_product = p1 * p2;
                        fprintf(out_f, "\t%g", cross_product);
                        any_valid = 1;
                    } else {
                        fprintf(out_f, "\tNA");
                    }
                }
                if (any_valid) found_both++;
                else missing++;
            } else {
                /* One or both individuals not found */
                for (uint32_t p = 0; p < n_phenos; p++) {
                    fprintf(out_f, "\tNA");
                }
                missing++;
            }
            fprintf(out_f, "\n");
        }
    }

    fclose(pair_f);
    fclose(out_f);

    fprintf(stderr, "Processed %u pairs:\n", processed);
    fprintf(stderr, "  %u pairs with at least one valid phenotype\n", found_both);
    fprintf(stderr, "  %u pairs with no valid phenotypes\n", missing);

    return 0;
}

/* ------------------------------------------------------------------ */
/* main                                                               */
/* ------------------------------------------------------------------ */
int main(int argc, char** argv)
{
    if (argc < 4) {
        fprintf(stderr, "usage: pheno_pairs <pair_file> <pheno_file> <out_file>\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Calculates cross products of phenotypic values for pairs.\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Input files:\n");
        fprintf(stderr, "  <pair_file>  -- pairs file (IID1 IID2 ...)\n");
        fprintf(stderr, "  <pheno_file> -- phenotype file (IID pheno1 pheno2 ...)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "Output:\n");
        fprintf(stderr, "  <out_file>   -- IID1 IID2 <pheno_names>.crossprod\n");
        return 1;
    }

    const char* pair_file  = argv[1];
    const char* pheno_file = argv[2];
    const char* out_file   = argv[3];

    /* Read phenotype data */
    uint32_t pheno_count = 0;
    char** pheno_names = NULL;
    uint32_t n_phenos = 0;
    pheno_entry_t* pheno_entries = read_pheno_file(pheno_file, &pheno_count, &pheno_names, &n_phenos);
    if (!pheno_entries || pheno_count == 0) {
        fprintf(stderr, "Error reading phenotype file %s\n", pheno_file);
        return 1;
    }

    /* Process pairs and calculate cross products */
    int ret = process_pairs(pair_file, pheno_entries, pheno_count, pheno_names, n_phenos, out_file);

    /* cleanup */
    for (uint32_t i = 0; i < pheno_count; i++) {
        free(pheno_entries[i].iid);
        free(pheno_entries[i].phenos);
    }
    free(pheno_entries);
    
    for (uint32_t i = 0; i < n_phenos; i++) {
        free(pheno_names[i]);
    }
    free(pheno_names);

    if (ret == 0) {
        fprintf(stderr, "Done. Output written to %s\n", out_file);
    }

    return ret;
}
