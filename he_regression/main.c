#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_LINE_LENGTH 2048
#define MAX_INTERVALS 20
#define MAX_PAIRS 10000000
#define DEFAULT_BOOTSTRAP_REPS 1000

// Function to safely parse double with missing data handling
double safe_atof(const char* str) {
    if (!str || strlen(str) == 0) return NAN;
    
    // Check for common missing value representations
    if (strcmp(str, "NA") == 0 || strcmp(str, "nan") == 0 || 
        strcmp(str, "NaN") == 0 || strcmp(str, ".") == 0 ||
        strcmp(str, "-") == 0 || strcmp(str, "NULL") == 0) {
        return NAN;
    }
    
    char* endptr;
    double val = strtod(str, &endptr);
    
    // Check if entire string was consumed (valid number)
    if (endptr == str || *endptr != '\0') {
        return NAN;
    }
    
    return val;
}

typedef struct {
    double grm_value;
    double phenotype_product;
} Pair;

typedef struct {
    double min_grm;
    double max_grm;
    char name[64];
} Interval;

typedef struct {
    double intercept;
    double slope;
    double r_squared;
    double mse;
    int n_pairs;
    double intercept_se;
    double slope_se;
    double intercept_ci_lower;
    double intercept_ci_upper;
    double slope_ci_lower;
    double slope_ci_upper;
} RegressionResult;

// Function to parse intervals from command line
int parse_intervals(const char* interval_str, Interval* intervals) {
    int count = 0;
    char* str_copy = strdup(interval_str);
    char* token = strtok(str_copy, ",");
    
    while (token != NULL && count < MAX_INTERVALS) {
        char* dash_pos = strchr(token, '-');
        if (dash_pos != NULL) {
            *dash_pos = '\0';
            intervals[count].min_grm = safe_atof(token);
            intervals[count].max_grm = safe_atof(dash_pos + 1);
            
            // Skip intervals with missing values
            if (isnan(intervals[count].min_grm) || isnan(intervals[count].max_grm)) {
                printf("Warning: Skipping interval with missing values: %s\n", token);
                continue;
            }
            
            snprintf(intervals[count].name, sizeof(intervals[count].name), 
                    "%.3f-%.3f", intervals[count].min_grm, intervals[count].max_grm);
            count++;
        }
        token = strtok(NULL, ",");
    }
    
    free(str_copy);
    return count;
}

// Function to filter pairs by GRM interval
int filter_pairs_by_interval(Pair* all_pairs, int total_pairs, 
                             const Interval* interval, Pair* filtered_pairs) {
    int count = 0;
    for (int i = 0; i < total_pairs; i++) {
        if (all_pairs[i].grm_value >= interval->min_grm && 
            all_pairs[i].grm_value < interval->max_grm) {
            filtered_pairs[count++] = all_pairs[i];
        }
    }
    return count;
}

// Simple linear regression
RegressionResult perform_regression(Pair* pairs, int n_pairs) {
    RegressionResult result = {0};
    
    if (n_pairs < 2) {
        result.n_pairs = n_pairs;
        return result;
    }
    
    // Calculate means
    double sum_x = 0, sum_y = 0;
    for (int i = 0; i < n_pairs; i++) {
        sum_x += pairs[i].grm_value;
        sum_y += pairs[i].phenotype_product;
    }
    double mean_x = sum_x / n_pairs;
    double mean_y = sum_y / n_pairs;
    
    // Calculate slope and intercept
    double numerator = 0, denominator = 0;
    for (int i = 0; i < n_pairs; i++) {
        double dx = pairs[i].grm_value - mean_x;
        double dy = pairs[i].phenotype_product - mean_y;
        numerator += dx * dy;
        denominator += dx * dx;
    }
    
    if (denominator == 0) {
        result.n_pairs = n_pairs;
        return result;
    }
    
    result.slope = numerator / denominator;
    result.intercept = mean_y - result.slope * mean_x;
    result.n_pairs = n_pairs;
    
    // Calculate R-squared and MSE
    double ss_tot = 0, ss_res = 0;
    for (int i = 0; i < n_pairs; i++) {
        double predicted = result.intercept + result.slope * pairs[i].grm_value;
        double residual = pairs[i].phenotype_product - predicted;
        ss_res += residual * residual;
        ss_tot += (pairs[i].phenotype_product - mean_y) * (pairs[i].phenotype_product - mean_y);
    }
    
    result.r_squared = (ss_tot > 0) ? 1.0 - (ss_res / ss_tot) : 0.0;
    result.mse = ss_res / (n_pairs - 2);
    
    // Calculate standard errors
    if (n_pairs > 2 && denominator > 0) {
        double s_xx = denominator;
        double se_slope_sq = result.mse / s_xx;
        result.slope_se = sqrt(se_slope_sq);
        
        double se_intercept_sq = result.mse * (1.0/n_pairs + (mean_x * mean_x) / s_xx);
        result.intercept_se = sqrt(se_intercept_sq);
    }
    
    return result;
}

// Bootstrap sampling
void bootstrap_sample(Pair* original_pairs, int n_pairs, Pair* boot_pairs) {
    for (int i = 0; i < n_pairs; i++) {
        int idx = rand() % n_pairs;
        boot_pairs[i] = original_pairs[idx];
    }
}

// Comparison function for qsort
int compare_double(const void* a, const void* b) {
    double da = *(double*)a;
    double db = *(double*)b;
    return (da > db) - (da < db);
}

// Perform bootstrapped regression
void bootstrap_regression(Pair* pairs, int n_pairs, int n_bootstrap, 
                         RegressionResult* result) {
    if (n_pairs < 2) {
        return;
    }
    
    double* boot_intercepts = malloc(n_bootstrap * sizeof(double));
    double* boot_slopes = malloc(n_bootstrap * sizeof(double));
    Pair* boot_pairs = malloc(n_pairs * sizeof(Pair));
    
    // Perform bootstrap replicates
    for (int b = 0; b < n_bootstrap; b++) {
        bootstrap_sample(pairs, n_pairs, boot_pairs);
        RegressionResult boot_result = perform_regression(boot_pairs, n_pairs);
        boot_intercepts[b] = boot_result.intercept;
        boot_slopes[b] = boot_result.slope;
    }
    
    // Sort for confidence intervals
    
    qsort(boot_intercepts, n_bootstrap, sizeof(double), compare_double);
    qsort(boot_slopes, n_bootstrap, sizeof(double), compare_double);
    
    // Calculate 95% confidence intervals (2.5% and 97.5% percentiles)
    int lower_idx = (int)(0.025 * n_bootstrap);
    int upper_idx = (int)(0.975 * n_bootstrap);
    
    result->intercept_ci_lower = boot_intercepts[lower_idx];
    result->intercept_ci_upper = boot_intercepts[upper_idx];
    result->slope_ci_lower = boot_slopes[lower_idx];
    result->slope_ci_upper = boot_slopes[upper_idx];
    
    free(boot_intercepts);
    free(boot_slopes);
    free(boot_pairs);
}

// Auto-detect file format
int detect_format(const char* filename) {
    FILE* fp = fopen(filename, "r");
    if (!fp) return 0;
    
    char line[MAX_LINE_LENGTH];
    if (fgets(line, sizeof(line), fp)) {
        int field_count = 0;
        char* token = strtok(line, " \t\n");
        while (token != NULL) {
            field_count++;
            token = strtok(NULL, " \t\n");
        }
        fclose(fp);
        return field_count;
    }
    fclose(fp);
    return 0;
}

// Read pairs from joined file
int read_pairs(const char* filename, Pair* pairs) {
    FILE* fp = fopen(filename, "r");
    if (!fp) {
        printf("Error: Cannot open file %s\n", filename);
        return -1;
    }
    
    char line[MAX_LINE_LENGTH];
    int count = 0;
    int grm_col = -1, pheno_col = -1;
    
    // Skip header and detect format
    if (fgets(line, sizeof(line), fp)) {
        // Auto-detect GRM column (contains "grm" or "relatedness") and phenotype column
        char* token = strtok(line, " \t\n");
        int col = 0;
        while (token != NULL && col < 20) {
            if (strstr(token, "grm") != NULL || strstr(token, "relatedness") != NULL || 
                strstr(token, "kinship") != NULL) {
                grm_col = col;
            }
            if (strstr(token, "pheno") != NULL || strstr(token, "product") != NULL ||
                strstr(token, "cross") != NULL) {
                pheno_col = col;
            }
            token = strtok(NULL, " \t\n");
            col++;
        }
        
        // Default columns if not auto-detected
        if (grm_col == -1) grm_col = 2;  // Assume 3rd column
        if (pheno_col == -1) pheno_col = 3; // Assume 4th column
        
        printf("Using GRM column %d and phenotype column %d\n", grm_col + 1, pheno_col + 1);
    }
    
    // Read data
    while (fgets(line, sizeof(line), fp) && count < MAX_PAIRS) {
        char* token = strtok(line, " \t\n");
        int col = 0;
        double grm_val = 0, pheno_val = 0;
        int grm_found = 0, pheno_found = 0;
        
        while (token != NULL) {
            if (col == grm_col) {
                grm_val = safe_atof(token);
                grm_found = 1;
            }
            if (col == pheno_col) {
                pheno_val = safe_atof(token);
                pheno_found = 1;
            }
            token = strtok(NULL, " \t\n");
            col++;
        }
        
        // Only include pairs with valid (non-missing) data
        if (grm_found && pheno_found && !isnan(grm_val) && !isnan(pheno_val)) {
            pairs[count].grm_value = grm_val;
            pairs[count].phenotype_product = pheno_val;
            count++;
        }
    }
    
    fclose(fp);
    
    int total_lines = count + 1; // +1 for header
    int valid_pairs = count;
    printf("Read %d pairs from %s\n", valid_pairs, filename);
    if (total_lines > valid_pairs + 1) {
        printf("  Note: %d lines with missing data were excluded\n", total_lines - valid_pairs - 1);
    }
    return count;
}

void print_usage(const char* program_name) {
    printf("Usage: %s -f <input_file> -i <intervals> [options]\n", program_name);
    printf("Options:\n");
    printf("  -f <file>       Input file with GRM values and phenotype cross-products\n");
    printf("  -i <intervals>  Comma-separated intervals (e.g., \"0-0.05,0.05-0.35,0.35-0.75\")\n");
    printf("  -b <number>     Bootstrap replicates (default: 1000)\n");
    printf("  -o <file>       Output file (default: he_regression_results.txt)\n");
    printf("  -s <seed>       Random seed (default: time-based)\n");
    printf("  -h              Show this help\n");
    printf("\nExample:\n");
    printf("  %s -f joined_data.txt -i \"0-0.05,0.05-0.35,0.35-0.75\" -b 1000\n", program_name);
}

int main(int argc, char* argv[]) {
    char* input_file = NULL;
    char* interval_str = NULL;
    char* output_file = "he_regression_results.txt";
    int n_bootstrap = DEFAULT_BOOTSTRAP_REPS;
    unsigned int seed = 0;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
            input_file = argv[++i];
        } else if (strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
            interval_str = argv[++i];
        } else if (strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            n_bootstrap = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            output_file = argv[++i];
        } else if (strcmp(argv[i], "-s") == 0 && i + 1 < argc) {
            seed = (unsigned int)atoi(argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0) {
            print_usage(argv[0]);
            return 0;
        }
    }
    
    if (!input_file || !interval_str) {
        printf("Error: Input file and intervals are required\n\n");
        print_usage(argv[0]);
        return 1;
    }
    
    // Set random seed
    if (seed == 0) {
        seed = (unsigned int)time(NULL);
    }
    srand(seed);
    printf("Using random seed: %u\n", seed);
    
    // Parse intervals
    Interval intervals[MAX_INTERVALS];
    int n_intervals = parse_intervals(interval_str, intervals);
    if (n_intervals == 0) {
        printf("Error: No valid intervals parsed\n");
        return 1;
    }
    
    printf("Parsed %d intervals:\n", n_intervals);
    for (int i = 0; i < n_intervals; i++) {
        printf("  %s: [%.3f, %.3f)\n", intervals[i].name, intervals[i].min_grm, intervals[i].max_grm);
    }
    
    // Read pairs
    Pair* all_pairs = malloc(MAX_PAIRS * sizeof(Pair));
    if (!all_pairs) {
        printf("Error: Cannot allocate memory for pairs\n");
        return 1;
    }
    
    int total_pairs = read_pairs(input_file, all_pairs);
    if (total_pairs <= 0) {
        printf("Error: No pairs read from file\n");
        free(all_pairs);
        return 1;
    }
    
    // Open output file
    FILE* out_fp = fopen(output_file, "w");
    if (!out_fp) {
        printf("Error: Cannot open output file %s\n", output_file);
        free(all_pairs);
        return 1;
    }
    
    // Write header
    fprintf(out_fp, "interval\tn_pairs\tintercept\tslope\tr_squared\tmse\t");
    fprintf(out_fp, "intercept_se\tslope_se\t");
    fprintf(out_fp, "intercept_ci_lower\tintercept_ci_upper\t");
    fprintf(out_fp, "slope_ci_lower\tslope_ci_upper\n");
    
    // Process each interval
    Pair* filtered_pairs = malloc(MAX_PAIRS * sizeof(Pair));
    if (!filtered_pairs) {
        printf("Error: Cannot allocate memory for filtered pairs\n");
        fclose(out_fp);
        free(all_pairs);
        return 1;
    }
    
    for (int i = 0; i < n_intervals; i++) {
        printf("\nProcessing interval %s...\n", intervals[i].name);
        
        // Filter pairs for this interval
        int n_filtered = filter_pairs_by_interval(all_pairs, total_pairs, 
                                                  &intervals[i], filtered_pairs);
        
        if (n_filtered < 2) {
            printf("  Warning: Only %d pairs in interval, skipping\n", n_filtered);
            fprintf(out_fp, "%s\t%d\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n", 
                    intervals[i].name, n_filtered);
            continue;
        }
        
        printf("  Found %d pairs in interval\n", n_filtered);
        
        // Perform basic regression
        RegressionResult result = perform_regression(filtered_pairs, n_filtered);
        
        // Perform bootstrap regression for confidence intervals
        printf("  Performing %d bootstrap replicates...\n", n_bootstrap);
        bootstrap_regression(filtered_pairs, n_filtered, n_bootstrap, &result);
        
        // Write results
        fprintf(out_fp, "%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
                intervals[i].name, result.n_pairs,
                result.intercept, result.slope, result.r_squared, result.mse,
                result.intercept_se, result.slope_se,
                result.intercept_ci_lower, result.intercept_ci_upper,
                result.slope_ci_lower, result.slope_ci_upper);
        
        printf("  Results: intercept=%.4f, slope=%.4f, R²=%.4f, n=%d\n",
               result.intercept, result.slope, result.r_squared, result.n_pairs);
        printf("  Bootstrap 95%% CI: slope [%.4f, %.4f]\n",
               result.slope_ci_lower, result.slope_ci_upper);
    }
    
    printf("\nResults written to %s\n", output_file);
    
    free(filtered_pairs);
    free(all_pairs);
    fclose(out_fp);
    
    return 0;
}
