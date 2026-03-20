#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    long long n_indiv;
    int n_bins;
    float *bins;
    float *y;
    int *na;
    char **ids;
    long long *covariance_n;
    double *covariance;
    double *avg_bins;
} Pheno_Data;

void print_usage(char *program_name) {
    fprintf(stderr, "Usage: %s nIndiv fileNamePhenoPrefix fileNameGrmPrefix\n", program_name);
    fprintf(stderr, "make sure pheno file and GRM are sorted in same order, this program does not check!\n");
    fprintf(stderr, "pheno file needs .txt suffix & grm needs .grm.id and .grm.bin suffix\n");
    fprintf(stderr, "missing values are coded -999!\n");
}

void setup_bins(Pheno_Data *data) {
    float bin_width = 0.001f;
    float bin_width_high = 0.005f;
    float min_val = -0.3f;
    float max_val = 1.05f;
    float threshold = 0.02f;
    
    // Allocate extra space for bins
    data->bins = (float *)malloc(500 * sizeof(float));
    
    int i = 0;
    data->bins[0] = min_val;
    data->bins[1] = -0.02f;
    i = 1;
    
    // Add uniform bins from -0.02 to 0.02 with width 0.001
    while (data->bins[i] <= threshold) {
        i++;
        data->bins[i] = data->bins[i-1] + bin_width;
    }
    
    // Add bins above 0.02 with width 0.005
    float current = threshold;
    while (current < max_val) {
        i++;
        current += bin_width_high;
        data->bins[i] = current;
    }
    
    // Ensure max_val is included as the upper bound
    if (data->bins[i] < max_val) {
        i++;
        data->bins[i] = max_val;
    } else {
        data->bins[i] = max_val;
    }
    
    data->n_bins = i;
    
    printf("genomic relationship: there are %d bins starting at %.4g until %.4g\n", 
           data->n_bins, data->bins[0], data->bins[data->n_bins]);
}

int read_phenotype(Pheno_Data *data, char *pheno_prefix) {
    char filename[500];
    sprintf(filename, "%s.txt", pheno_prefix);
    
    FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Error: cannot open file %s\n", filename);
        return 1;
    }
    
    float tmp;
    char id[500];
    
    for (long long i = 0; i < data->n_indiv; i++) {
        if (fscanf(f, "%s %f\n", id, &tmp) != 2) {
            fprintf(stderr, "Error reading phenotype file at line %lld\n", i+1);
            fclose(f);
            return 1;
        }
        
        if (tmp == -999.0f) {
            data->na[i] = 0;
        } else {
            data->y[i] = tmp;
            data->na[i] = 1;
        }
    }
    
    fclose(f);
    
    // Calculate mean and variance
    double mu = 0.0;
    long long valid_count = 0;
    
    for (long long i = 0; i < data->n_indiv; i++) {
        if (data->na[i]) {
            mu += data->y[i];
            valid_count++;
        }
    }
    mu /= valid_count;
    
    double var = 0.0;
    for (long long i = 0; i < data->n_indiv; i++) {
        if (data->na[i]) {
            var += (data->y[i] - mu) * (data->y[i] - mu);
        }
    }
    var /= valid_count;
    
    // Standardize
    double sd = sqrt(var);
    for (long long i = 0; i < data->n_indiv; i++) {
        data->y[i] = (data->y[i] - mu) / sd;
    }
    
    return 0;
}

long long get_cell_index(long long i, long long j) {
    // Calculate cell index from lower triangle
    // Equivalent to: i2 = i*(i-1)/2, then cell = i2 + j
    long long i2;
    if (i % 2 == 0) {
        i2 = (i / 2) * (i - 1);
    } else {
        i2 = ((i - 1) / 2 + 1) * ((i - 1) / 2);
    }
    return i2 + j;
}

int read_grm_and_calculate(Pheno_Data *data, char *grm_prefix) {
    char grm_bin[500], grm_id[500];
    sprintf(grm_bin, "%s.grm.bin", grm_prefix);
    sprintf(grm_id, "%s.grm.id", grm_prefix);
    
    FILE *f = fopen(grm_bin, "rb");
    if (!f) {
        fprintf(stderr, "Error: cannot open file %s\n", grm_bin);
        return 1;
    }
    
    // Initialize statistics
    for (int k = 0; k < data->n_bins; k++) {
        data->covariance_n[k] = 0;
        data->covariance[k] = 0.0;
        data->avg_bins[k] = 0.0;
    }
    
    // Calculate covariances
    for (long long i = 0; i < data->n_indiv; i++) {
        if (!data->na[i]) continue;
        
        for (long long j = 0; j < i; j++) {
            if (!data->na[j]) continue;
            
            long long cell = get_cell_index(i+1, j+1);
            long long m = 4 * (cell - 1);
            
            float grm_val;
            if (fseek(f, m, SEEK_SET) != 0) {
                fprintf(stderr, "Error seeking in GRM file\n");
                fclose(f);
                return 1;
            }
            
            if (fread(&grm_val, sizeof(float), 1, f) != 1) {
                fprintf(stderr, "Error reading GRM value at position %lld\n", m);
                fclose(f);
                return 1;
            }
            
            // Find bin
            int flag = 0;
            int k = (grm_val > 0) ? 20 : 1;
            
            while (!flag) {
                if (k >= data->n_bins) {
                    fprintf(stderr, "Error: GRM value %f outside bin range\n", grm_val);
                    break;
                }
                
                if (grm_val < data->bins[k+1]) {
                    flag = 1;
                    data->covariance_n[k]++;
                    data->covariance[k] += data->y[i] * data->y[j];
                    data->avg_bins[k] += grm_val;
                }
                k++;
            }
        }
    }
    
    fclose(f);
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }
    
    Pheno_Data data;
    
    // Parse arguments
    data.n_indiv = atoll(argv[1]);
    char *pheno_prefix = argv[2];
    char *grm_prefix = argv[3];
    
    printf("phenotype file: %s\n", pheno_prefix);
    printf("grm file:       %s\n", grm_prefix);
    
    // Allocate memory
    data.y = (float *)malloc(data.n_indiv * sizeof(float));
    data.na = (int *)malloc(data.n_indiv * sizeof(int));
    data.ids = (char **)malloc(data.n_indiv * sizeof(char *));
    
    for (long long i = 0; i < data.n_indiv; i++) {
        data.ids[i] = (char *)malloc(256 * sizeof(char));
    }
    
    setup_bins(&data);
    
    data.covariance_n = (long long *)malloc(data.n_bins * sizeof(long long));
    data.covariance = (double *)malloc(data.n_bins * sizeof(double));
    data.avg_bins = (double *)malloc(data.n_bins * sizeof(double));
    
    if (read_phenotype(&data, pheno_prefix) != 0) {
        return 1;
    }
    
    if (read_grm_and_calculate(&data, grm_prefix) != 0) {
        return 1;
    }
    
    // Calculate averages and output
    long long total_comparisons = 0;
    for (int k = 0; k < data.n_bins; k++) {
        total_comparisons += data.covariance_n[k];
    }
    
    printf("...a total of %lld comparisons\n", total_comparisons);
    
    char output_file[500];
    sprintf(output_file, "phenCovariance_%s.txt", pheno_prefix);
    FILE *out = fopen(output_file, "w");
    
    for (int k = 0; k < data.n_bins; k++) {
        if (data.covariance_n[k] > 0) {
            double avg_cov = data.covariance[k] / data.covariance_n[k];
            double avg_bin = data.avg_bins[k] / data.covariance_n[k];
            
            printf("%5d %10.4g %10.4g %10.4g %16.10g %16lld\n", 
                   k, data.bins[k], data.bins[k+1], avg_bin, avg_cov, data.covariance_n[k]);
            fprintf(out, "%5d %10.4g %10.4g %10.4g %16.10g %16lld\n", 
                   k, data.bins[k], data.bins[k+1], avg_bin, avg_cov, data.covariance_n[k]);
        }
    }
    
    fclose(out);
    
    // Free memory
    free(data.y);
    free(data.na);
    free(data.covariance_n);
    free(data.covariance);
    free(data.avg_bins);
    free(data.bins);
    
    for (long long i = 0; i < data.n_indiv; i++) {
        free(data.ids[i]);
    }
    free(data.ids);
    
    return 0;
}
