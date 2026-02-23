/* main.c -- test harness for calc_rel_pairs()
 *
 * Reads .fam, .frq files directly. No subprocess calls.
 *
 * Usage: grm_pairs <bfile_prefix> <pair_file> <out_file>
 *
 * Expects:
 *   <bfile>.fam  -- plink fam file (FID IID PAT MAT SEX PHENO)
 *   <bfile>.bed  -- plink bed file
 *   <bfile>.frq  -- plink --freq output:
 *                   CHR SNP A1 A2 MAF NCHROBS
 *
 * Notes:
 * - For plink1 --freq, A1 is the (minor) allele and MAF is freq(A1) (<= 0.5).
 * - In calc_rel_pairs.c we decode genotypes into "A1 allele count" (0/1/2),
 *   so centering by 2*p where p=freq(A1) is consistent.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* forward declaration of the GRM function */
int calc_rel_pairs(const char* bed_fname,
                   const char* pair_fname,
                   char**      sample_ids,
                   uint32_t    sample_ct,
                   double*     allele_freqs,
                   uint32_t    snp_ct,
                   const char* out_fname);

/* ------------------------------------------------------------------ */
/* Read .fam -- returns array of IID strings, sets sample_ct          */
/* ------------------------------------------------------------------ */
static char** read_fam(const char* fam_fname, uint32_t* sample_ct_out)
{
    FILE* f = fopen(fam_fname, "r");
    if (!f) { perror(fam_fname); return NULL; }

    /* count lines */
    uint32_t n = 0;
    char line[1024];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '\0' && line[0] != '\n') n++;
    }
    rewind(f);

    char** ids = malloc(n * sizeof(char*));
    if (!ids) { fclose(f); return NULL; }

    uint32_t i = 0;
    char fid[256], iid[256];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '\n') continue;
        if (sscanf(line, "%255s %255s", fid, iid) < 2) continue;
        ids[i] = malloc(strlen(iid) + 1);
        if (!ids[i]) { fclose(f); return NULL; }
        strcpy(ids[i], iid);
        i++;
    }
    fclose(f);

    *sample_ct_out = i;
    fprintf(stderr, "Read %u samples from %s\n", i, fam_fname);
    return ids;
}

/* ------------------------------------------------------------------ */
/* Read .frq (plink --freq) -- returns array of freq(A1)=MAF, sets n  */
/* ------------------------------------------------------------------ */
static double* read_frq(const char* frq_fname, uint32_t* snp_ct_out)
{
    FILE* f = fopen(frq_fname, "r");
    if (!f) { perror(frq_fname); return NULL; }

    char line[1024];

    /* skip header (check return to avoid -Wunused-result warning) */
    if (fgets(line, sizeof(line), f) == NULL) {
        fprintf(stderr, "Error: %s appears empty (no header)\n", frq_fname);
        fclose(f);
        return NULL;
    }

    /* count remaining nonblank lines */
    uint32_t n = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '\0' && line[0] != '\n') n++;
    }

    rewind(f);
    if (fgets(line, sizeof(line), f) == NULL) { /* skip header again */
        fprintf(stderr, "Error: %s appears empty (no header on rewind)\n", frq_fname);
        fclose(f);
        return NULL;
    }

    double* freqs = malloc(n * sizeof(double));
    if (!freqs) { fclose(f); return NULL; }

    uint32_t i = 0;
    char chr[32], snpid[64], a1[16], a2[16];
    double maf;
    int nchrobs;

    /* Parse: CHR SNP A1 A2 MAF NCHROBS */
    while (fscanf(f, "%31s %63s %15s %15s %lf %d",
                  chr, snpid, a1, a2, &maf, &nchrobs) == 6) {
        freqs[i++] = maf;
        if (i == n) break;
    }

    fclose(f);

    *snp_ct_out = i;
    fprintf(stderr, "Read %u SNPs from %s\n", i, frq_fname);
    return freqs;
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */
int main(int argc, char** argv)
{
    if (argc < 4) {
        fprintf(stderr, "usage: grm_pairs <bfile> <pair_file> <out>\n");
        fprintf(stderr, "  expects <bfile>.fam, <bfile>.bed, <bfile>.frq\n");
        return 1;
    }

    const char* bfile     = argv[1];
    const char* pair_file = argv[2];
    const char* out_fname = argv[3];

    /* build file paths from bfile prefix */
    char bed_fname[512], fam_fname[512], frq_fname[512];
    snprintf(bed_fname, sizeof(bed_fname), "%s.bed", bfile);
    snprintf(fam_fname, sizeof(fam_fname), "%s.fam", bfile);
    snprintf(frq_fname, sizeof(frq_fname), "%s.frq", bfile);

    /* read sample IDs from .fam */
    uint32_t sample_ct = 0;
    char** sample_ids = read_fam(fam_fname, &sample_ct);
    if (!sample_ids || sample_ct == 0) {
        fprintf(stderr, "Error reading %s\n", fam_fname);
        return 1;
    }

    /* read allele frequencies from .frq */
    uint32_t snp_ct = 0;
    double* allele_freqs = read_frq(frq_fname, &snp_ct);
    if (!allele_freqs || snp_ct == 0) {
        fprintf(stderr, "Error reading %s\n", frq_fname);
        for (uint32_t i = 0; i < sample_ct; i++) free(sample_ids[i]);
        free(sample_ids);
        return 1;
    }

    fprintf(stderr, "Running calc_rel_pairs: %u samples, %u SNPs\n",
            sample_ct, snp_ct);

    int ret = calc_rel_pairs(bed_fname, pair_file, sample_ids,
                             sample_ct, allele_freqs, snp_ct, out_fname);

    /* cleanup */
    for (uint32_t i = 0; i < sample_ct; i++) free(sample_ids[i]);
    free(sample_ids);
    free(allele_freqs);

    if (ret == 0) fprintf(stderr, "Done. Output written to %s\n", out_fname);
    return ret;
}