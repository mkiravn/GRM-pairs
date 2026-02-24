/*
 * main.c -- test harness for calc_rel_pairs()
 *
 * Reads .fam, .bim, and a frequency file directly. No subprocess calls.
 *
 * Usage: grm_pairs <bfile_prefix> <pair_file> <out_file>
 *
 * Expects:
 *   <bfile>.fam   -- plink fam file (FID IID PAT MAT SEX PHENO)
 *   <bfile>.bed   -- plink bed file
 *   <bfile>.bim   -- plink bim file (CHR SNP CM BP A1 A2)
 *   <bfile>.afreq -- plink2 frequency file from --freq (tried first)
 *                    (#CHROM ID REF ALT1 ALT1_FREQ OBS_CT)
 * or
 *   <bfile>.frq   -- plink1.9 frequency file from --freq
 *                    (CHR SNP A1 A2 MAF NCHROBS)
 *
 * The allele reported in the frequency file is compared against the A1
 * column of the .bim.  If they differ, the frequency is flipped (1-p)
 * so that p always equals freq(bim_A1), consistent with the .bed encoding
 * (x = 0/1/2 copies of A1).
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
char** read_fam(const char* fam_fname, uint32_t* sample_ct_out)
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
    char fid[256], iid[256], rest[512];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '\n') continue;
        if (sscanf(line, "%255s %255s %511[^\n]", fid, iid, rest) < 2) continue;
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
/* Read .bim -- returns array of A1 allele strings, sets snp_ct       */
/*                                                                      */
/* .bim format: CHR SNP CM BP A1 A2                                   */
/* A1 is the allele whose copy count (0/1/2) is stored in the .bed.   */
/* ------------------------------------------------------------------ */
static char** read_bim_a1(const char* bim_fname, uint32_t* snp_ct_out)
{
    FILE* f = fopen(bim_fname, "r");
    if (!f) { perror(bim_fname); return NULL; }

    uint32_t n = 0;
    char line[1024];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '\0' && line[0] != '\n') n++;
    }
    rewind(f);

    char** a1s = malloc(n * sizeof(char*));
    if (!a1s) { fclose(f); return NULL; }

    uint32_t i = 0;
    char chr[32], snpid[64], cm[32], bp[32], a1[16], a2[16];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '\n') continue;
        if (sscanf(line, "%31s %63s %31s %31s %15s %15s",
                   chr, snpid, cm, bp, a1, a2) < 6) continue;
        a1s[i] = malloc(strlen(a1) + 1);
        if (!a1s[i]) {
            for (uint32_t j = 0; j < i; j++) free(a1s[j]);
            free(a1s); fclose(f); return NULL;
        }
        strcpy(a1s[i], a1);
        i++;
    }
    fclose(f);

    *snp_ct_out = i;
    fprintf(stderr, "Read %u SNPs from %s\n", i, bim_fname);
    return a1s;
}

/* ------------------------------------------------------------------ */
/* Read a frequency file -- returns per-SNP frequencies and the allele */
/* name for which the frequency was reported.                           */
/*                                                                      */
/* Format is auto-detected from the header:                            */
/*   .frq  (plink1.9): CHR SNP A1 A2 MAF NCHROBS                       */
/*     header does NOT start with '#'; the frequency column is MAF     */
/*     and the allele string stored in *alleles_out[k] is A1.           */
/*   .afreq (plink2):  #CHROM ID REF ALT1 ALT1_FREQ OBS_CT            */
/*     header starts with '#'; the frequency column is ALT1_FREQ       */
/*     and the allele string stored in *alleles_out[k] is ALT1.         */
/*                                                                      */
/* The caller should compare *alleles_out[k] with bim A1 and flip      */
/* (1-p) for SNPs where they differ; see check_and_flip_freqs().       */
/* ------------------------------------------------------------------ */
static double* read_freq_file(const char* freq_fname,
                               char***     alleles_out,
                               uint32_t*   snp_ct_out)
{
    FILE* f = fopen(freq_fname, "r");
    if (!f) { perror(freq_fname); return NULL; }

    char line[1024];
    if (!fgets(line, sizeof(line), f)) {
        fprintf(stderr, "Error: %s appears empty (no header)\n", freq_fname);
        fclose(f);
        return NULL;
    }

    /* plink2 .afreq header starts with '#'; plink1.9 .frq does not */
    int is_afreq = (line[0] == '#');

    /* count data lines */
    uint32_t n = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '\0' && line[0] != '\n') n++;
    }
    rewind(f);
    if (!fgets(line, sizeof(line), f)) { fclose(f); return NULL; } /* skip header */

    double* freqs   = malloc(n * sizeof(double));
    char**  alleles = malloc(n * sizeof(char*));
    if (!freqs || !alleles) {
        free(freqs); free(alleles); fclose(f); return NULL;
    }

    uint32_t i = 0;
    char col1[32], col2[64], a1_col[16], a2_col[16];
    double freq;
    int    obs;

    if (is_afreq) {
        /* .afreq format: #CHROM ID REF ALT PROVISIONAL_REF? ALT_FREQS OBS_CT */
        char provisional[8];
        while (i < n &&
               fscanf(f, "%31s %63s %15s %15s %7s %lf %d",
                      col1, col2, a1_col, a2_col, provisional, &freq, &obs) == 7) {
            freqs[i] = freq;
            /* frequency is for ALT (a2_col) */
            const char* counted = a2_col;
            alleles[i] = malloc(strlen(counted) + 1);
            if (!alleles[i]) {
                for (uint32_t j = 0; j < i; j++) free(alleles[j]);
                free(alleles); free(freqs); fclose(f); return NULL;
            }
            strcpy(alleles[i], counted);
            i++;
        }
    } else {
        /* .frq format: CHR SNP A1 A2 MAF NCHROBS */
        while (i < n &&
               fscanf(f, "%31s %63s %15s %15s %lf %d",
                      col1, col2, a1_col, a2_col, &freq, &obs) == 6) {
            freqs[i] = freq;
            /* frequency is for A1 (a1_col) */
            const char* counted = a1_col;
            alleles[i] = malloc(strlen(counted) + 1);
            if (!alleles[i]) {
                for (uint32_t j = 0; j < i; j++) free(alleles[j]);
                free(alleles); free(freqs); fclose(f); return NULL;
            }
            strcpy(alleles[i], counted);
            i++;
        }
    }
    fclose(f);

    fprintf(stderr, "Read %u SNPs from %s (%s format)\n",
            i, freq_fname, is_afreq ? "plink2 .afreq" : "plink1.9 .frq");
    *alleles_out = alleles;
    *snp_ct_out  = i;
    return freqs;
}

/* ------------------------------------------------------------------ */
/* check_and_flip_freqs -- compare frequency-file alleles against bim  */
/* A1 and flip (p -> 1-p) for any SNP where they differ.              */
/*                                                                      */
/* The .bed encodes x = 0/1/2 copies of bim_A1.  The GRM formula      */
/* needs p = freq(bim_A1).  If the frequency file reported freq for    */
/* bim_A2 instead, p must be flipped so the centering is correct.      */
/* ------------------------------------------------------------------ */
static void check_and_flip_freqs(double*   freqs,
                                  char**    freq_alleles,
                                  char**    bim_a1,
                                  uint32_t  snp_ct)
{
    uint32_t n_flipped = 0;
    for (uint32_t k = 0; k < snp_ct; k++) {
        if (strcmp(freq_alleles[k], bim_a1[k]) != 0) {
            freqs[k] = 1.0 - freqs[k];
            n_flipped++;
        }
    }
    if (n_flipped > 0)
        fprintf(stderr, "Flipped allele coding for %u / %u SNPs\n",
                n_flipped, snp_ct);
}

/* ------------------------------------------------------------------ */
/* main                                                                 */
/* ------------------------------------------------------------------ */
int main(int argc, char** argv)
{
    if (argc < 4) {
        fprintf(stderr, "usage: grm_pairs <bfile> <pair_file> <out>\n");
        fprintf(stderr,
                "  expects <bfile>.bed, <bfile>.bim, <bfile>.fam and\n"
                "          <bfile>.afreq (plink2) or <bfile>.frq (plink1.9)\n");
        return 1;
    }

    const char* bfile     = argv[1];
    const char* pair_file = argv[2];
    const char* out_fname = argv[3];

    /* build file paths from bfile prefix */
    char bed_fname[512], fam_fname[512], bim_fname[512];
    char afreq_fname[512], frq_fname[512];
    snprintf(bed_fname,   sizeof(bed_fname),   "%s.bed",   bfile);
    snprintf(fam_fname,   sizeof(fam_fname),   "%s.fam",   bfile);
    snprintf(bim_fname,   sizeof(bim_fname),   "%s.bim",   bfile);
    snprintf(afreq_fname, sizeof(afreq_fname), "%s.afreq", bfile);
    snprintf(frq_fname,   sizeof(frq_fname),   "%s.frq",   bfile);

    /* read sample IDs from .fam */
    uint32_t sample_ct = 0;
    char** sample_ids = read_fam(fam_fname, &sample_ct);
    if (!sample_ids || sample_ct == 0) {
        fprintf(stderr, "Error reading %s\n", fam_fname);
        return 1;
    }

    /* auto-detect frequency file: .afreq (plink2) takes precedence over .frq */
    const char* freq_fname = NULL;
    {
        FILE* probe = fopen(afreq_fname, "r");
        if (probe) { fclose(probe); freq_fname = afreq_fname; }
        else        {               freq_fname = frq_fname;   }
    }

    uint32_t snp_ct       = 0;
    char**   freq_alleles = NULL;
    double*  allele_freqs = read_freq_file(freq_fname, &freq_alleles, &snp_ct);
    if (!allele_freqs || snp_ct == 0) {
        fprintf(stderr, "Error reading frequency file %s\n", freq_fname);
        for (uint32_t i = 0; i < sample_ct; i++) free(sample_ids[i]);
        free(sample_ids);
        return 1;
    }

    /* read .bim A1 alleles for the allele-flip check */
    uint32_t bim_snp_ct = 0;
    char**   bim_a1     = read_bim_a1(bim_fname, &bim_snp_ct);
    if (!bim_a1 || bim_snp_ct != snp_ct) {
        fprintf(stderr,
                "Error reading %s or SNP count mismatch "
                "(%u in .bim vs %u in freq file)\n",
                bim_fname, bim_snp_ct, snp_ct);
        for (uint32_t i = 0; i < sample_ct; i++) free(sample_ids[i]);
        free(sample_ids);
        for (uint32_t i = 0; i < snp_ct; i++) free(freq_alleles[i]);
        free(freq_alleles);
        free(allele_freqs);
        if (bim_a1) {
            for (uint32_t i = 0; i < bim_snp_ct; i++) free(bim_a1[i]);
            free(bim_a1);
        }
        return 1;
    }

    /* flip p -> 1-p for any SNP where the frequency file's allele != bim A1 */
    check_and_flip_freqs(allele_freqs, freq_alleles, bim_a1, snp_ct);

    /* free temporary allele string arrays */
    for (uint32_t i = 0; i < snp_ct; i++) {
        free(freq_alleles[i]);
        free(bim_a1[i]);
    }
    free(freq_alleles);
    free(bim_a1);

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
