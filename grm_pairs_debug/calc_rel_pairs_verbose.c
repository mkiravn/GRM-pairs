/* calc_rel_pairs.c -- sparse GRM for a user-specified list of individual pairs
 *
 * The PLINK1 .bed format stores genotypes as 2 bits per sample, 4 samples per byte,
 * in variant-major order (one row per SNP, all samples). The 2-bit encoding is:
 *   00 = homozygous A1
 *   01 = missing
 *   10 = heterozygous
 *   11 = homozygous A2
 *
 * We decode each genotype into "A1 allele count" x in {0,1,2}:
 *   hom A1 -> 2, het -> 1, hom A2 -> 0, missing -> MISSING_GENO
 *
 * GRM entry for pair (i,j):
 *   A_ij = [sum_k w_ik * w_jk] / n_valid_ij
 * where
 *   w_ik = (x_ik - 2*p_k) / sqrt(2 * p_k * (1 - p_k))
 * and n_valid_ij counts SNPs where neither individual is missing.
 *
 * Debugging:
 *   This file includes an optional verbose debug mode which prints, for the first
 *   5 SNPs and first 5 individuals:
 *     - decoded genotypes
 *     - p, denom
 *     - weights
 *     - contributions for the first 5 pairs
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define MISSING_GENO 255

/* set to 1 to print debug for first 5 SNPs/inds/pairs */
#ifndef GRM_VERBOSE_DEBUG
#define GRM_VERBOSE_DEBUG 1
#endif

typedef struct {
    uint32_t idx1;
    uint32_t idx2;
    double   numerator;
    uint32_t n_valid;
    double   grm;
} RelPair;

/* ------------------------------------------------------------------ */
/* Pair file loader                                                    */
/*   Accepts either:                                                   */
/*     IID1 IID2                                                       */
/*   or:                                                               */
/*     FID1 IID1 FID2 IID2                                             */
/* We match on IID only (consistent with main.c's sample_ids).         */
/* ------------------------------------------------------------------ */
static int load_pair_file(const char* pair_fname,
                          char**      sample_ids,
                          uint32_t    sample_ct,
                          RelPair**   pairs_out)
{
    FILE* f = fopen(pair_fname, "r");
    if (!f) { perror(pair_fname); return -1; }

    /* count non-comment nonblank lines */
    int pair_ct = 0;
    char line[1024];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '\n' || line[0] == '#') continue;
        pair_ct++;
    }
    rewind(f);

    RelPair* pairs = calloc((size_t)pair_ct, sizeof(RelPair));
    if (!pairs) { fclose(f); return -1; }

    /* helper: resolve IID to sample index (O(N) linear scan) */
    #define FIND_SAMPLE(iid_name, out_idx) do {                         \
        int found = 0;                                                  \
        for (uint32_t s = 0; s < sample_ct; s++) {                      \
            if (strcmp(sample_ids[s], (iid_name)) == 0) {               \
                (out_idx) = s; found = 1; break;                        \
            }                                                           \
        }                                                               \
        if (!found) {                                                   \
            fprintf(stderr, "Error: sample IID '%s' not in .fam\n", (iid_name)); \
            fclose(f); free(pairs); return -1;                          \
        }                                                               \
    } while(0)

    int p = 0;
    char t1[256], t2[256], t3[256], t4[256];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '\n' || line[0] == '#') continue;

        int n_tok = sscanf(line, "%255s %255s %255s %255s", t1, t2, t3, t4);
        if (n_tok < 2) continue;

        /* If 4 tokens: FID1 IID1 FID2 IID2 -> use IID1=t2, IID2=t4
           If 2 tokens: IID1 IID2 -> use t1,t2 */
        const char* id1 = (n_tok >= 4) ? t2 : t1;
        const char* id2 = (n_tok >= 4) ? t4 : t2;

        uint32_t idx1, idx2;
        FIND_SAMPLE(id1, idx1);
        FIND_SAMPLE(id2, idx2);

        if (idx1 > idx2) { uint32_t tmp = idx1; idx1 = idx2; idx2 = tmp; }

        pairs[p++] = (RelPair){ idx1, idx2, 0.0, 0, 0.0 };
    }
    fclose(f);

    *pairs_out = pairs;
    return pair_ct;
}

static void build_needed_mask(const RelPair* pairs, uint32_t pair_ct,
                              uint32_t sample_ct, uint8_t* needed)
{
    memset(needed, 0, sample_ct);
    for (uint32_t p = 0; p < pair_ct; p++) {
        needed[pairs[p].idx1] = 1;
        needed[pairs[p].idx2] = 1;
    }
}

static void decode_bed_row(const uint8_t* bed_row,
                           uint32_t       sample_ct,
                           const uint8_t* needed,
                           uint8_t*       geno)
{
    /* Map PLINK2-bit code -> A1 allele count */
    static const uint8_t bed_to_a1count[4] = {
        2,            /* 00 -> hom A1 -> 2 copies of A1 */
        MISSING_GENO, /* 01 -> missing */
        1,            /* 10 -> het     -> 1 copy of A1 */
        0             /* 11 -> hom A2  -> 0 copies of A1 */
    };

    for (uint32_t i = 0; i < sample_ct; i++) {
        if (!needed[i]) continue;
        uint8_t byte  = bed_row[i / 4];
        uint8_t shift = (uint8_t)((i % 4) * 2);
        geno[i] = bed_to_a1count[(byte >> shift) & 0x3];
    }
}

/* debug helpers */
#if GRM_VERBOSE_DEBUG
static const char* geno_str(uint8_t g, char buf[8]) {
    if (g == MISSING_GENO) return "NA";
    snprintf(buf, 8, "%u", (unsigned)g);
    return buf;
}
#endif

int calc_rel_pairs(const char*  bed_fname,
                   const char*  pair_fname,
                   char**       sample_ids,
                   uint32_t     sample_ct,
                   double*      allele_freqs,
                   uint32_t     snp_ct,
                   const char*  out_fname)
{
    RelPair* pairs = NULL;
    int pair_ct = load_pair_file(pair_fname, sample_ids, sample_ct, &pairs);
    if (pair_ct < 0) return -1;
    fprintf(stderr, "Loaded %d pairs\n", pair_ct);

    uint8_t* needed = calloc(sample_ct, 1);
    uint8_t* geno   = calloc(sample_ct, 1);
    double*  w      = malloc(sample_ct * sizeof(double));
    if (!needed || !geno || !w) {
        fprintf(stderr, "Error: out of memory\n");
        free(pairs); free(needed); free(geno); free(w);
        return -1;
    }
    build_needed_mask(pairs, (uint32_t)pair_ct, sample_ct, needed);

    uint32_t row_bytes = (sample_ct + 3) / 4;
    uint8_t* bed_row   = malloc(row_bytes);
    if (!bed_row) {
        fprintf(stderr, "Error: out of memory\n");
        free(pairs); free(needed); free(geno); free(w);
        return -1;
    }

    FILE* bed = fopen(bed_fname, "rb");
    if (!bed) {
        perror(bed_fname);
        free(pairs); free(needed); free(geno); free(w); free(bed_row);
        return -1;
    }

    /* skip 3-byte .bed header */
    if (fseek(bed, 3, SEEK_SET) != 0) {
        fprintf(stderr, "Error: could not seek past .bed header\n");
        fclose(bed);
        free(pairs); free(needed); free(geno); free(w); free(bed_row);
        return -1;
    }

#if GRM_VERBOSE_DEBUG
    const uint32_t DBG_SNP  = (snp_ct < 5 ? snp_ct : 5);
    const uint32_t DBG_IND  = (sample_ct < 5 ? sample_ct : 5);
    const uint32_t DBG_PAIR = ((uint32_t)pair_ct < 5 ? (uint32_t)pair_ct : 5);

    fprintf(stderr, "\n[DEBUG] First %u individuals:\n", DBG_IND);
    for (uint32_t i = 0; i < DBG_IND; i++) {
        fprintf(stderr, "  i=%u id=%s\n", i, sample_ids[i]);
    }

    fprintf(stderr, "[DEBUG] First %u pairs:\n", DBG_PAIR);
    for (uint32_t p_idx = 0; p_idx < DBG_PAIR; p_idx++) {
        uint32_t i = pairs[p_idx].idx1, j = pairs[p_idx].idx2;
        fprintf(stderr, "  p=%u (%u,%u) %s %s\n", p_idx, i, j, sample_ids[i], sample_ids[j]);
    }
#endif

    for (uint32_t k = 0; k < snp_ct; k++) {
        if (fread(bed_row, 1, row_bytes, bed) != row_bytes) {
            fprintf(stderr, "Error: truncated .bed at SNP %u\n", k);
            fclose(bed);
            free(pairs); free(needed); free(geno); free(w); free(bed_row);
            return -1;
        }

        double p = allele_freqs[k];
        double denom = sqrt(2.0 * p * (1.0 - p));
        if (denom < 1e-8) {
#if GRM_VERBOSE_DEBUG
            if (k < DBG_SNP) {
                fprintf(stderr, "\n[DEBUG] SNP k=%u p=%.8g denom=%.8g (SKIP denom small)\n", k, p, denom);
            }
#endif
            continue;
        }

        decode_bed_row(bed_row, sample_ct, needed, geno);

        for (uint32_t i = 0; i < sample_ct; i++) {
            if (!needed[i] || geno[i] == MISSING_GENO) w[i] = NAN;
            else w[i] = ((double)geno[i] - 2.0 * p) / denom;
        }

#if GRM_VERBOSE_DEBUG
        if (k < DBG_SNP) {
            fprintf(stderr, "\n[DEBUG] SNP k=%u  p=%.10g  denom=%.10g\n", k, p, denom);

            fprintf(stderr, "[DEBUG]   geno (first %u): ", DBG_IND);
            for (uint32_t i = 0; i < DBG_IND; i++) {
                char buf[8];
                fprintf(stderr, "%s%s", geno_str(geno[i], buf), (i + 1 == DBG_IND) ? "" : " ");
            }
            fprintf(stderr, "\n");

            fprintf(stderr, "[DEBUG]   w    (first %u): ", DBG_IND);
            for (uint32_t i = 0; i < DBG_IND; i++) {
                if (isnan(w[i])) fprintf(stderr, "NA%s", (i + 1 == DBG_IND) ? "" : " ");
                else fprintf(stderr, "%.6g%s", w[i], (i + 1 == DBG_IND) ? "" : " ");
            }
            fprintf(stderr, "\n");

            fprintf(stderr, "[DEBUG]   contributions (first %u pairs):\n", DBG_PAIR);
            for (uint32_t p_idx = 0; p_idx < DBG_PAIR; p_idx++) {
                uint32_t i = pairs[p_idx].idx1;
                uint32_t j = pairs[p_idx].idx2;
                int valid = (!isnan(w[i]) && !isnan(w[j]));
                if (!valid) {
                    fprintf(stderr, "    p=%u (%s,%s): missing -> skip\n",
                            p_idx, sample_ids[i], sample_ids[j]);
                } else {
                    double contrib = w[i] * w[j];
                    fprintf(stderr, "    p=%u (%s,%s): w_i=%.8g w_j=%.8g contrib=%.8g\n",
                            p_idx, sample_ids[i], sample_ids[j], w[i], w[j], contrib);
                }
            }
        }
#endif

        for (uint32_t p_idx = 0; p_idx < (uint32_t)pair_ct; p_idx++) {
            uint32_t i = pairs[p_idx].idx1;
            uint32_t j = pairs[p_idx].idx2;
            if (isnan(w[i]) || isnan(w[j])) continue;
            pairs[p_idx].numerator += w[i] * w[j];
            pairs[p_idx].n_valid   += 1;
        }
    }

    fclose(bed);

    for (int p_idx = 0; p_idx < pair_ct; p_idx++) {
        pairs[p_idx].grm = (pairs[p_idx].n_valid > 0)
            ? (pairs[p_idx].numerator / (double)pairs[p_idx].n_valid)
            : NAN;
    }

#if GRM_VERBOSE_DEBUG
    {
        const uint32_t DBG_PAIR = ((uint32_t)pair_ct < 5 ? (uint32_t)pair_ct : 5);
        fprintf(stderr, "\n[DEBUG] Final summary (first %u pairs):\n", DBG_PAIR);
        for (uint32_t p_idx = 0; p_idx < DBG_PAIR; p_idx++) {
            fprintf(stderr, "  p=%u (%s,%s): numerator=%.12g n_valid=%u grm=%.12g\n",
                    p_idx,
                    sample_ids[pairs[p_idx].idx1],
                    sample_ids[pairs[p_idx].idx2],
                    pairs[p_idx].numerator,
                    pairs[p_idx].n_valid,
                    pairs[p_idx].grm);
        }
    }
#endif

    FILE* out = fopen(out_fname, "w");
    if (!out) {
        perror(out_fname);
        free(pairs); free(needed); free(geno); free(w); free(bed_row);
        return -1;
    }

    fprintf(out, "IID1\tIID2\tN_SNPs\tGRM\n");
    for (int p_idx = 0; p_idx < pair_ct; p_idx++) {
        fprintf(out, "%s\t%s\t%u\t%.8g\n",
                sample_ids[pairs[p_idx].idx1],
                sample_ids[pairs[p_idx].idx2],
                pairs[p_idx].n_valid,
                pairs[p_idx].grm);
    }
    fclose(out);

    free(pairs);
    free(needed);
    free(geno);
    free(w);
    free(bed_row);
    return 0;
}