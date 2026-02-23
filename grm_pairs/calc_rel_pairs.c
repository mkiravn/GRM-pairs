/*
 * calc_rel_pairs() -- sparse GRM for a user-specified list of individual pairs
 *
 * Structural model:
 *   - SNP streaming loop from plink 1.9's calc_rel() in plink_rel.c
 *   - Pair-list dispatch pattern from plink 2.0's CalcKingTableSubset()
 *     in plink2_matrix_calc.cc
 *
 * The .bed format stores genotypes as 2 bits per sample, 4 samples per byte,
 * in variant-major order (one row per SNP, all samples). The encoding is:
 *   00 = homozygous A1  (coded as 0)
 *   01 = missing        (coded as MISSING sentinel)
 *   10 = heterozygous   (coded as 1)
 *   11 = homozygous A2  (coded as 2)
 *
 * GRM entry for pair (i,j):
 *   A_ij = [sum_k w_ik * w_jk] / n_valid_ij
 * where
 *   w_ik = (x_ik - 2*p_k) / sqrt(2 * p_k * (1 - p_k))
 * and n_valid_ij counts SNPs where neither individual is missing.
 * This is the Yang et al. 2010 (GCTA) variance-standardized estimator.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define MISSING_GENO  255   /* sentinel for missing genotype after decoding */

/* ------------------------------------------------------------------ */
/* Data structures                                                      */
/* ------------------------------------------------------------------ */

typedef struct {
    uint32_t idx1;       /* 0-based index into sample array */
    uint32_t idx2;
    double   numerator;  /* running sum of w_i * w_j across SNPs */
    uint32_t n_valid;    /* SNPs where neither sample is missing  */
    double   grm;        /* final value filled in after SNP loop  */
} RelPair;

/* ------------------------------------------------------------------ */
/* Step 1: load the pair file                                           */
/*                                                                      */
/* Pair file format (whitespace-delimited):                             */
/*   IID1  IID2                                                         */
/* or with family IDs:                                                  */
/*   FID1  IID1  FID2  IID2                                             */
/*                                                                      */
/* Returns number of pairs loaded, or -1 on error.                      */
/* ------------------------------------------------------------------ */
int load_pair_file(const char*   pair_fname,
                   char**        sample_ids,   /* array of IID strings  */
                   uint32_t      sample_ct,
                   RelPair**     pairs_out)
{
    FILE* f = fopen(pair_fname, "r");
    if (!f) { perror(pair_fname); return -1; }

    /* two-pass: first count lines, then allocate */
    int pair_ct = 0;
    char line[1024];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '\0' && line[0] != '#') pair_ct++;
    }
    rewind(f);

    RelPair* pairs = calloc(pair_ct, sizeof(RelPair));
    if (!pairs) { fclose(f); return -1; }

    /* helper: linear scan to resolve an IID to its sample index */
    /* In production plink code this uses a hash map (id_htable)  */
    #define FIND_SAMPLE(name, out_idx) do {                         \
        int found = 0;                                              \
        for (uint32_t s = 0; s < sample_ct; s++) {                 \
            if (strcmp(sample_ids[s], (name)) == 0) {              \
                (out_idx) = s; found = 1; break;                   \
            }                                                       \
        }                                                           \
        if (!found) {                                               \
            fprintf(stderr, "Error: sample '%s' not in .fam\n", (name)); \
            fclose(f); free(pairs); return -1;                     \
        }                                                           \
    } while(0)

    int p = 0;
    char id1[256], id2[256];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '\0' || line[0] == '#') continue;
        if (sscanf(line, "%255s %255s", id1, id2) != 2) continue;

        uint32_t idx1, idx2;
        FIND_SAMPLE(id1, idx1);
        FIND_SAMPLE(id2, idx2);

        /* canonicalise: always store lower index first */
        if (idx1 > idx2) { uint32_t tmp = idx1; idx1 = idx2; idx2 = tmp; }

        pairs[p++] = (RelPair){ idx1, idx2, 0.0, 0, 0.0 };
    }
    fclose(f);

    *pairs_out = pairs;
    return pair_ct;
}

/* ------------------------------------------------------------------ */
/* Step 2: precompute the set of individuals actually needed            */
/*                                                                      */
/* Avoids decoding every individual's genotype from the .bed row when  */
/* only a small fraction of samples appear in any pair.                 */
/* ------------------------------------------------------------------ */
void build_needed_mask(const RelPair* pairs, uint32_t pair_ct,
                       uint32_t sample_ct, uint8_t* needed)
{
    memset(needed, 0, sample_ct);
    for (uint32_t p = 0; p < pair_ct; p++) {
        needed[pairs[p].idx1] = 1;
        needed[pairs[p].idx2] = 1;
    }
}

/* ------------------------------------------------------------------ */
/* Step 3: decode one SNP's .bed row into a uint8_t genotype array     */
/*                                                                      */
/* .bed encoding (2 bits per sample):                                   */
/*   00 -> 0 (hom A1), 01 -> MISSING, 10 -> 1 (het), 11 -> 2 (hom A2) */
/*                                                                      */
/* Only decodes samples flagged in needed[]; others left as 0.          */
/* ------------------------------------------------------------------ */
void decode_bed_row(const uint8_t* bed_row,  /* raw bytes for this SNP */
                    uint32_t       sample_ct,
                    const uint8_t* needed,
                    uint8_t*       geno)      /* output, length sample_ct */
{
        static const uint8_t bed_to_geno[4] = {
        0,            /* 00 -> hom A1 -> 2 copies of A1 */
        MISSING_GENO, /* 01 -> missing */
        1,            /* 10 -> het     -> 1 copy of A1  */
        2             /* 11 -> hom A2  -> 0 copies of A1 */
    };

    for (uint32_t i = 0; i < sample_ct; i++) {
        if (!needed[i]) continue;
        uint8_t byte  = bed_row[i / 4];
        uint8_t shift = (i % 4) * 2;
        geno[i] = bed_to_geno[(byte >> shift) & 0x3];
    }
}

/* ------------------------------------------------------------------ */
/* Step 4: the main function                                            */
/* ------------------------------------------------------------------ */
int calc_rel_pairs(const char*  bed_fname,
                   const char*  pair_fname,
                   char**       sample_ids,   /* from .fam, length sample_ct */
                   uint32_t     sample_ct,
                   double*      allele_freqs, /* per-SNP freq, length snp_ct */
                   uint32_t     snp_ct,
                   const char*  out_fname)
{
    /* --- load pair list ------------------------------------------- */
    RelPair* pairs = NULL;
    int pair_ct = load_pair_file(pair_fname, sample_ids, sample_ct, &pairs);
    if (pair_ct < 0) return -1;
    fprintf(stderr, "Loaded %d pairs\n", pair_ct);

    /* --- precompute needed-sample mask ----------------------------- */
    uint8_t* needed = calloc(sample_ct, 1);
    build_needed_mask(pairs, pair_ct, sample_ct, needed);

    /* --- allocate per-SNP genotype buffer -------------------------- */
    uint8_t* geno = calloc(sample_ct, 1);

    /* bytes per SNP row in .bed (rounded up to whole bytes) */
    uint32_t row_bytes = (sample_ct + 3) / 4;
    uint8_t* bed_row   = malloc(row_bytes);

    /* --- open .bed ------------------------------------------------- */
    FILE* bed = fopen(bed_fname, "rb");
    if (!bed) { perror(bed_fname); return -1; }

    /* skip 3-byte .bed magic header */
    fseek(bed, 3, SEEK_SET);

    /* ----------------------------------------------------------------
     * Main SNP loop
     *
     * For each SNP k:
     *   1. Read the raw .bed row (all N samples packed as 2-bit codes)
     *   2. Decode genotypes for needed individuals only
     *   3. Compute the variance-standardised weight w_ik for each sample
     *   4. Accumulate w_i * w_j for every pair
     * ---------------------------------------------------------------- */
    for (uint32_t k = 0; k < snp_ct; k++) {

        if (fread(bed_row, 1, row_bytes, bed) != row_bytes) {
            fprintf(stderr, "Error: truncated .bed at SNP %u\n", k);
            return -1;
        }

        double p  = allele_freqs[k];
        double denom = sqrt(2.0 * p * (1.0 - p));

        /* skip monomorphic or degenerate SNPs */
        if (denom < 1e-8) continue;

        decode_bed_row(bed_row, sample_ct, needed, geno);

        /* precompute normalised weights for needed individuals */
        double* w = malloc(sample_ct * sizeof(double));
        for (uint32_t i = 0; i < sample_ct; i++) {
            if (!needed[i] || geno[i] == MISSING_GENO) {
                w[i] = NAN;
            } else {
                w[i] = (geno[i] - 2.0 * p) / denom;
            }
        }

        /* accumulate for each pair */
        for (uint32_t p_idx = 0; p_idx < (uint32_t)pair_ct; p_idx++) {
            uint32_t i = pairs[p_idx].idx1;
            uint32_t j = pairs[p_idx].idx2;

            if (isnan(w[i]) || isnan(w[j])) continue;
            pairs[p_idx].numerator += w[i] * w[j];
            pairs[p_idx].n_valid   += 1;
        }

        free(w);

    }
    fclose(bed);

    /* --- finalise GRM entries -------------------------------------- */
    for (int p_idx = 0; p_idx < pair_ct; p_idx++) {
        if (pairs[p_idx].n_valid > 0) {
            pairs[p_idx].grm = pairs[p_idx].numerator / pairs[p_idx].n_valid;
        } else {
            pairs[p_idx].grm = NAN;  /* no valid SNPs for this pair */
        }
    }

    /* --- write output ---------------------------------------------- */
    /*
     * Output format mirrors plink --make-grm-list:
     *   IID1  IID2  N_valid_SNPs  GRM
     */
    FILE* out = fopen(out_fname, "w");
    if (!out) { perror(out_fname); return -1; }

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
    free(bed_row);
    return 0;
}

/*
 * NOTES ON PRODUCTION HARDENING
 * ==============================
 *
 * 1. PERFORMANCE: the w[] malloc/free inside the SNP loop is wasteful.
 *    Allocate once outside the loop. Also, precompute a compact array of
 *    needed indices (uint32_t needed_idxs[]) so the w-computation loop
 *    iterates over maybe 1000 individuals instead of all 500,000.
 *
 * 2. THREADING: the pair accumulation loop is embarrassingly parallel.
 *    In plink 1.9 style, split pairs across threads each with their own
 *    partial numerator/n_valid arrays, then sum at the end. The .bed
 *    reading must stay single-threaded (or use pgen reader in 2.0).
 *
 * 3. ID RESOLUTION: replace the O(N) linear scan in FIND_SAMPLE with
 *    plink's id_htable hash map (plink_common.c: populate_id_htable()).
 *    Also handle FID+IID pairs properly for samples with duplicate IIDs.
 *
 * 4. ALLELE FREQUENCIES: allele_freqs[] is assumed precomputed (plink
 *    does this in a prior pass for all SNPs). Optionally allow per-pair
 *    frequencies computed only from the two individuals (not standard).
 *
 * 5. MISSING DATA: plink's default for --make-rel is mean imputation
 *    (substitute 0 for missing, i.e., the mean of the centred column).
 *    The sketch above instead excludes missing SNPs per pair and adjusts
 *    the denominator, which matches --make-grm-list's reported N column.
 *
 * 6. X CHROMOSOME: plink applies hemizygosity correction for males on
 *    chrX. Skip chrX SNPs or add the correction if needed.
 */
