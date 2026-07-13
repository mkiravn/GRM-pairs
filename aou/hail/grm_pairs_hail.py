#!/usr/bin/env python3
"""Rough, UNVERIFIED sketch of grm_pairs for Hail MatrixTable input.

Not yet run against real AoU data. Mirrors ../plink/run_chunk.sh: given one
chunk of a pair list, compute the Yang et al. 2010 (GCTA) variance-
standardised GRM estimator

    A_ij = (1/M) * sum_k [ (x_ik - 2*p_k)(x_jk - 2*p_k) / (2*p_k*(1-p_k)) ]

for just the individuals referenced by that chunk, without materialising the
full-cohort GRM -- same "needed samples only" idea as
../../grm_pairs/calc_rel_pairs.c build_needed_mask().

Approach: subset the MatrixTable to the samples the chunk actually needs
(usually a small fraction of the cohort), variant-standardise genotypes,
convert to a Hail BlockMatrix, and take the dense sample-by-sample product
X^T X / M restricted to that subset -- then look up each requested pair in
the resulting small dense matrix. This still computes a dense (small) GRM
for the chunk's needed samples rather than only the exact requested pairs;
for chunks where the needed-sample set is close to the full pair count this
is fine, but if a chunk references a small number of pairs scattered across
many distinct samples this wastes work relative to the plink path's
per-SNP, per-pair accumulation. Worth revisiting once real timings exist.

Known divergence from ../../grm_pairs (not yet reconciled): missing
genotypes are mean-imputed here (BlockMatrix.from_entry_expr mean_impute),
so every pair reports the same N_SNPs (all variants after the af filter),
whereas grm_pairs excludes missing calls per pair and reports the true
per-pair N_valid. Numbers from this script should not be assumed to match
grm_pairs output until that's reconciled.

Usage (illustrative):
    python3 grm_pairs_hail.py <mt_path> <pair_chunk_file> <out_file>
"""
import sys

import hail as hl


def load_pairs(path: str) -> list[tuple[str, str]]:
    pairs = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            iid1, iid2 = line.split()[:2]
            pairs.append((iid1, iid2))
    return pairs


def main() -> None:
    if len(sys.argv) != 4:
        sys.exit(__doc__)

    mt_path, pair_chunk_file, out_file = sys.argv[1:4]

    pairs = load_pairs(pair_chunk_file)
    needed_ids = sorted({iid for pair in pairs for iid in pair})

    hl.init(quiet=True)
    mt = hl.read_matrix_table(mt_path)
    mt = mt.filter_cols(hl.literal(set(needed_ids)).contains(mt.s))

    # Variant-standardise: (x - 2p) / sqrt(2p(1-p)), skip monomorphic sites.
    mt = mt.annotate_rows(af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2)
    mt = mt.filter_rows((mt.af > 0) & (mt.af < 1))
    mt = mt.annotate_entries(
        std_gt=(mt.GT.n_alt_alleles() - 2 * mt.af)
        / hl.sqrt(2 * mt.af * (1 - mt.af))
    )

    bm = hl.linalg.BlockMatrix.from_entry_expr(mt.std_gt, mean_impute=True)
    n_variants = bm.shape[0]
    grm_bm = (bm.T @ bm) / n_variants

    grm = grm_bm.to_numpy()
    col_ids = mt.s.collect()
    idx_of = {iid: i for i, iid in enumerate(col_ids)}

    with open(out_file, "w") as out:
        out.write("IID1\tIID2\tN_SNPs\tGRM\n")
        for iid1, iid2 in pairs:
            i, j = idx_of[iid1], idx_of[iid2]
            out.write(f"{iid1}\t{iid2}\t{n_variants}\t{grm[i, j]:.8g}\n")


if __name__ == "__main__":
    main()
