#!/usr/bin/env python3
"""Manufacture a dense .grm.bin/.grm.id from grm_pairs pairwise output.

We don't have a working plink/gcta build on this machine to produce a real
dense GRM, so this stitches one together from the already-validated
grm_pairs sparse pairwise tool: ask it for every pair among a small subset
of individuals, then lay the results out in the lower-triangular float32
layout blocked_grm_bin_means / grm_shard_tool expect. Diagonal entries are
never read by either tool (they skip i==j), so they're filled with a
placeholder.

Usage:
    build_fake_grm.py <grm_pairs_output.tsv> <ids_file> <bfile_prefix> <out_prefix>

<ids_file>: one IID per line, in the desired row/col order for the matrix.
<bfile_prefix>: used only to read real FIDs from <bfile_prefix>.fam, so
grm.id's (FID, IID) pairs match what simulate_pheno.py writes into the
phenotype file -- read_pheno_aligned() looks samples up by (FID, IID), so
a mismatched FID (e.g. reusing IID as FID here) silently drops every
phenotype.
Writes <out_prefix>.grm.bin and <out_prefix>.grm.id.
"""
import struct
import sys


def main() -> None:
    if len(sys.argv) != 5:
        sys.exit(__doc__)

    grm_pairs_out, ids_file, bfile_prefix, out_prefix = sys.argv[1:5]

    with open(ids_file) as f:
        ids = [line.strip() for line in f if line.strip()]
    n = len(ids)

    iid_to_fid = {}
    with open(f"{bfile_prefix}.fam") as f:
        for line in f:
            fid, iid = line.split()[:2]
            iid_to_fid[iid] = fid

    grm = {}
    with open(grm_pairs_out) as f:
        header = f.readline()
        for line in f:
            iid1, iid2, _n_snps, val = line.split()
            grm[frozenset((iid1, iid2))] = float(val)

    with open(f"{out_prefix}.grm.id", "w") as f:
        for iid in ids:
            f.write(f"{iid_to_fid[iid]}\t{iid}\n")

    with open(f"{out_prefix}.grm.bin", "wb") as f:
        for i in range(n):
            for j in range(i + 1):
                if i == j:
                    val = 1.0
                else:
                    val = grm[frozenset((ids[i], ids[j]))]
                f.write(struct.pack("<f", val))

    n_entries = n * (n + 1) // 2
    print(f"wrote {out_prefix}.grm.bin ({n_entries} entries, {n_entries * 4} bytes) "
          f"and {out_prefix}.grm.id ({n} ids)")


if __name__ == "__main__":
    main()
