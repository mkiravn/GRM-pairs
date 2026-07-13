#!/usr/bin/env python3
"""Simulate a simple additive phenotype directly from a plink .bed file.

Stands in for `plink2 --score` (not runnable on this machine -- the local
plink2 build needs AVX2). Decodes genotypes for a chosen sample subset and a
random SNP subset, samples per-SNP effects ~ N(0, h2/M), computes a
variance-standardized additive score, and adds noise for the target h2.

Usage:
    simulate_pheno.py <bfile_prefix> <ids_file> <n_snps> <h2> <seed> <out_pheno>

<ids_file>: one IID per line, defining which samples to simulate for (order
does not need to match anything -- phenotypes are looked up by ID downstream).
"""
import struct
import sys

import numpy as np


def read_fam(bfile: str):
    fids, iids = [], []
    with open(f"{bfile}.fam") as f:
        for line in f:
            parts = line.split()
            fids.append(parts[0])
            iids.append(parts[1])
    return fids, iids


def read_bed(bfile: str, n_samples: int, snp_indices, needed_rows):
    """Decode genotype dosage (count of A2 allele, 0/1/2, NaN if missing) for
    the given SNP row indices and sample column indices."""
    row_bytes = (n_samples + 3) // 4
    code_to_dosage = np.array([0.0, np.nan, 1.0, 2.0])

    geno = np.empty((len(snp_indices), len(needed_rows)))
    with open(f"{bfile}.bed", "rb") as f:
        magic = f.read(3)
        if magic != b"\x6c\x1b\x01":
            raise ValueError("Not a SNP-major .bed file")

        for out_row, snp_idx in enumerate(snp_indices):
            f.seek(3 + snp_idx * row_bytes)
            raw = np.frombuffer(f.read(row_bytes), dtype=np.uint8)
            codes = np.empty(n_samples, dtype=np.uint8)
            for shift in range(4):
                codes[shift::4] = (raw[: (n_samples - shift + 3) // 4] >> (2 * shift)) & 0x3
            geno[out_row, :] = code_to_dosage[codes[needed_rows]]
    return geno


def main() -> None:
    if len(sys.argv) != 7:
        sys.exit(__doc__)

    bfile, ids_file, n_snps_str, h2_str, seed_str, out_pheno = sys.argv[1:7]
    n_snps = int(n_snps_str)
    h2 = float(h2_str)
    rng = np.random.default_rng(int(seed_str))

    fids, iids = read_fam(bfile)
    iid_to_row = {iid: i for i, iid in enumerate(iids)}

    with open(ids_file) as f:
        target_iids = [line.strip() for line in f if line.strip()]
    rows = np.array([iid_to_row[iid] for iid in target_iids])

    with open(f"{bfile}.bim") as f:
        n_total_snps = sum(1 for _ in f)
    snp_indices = rng.choice(n_total_snps, size=n_snps, replace=False)
    snp_indices.sort()

    geno = read_bed(bfile, len(iids), snp_indices, rows)  # (n_snps, n_target)

    freq = np.nanmean(geno, axis=1) / 2.0
    keep = (freq > 0.01) & (freq < 0.99)
    geno = geno[keep]
    freq = freq[keep]
    m = geno.shape[0]

    col_mean = 2.0 * freq[:, None]
    denom = np.sqrt(2.0 * freq * (1.0 - freq))[:, None]
    std_geno = (geno - col_mean) / denom
    std_geno = np.nan_to_num(std_geno, nan=0.0)  # mean-impute missing calls

    beta = rng.normal(0.0, np.sqrt(h2 / m), size=m)
    g = std_geno.T @ beta  # (n_target,)
    g = g / g.std() * np.sqrt(h2)  # renormalize: small m makes realized var(g) noisy

    e = rng.normal(0.0, np.sqrt(1.0 - h2), size=len(target_iids))
    y = g + e

    with open(out_pheno, "w") as out:
        for iid, val in zip(target_iids, y):
            fid = fids[iid_to_row[iid]]
            out.write(f"{fid}\t{iid}\t{val:.6f}\n")

    print(f"wrote {len(target_iids)} phenotypes using {m} SNPs (h2={h2}) to {out_pheno}")


if __name__ == "__main__":
    main()
