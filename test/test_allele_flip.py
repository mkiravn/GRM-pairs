#!/usr/bin/env python3
"""
test_allele_flip.py -- minimal end-to-end test for allele-flip robustness.

Verifies that grm_pairs produces identical GRM values regardless of whether
the frequency file uses plink1.9 .frq or plink2 .afreq format, and regardless
of which allele is labelled A1/ALT1 (i.e. whether a flip is needed).

Test layout
-----------
4 samples, 3 SNPs.

.bim (A1 is the "counted" allele in the .bed):
    1  rs0  0  1000  A  G
    1  rs1  0  2000  C  T
    1  rs2  0  3000  G  A

Genotypes encoded in .bed (x = 0=HOM_A1, 1=HET, 2=HOM_A2):
         rs0  rs1  rs2
    s0:   0    2    1
    s1:   1    0    2
    s2:   2    1    0
    s3:   0    2    1

Allele frequencies fed to the GRM formula (p = freq(bim_A1)):
    rs0: p = 0.30
    rs1: p = 0.50
    rs2: p = 0.40

Four frequency-file variants are tested -- all must produce the same GRM:
    1. frq_standard   -- .frq  with A1 = bim_A1 (no flip needed)
    2. frq_swapped    -- .frq  with A1 = bim_A2 (flip needed; MAF = 1-p)
    3. afreq_alt_a1   -- .afreq with ALT1 = bim_A1 (no flip needed)
    4. afreq_alt_a2   -- .afreq with ALT1 = bim_A2 (flip needed; ALT1_FREQ = 1-p)

Expected GRM values (computed analytically from the genotypes and p above,
using w_ik = (2 - x_ik - 2*p_k) / sqrt(2*p_k*(1-p_k)) where x_ik is the
.bed-decoded value, i.e. 0=HOM_A1, 1=HET, 2=HOM_A2):
    GRM(s0, s1) = -1/3   ≈ -0.33333333
    GRM(s0, s2) = -1/2   = -0.50000000
"""

import os
import shutil
import struct
import subprocess
import sys
import tempfile

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BINARY    = os.path.join(REPO_ROOT, "grm_pairs", "grm_pairs")

# ---------------------------------------------------------------------------
# Expected GRM values (exact fractions → floats)
# ---------------------------------------------------------------------------
EXPECTED = {
    ("s0", "s1"): -1.0 / 3.0,   # ≈ -0.33333
    ("s0", "s2"): -1.0 / 2.0,   # = -0.5
}
TOLERANCE = 1e-6


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_fam(path):
    """4 samples, no phenotype."""
    lines = [
        "0 s0 0 0 1 -9",
        "0 s1 0 0 1 -9",
        "0 s2 0 0 1 -9",
        "0 s3 0 0 1 -9",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def write_bim(path):
    """3 SNPs.  A1 is the allele counted as x=0/1/2 in the .bed."""
    lines = [
        "1\trs0\t0\t1000\tA\tG",
        "1\trs1\t0\t2000\tC\tT",
        "1\trs2\t0\t3000\tG\tA",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _geno_to_bed_code(x):
    """Map genotype count (0=HOM_A1, 1=HET, 2=HOM_A2) to 2-bit bed code."""
    return {0: 0b00, 1: 0b10, 2: 0b11}[x]


def write_bed(path):
    """
    Genotypes (x = copies of A1):
         rs0  rs1  rs2
    s0:   0    2    1
    s1:   1    0    2
    s2:   2    1    0
    s3:   0    2    1

    4 samples fit in 1 byte per SNP row (4 samples × 2 bits = 8 bits).
    Bit layout: bits[2i : 2i+2] = sample i's code.
    """
    genotypes = [
        [0, 2, 1],   # rs0: s0, s1, s2, s3
        [1, 0, 2],
        [2, 1, 0],
        [0, 2, 1],
    ]
    # Build SNP-major order: snp_genos[k] = list of 4 sample genotypes for SNP k
    snp_genos = [
        [genotypes[s][snp] for s in range(4)] for snp in range(3)
    ]

    with open(path, "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))  # plink .bed magic
        for row in snp_genos:
            byte = 0
            for s, x in enumerate(row):
                byte |= _geno_to_bed_code(x) << (2 * s)
            f.write(struct.pack("B", byte))


def write_pairs(path):
    with open(path, "w") as f:
        f.write("s0 s1\n")
        f.write("s0 s2\n")


# ---------------------------------------------------------------------------
# Frequency-file writers
# ---------------------------------------------------------------------------

# p = freq(bim_A1) that the formula should receive
P = [0.30, 0.50, 0.40]

BIM_A1 = ["A", "C", "G"]
BIM_A2 = ["G", "T", "A"]
SNPIDS = ["rs0", "rs1", "rs2"]
OBS_CT = 8


def write_frq_standard(path):
    """plink1.9 .frq: A1 = bim_A1, MAF = p  →  no flip needed."""
    with open(path, "w") as f:
        f.write(" CHR          SNP   A1   A2          MAF  NCHROBS\n")
        for snp, a1, a2, p in zip(SNPIDS, BIM_A1, BIM_A2, P):
            f.write(f"   1  {snp:>12s}    {a1}    {a2}     {p:.5f}        {OBS_CT}\n")


def write_frq_swapped(path):
    """plink1.9 .frq: A1 = bim_A2, MAF = 1-p  →  flip needed for all SNPs."""
    with open(path, "w") as f:
        f.write(" CHR          SNP   A1   A2          MAF  NCHROBS\n")
        for snp, a1, a2, p in zip(SNPIDS, BIM_A2, BIM_A1, P):
            f.write(f"   1  {snp:>12s}    {a1}    {a2}     {1-p:.5f}        {OBS_CT}\n")


def write_afreq_alt_is_a1(path):
    """plink2 .afreq: ALT1 = bim_A1, ALT1_FREQ = p  →  no flip needed."""
    with open(path, "w") as f:
        f.write("#CHROM\tID\tREF\tALT1\tALT1_FREQ\tOBS_CT\n")
        for snp, ref, alt1, p in zip(SNPIDS, BIM_A2, BIM_A1, P):
            f.write(f"1\t{snp}\t{ref}\t{alt1}\t{p:.5f}\t{OBS_CT}\n")


def write_afreq_alt_is_a2(path):
    """plink2 .afreq: ALT1 = bim_A2, ALT1_FREQ = 1-p  →  flip needed for all SNPs."""
    with open(path, "w") as f:
        f.write("#CHROM\tID\tREF\tALT1\tALT1_FREQ\tOBS_CT\n")
        for snp, ref, alt1, p in zip(SNPIDS, BIM_A1, BIM_A2, P):
            f.write(f"1\t{snp}\t{ref}\t{alt1}\t{1-p:.5f}\t{OBS_CT}\n")


# ---------------------------------------------------------------------------
# Output reader
# ---------------------------------------------------------------------------

def read_grm_output(path):
    """Return dict (id1, id2) -> grm float."""
    result = {}
    with open(path) as f:
        header = f.readline()   # IID1  IID2  N_SNPs  GRM
        assert "IID1" in header, f"unexpected header: {header!r}"
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            result[(parts[0], parts[1])] = float(parts[3])
    return result


# ---------------------------------------------------------------------------
# Main test
# ---------------------------------------------------------------------------

def main():
    if not os.path.isfile(BINARY):
        print(f"Binary not found: {BINARY}")
        print("Please build first:  cd grm_pairs && make")
        sys.exit(1)

    with tempfile.TemporaryDirectory() as tmpdir:
        bfile  = os.path.join(tmpdir, "test")
        pairs  = os.path.join(tmpdir, "pairs.txt")

        # Write shared plink files
        write_fam(bfile + ".fam")
        write_bim(bfile + ".bim")
        write_bed(bfile + ".bed")
        write_pairs(pairs)

        # Write the four frequency-file variants into the tmp dir
        frq_std   = os.path.join(tmpdir, "frq_standard.frq")
        frq_swap  = os.path.join(tmpdir, "frq_swapped.frq")
        af_a1     = os.path.join(tmpdir, "afreq_alt_a1.afreq")
        af_a2     = os.path.join(tmpdir, "afreq_alt_a2.afreq")

        write_frq_standard(frq_std)
        write_frq_swapped(frq_swap)
        write_afreq_alt_is_a1(af_a1)
        write_afreq_alt_is_a2(af_a2)

        print("=== test_allele_flip: running grm_pairs with four frequency-file variants ===\n")

        # Run each variant sequentially (grm_pairs auto-detects .afreq before .frq)
        all_results = {}
        for label, freq_path in [
            ("frq_standard", frq_std),
            ("frq_swapped",  frq_swap),
            ("afreq_alt_a1", af_a1),
            ("afreq_alt_a2", af_a2),
        ]:
            ext  = os.path.splitext(freq_path)[1]
            dest = bfile + ext
            shutil.copy(freq_path, dest)

            out  = bfile + f"_{label}.out"
            cmd  = [BINARY, bfile, pairs, out]
            proc = subprocess.run(cmd, capture_output=True, text=True)

            # Remove the auto-detect file before next iteration
            if os.path.exists(dest):
                os.remove(dest)

            if proc.returncode != 0:
                print(f"FAIL [{label}]: grm_pairs exited with code {proc.returncode}")
                print(proc.stderr)
                sys.exit(1)

            grm = read_grm_output(out)
            all_results[label] = grm
            print(f"  [{label}]")
            for pair, val in sorted(grm.items()):
                print(f"    GRM({pair[0]}, {pair[1]}) = {val:.8f}")

        print()

        # ------------------------------------------------------------------
        # Assertion 1: all four variants produce identical GRM values
        # ------------------------------------------------------------------
        print("=== Assertion 1: all four variants produce the same GRM ===")
        reference = all_results["frq_standard"]
        failed = False
        for label in ("frq_swapped", "afreq_alt_a1", "afreq_alt_a2"):
            for pair, expected_val in reference.items():
                got = all_results[label].get(pair)
                if got is None:
                    print(f"  FAIL: {label} missing pair {pair}")
                    failed = True
                    continue
                diff = abs(got - expected_val)
                if diff > TOLERANCE:
                    print(f"  FAIL: {label} GRM({pair}) = {got:.8f}, "
                          f"expected {expected_val:.8f} (diff={diff:.2e})")
                    failed = True
                else:
                    print(f"  OK:   {label} GRM({pair}) matches "
                          f"frq_standard (diff={diff:.2e})")
        if not failed:
            print("  All variants agree.\n")

        # ------------------------------------------------------------------
        # Assertion 2: GRM values match analytically computed reference
        #   GRM(s0,s1) = -29/42  ≈ -0.6904762
        #   GRM(s0,s2) = -7/9   ≈ -0.7777778
        # ------------------------------------------------------------------
        print("=== Assertion 2: GRM values match analytic reference ===")
        for pair, expected_val in EXPECTED.items():
            got  = reference.get(pair)
            if got is None:
                print(f"  FAIL: frq_standard missing pair {pair}")
                failed = True
                continue
            diff = abs(got - expected_val)
            if diff > TOLERANCE:
                print(f"  FAIL: GRM({pair}) = {got:.8f}, "
                      f"expected {expected_val:.8f} (diff={diff:.2e})")
                failed = True
            else:
                print(f"  OK:   GRM({pair}) = {got:.8f} "
                      f"≈ {expected_val:.8f} (diff={diff:.2e})")

        print()
        if failed:
            print("RESULT: FAILED")
            sys.exit(1)
        else:
            print("RESULT: ALL TESTS PASSED")


if __name__ == "__main__":
    main()
