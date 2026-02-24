# Minimal PLINK Dummy Example

This directory contains a tiny PLINK dataset (3 samples, 2 SNPs) that
demonstrates how to run `grm_pairs` with both plink1.9 `.frq` and plink2
`.afreq` frequency files.

## Files

| File          | Description                                      |
|---------------|--------------------------------------------------|
| `dummy.fam`   | Sample information (3 individuals)               |
| `dummy.bim`   | Variant information (2 SNPs: rs1 A/G, rs2 C/T)  |
| `dummy.bed`   | Genotypes in PLINK binary format                 |
| `dummy.frq`   | Allele frequencies – plink1.9 `.frq` format      |
| `dummy.afreq` | Allele frequencies – plink2 `.afreq` format      |
| `pairs.txt`   | Pairs of individuals to compute GRM for          |

## Genotype layout

```
        rs1  rs2
ind1:    0    1      (0=HOM_A1, 1=HET, 2=HOM_A2)
ind2:    1    0
ind3:    2    1
```

## Running the example

First build the binary (from the repo root):

```bash
cd grm_pairs && make
```

### With plink1.9 `.frq`

```bash
# Copy only the .frq file (so auto-detection picks it up)
cp examples/dummy.frq /tmp/dummy.frq
cp examples/dummy.{fam,bim,bed} /tmp/
cp examples/pairs.txt /tmp/

grm_pairs/grm_pairs /tmp/dummy /tmp/pairs.txt /tmp/out_frq.txt
cat /tmp/out_frq.txt
```

### With plink2 `.afreq`

```bash
# .afreq takes precedence over .frq when both are present
cp examples/dummy.afreq /tmp/dummy.afreq
cp examples/dummy.{fam,bim,bed} /tmp/
cp examples/pairs.txt /tmp/

grm_pairs/grm_pairs /tmp/dummy /tmp/pairs.txt /tmp/out_afreq.txt
cat /tmp/out_afreq.txt
```

Both commands produce the same output (GRM values are printed with up to 8
significant figures via `%.8g`):

```
IID1    IID2    N_SNPs  GRM
ind1    ind2    2       0.500015
ind1    ind3    2       -0.87499437
```

## Notes

- `grm_pairs` auto-detects the frequency format: it tries `<bfile>.afreq`
  first (plink2); if not found it falls back to `<bfile>.frq` (plink1.9).
- The allele reported in the frequency file is compared against the A1 column
  in the `.bim`. If they differ, the frequency is flipped (`1-p`) so that `p`
  always equals `freq(bim_A1)`.
- The verbose/debug build (`grm_pairs_debug/`) produces identical numerical
  results and adds per-SNP debug output to stderr.
