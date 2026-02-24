# grm_pairs

Computes the genetic relationship matrix (GRM) for a user-specified list of
individual pairs, without building the full N×N matrix.

This is useful when you only need relatedness estimates for a sparse set of
pairs — for example, known relatives, or pairs identified by a prior IBD scan
— and N is too large to compute or store the full GRM.

## Method

Implements the Yang et al. 2010 (GCTA) variance-standardised estimator:

```
A_ij = (1/M) * sum_k [ (x_ik - 2*p_k)(x_jk - 2*p_k) / (2*p_k*(1-p_k)) ]
```

where `x_ik` is the count of A1 alleles (0, 1, 2) for individual `i` at
SNP `k`, `p_k` is the A1 allele frequency, and M is the number of SNPs
where neither individual is missing.

Allele frequencies are read from a plink frequency file and checked for
consistent allele coding against the `.bim` file (see **Allele-flip check**
below).

## Building

Requires gcc and make. No external libraries needed beyond the standard C
math library.

```bash
git clone <this-repo>
cd grm_pairs
make
```

On a DNAnexus UKB node (Ubuntu 24.04 x86_64) or any standard Linux cluster:

```bash
make
# optional: put the binary on your PATH
make install   # copies to ~/bin/
```

## Usage

```
grm_pairs <bfile> <pair_file> <out_file>
```

**Arguments:**

| Argument | Description |
|---|---|
| `bfile` | Plink binary fileset prefix — expects `<bfile>.bed`, `<bfile>.bim`, `<bfile>.fam` |
| `pair_file` | Whitespace-delimited file of IID pairs, one pair per line |
| `out_file` | Output path |

**Pair file format:**
```
IID1  IID2
per0  per1
per0  per5
per3  per9
```

**Output format** (tab-delimited):
```
IID1    IID2    N_SNPs  GRM
per0    per1    500000  0.0023451
per0    per5    500000  -0.0041223
per3    per9    500000  0.4981234
```

`N_SNPs` is the number of non-missing SNPs for that pair. `GRM` matches
the off-diagonal entries of `plink --make-rel`.

## Frequency files

The program reads allele frequencies from either a **plink2** or **plink1.9**
frequency file.  It tries `<bfile>.afreq` first; if that file does not exist
it falls back to `<bfile>.frq`.  The format is auto-detected from the header
line.

| Format | Produced by | Header | Columns |
|---|---|---|---|
| `.afreq` | `plink2 --freq` | starts with `#CHROM` | `#CHROM ID REF ALT1 ALT1_FREQ OBS_CT` |
| `.frq` | `plink1.9 --freq` | starts with `CHR` | `CHR SNP A1 A2 MAF NCHROBS` |

### Allele-flip check

The `.bed` file encodes genotypes as the number of copies of **bim A1**
(0 = hom-A1, 1 = het, 2 = hom-A2).  The GRM formula therefore requires
`p_k = freq(bim_A1)`.

Because `.afreq` files report `ALT1_FREQ` (which may correspond to either
bim A1 or bim A2 depending on how the reference was assigned), and because
`.frq` files can in principle list a different allele as A1, the program
compares the reported allele against the A1 column of the `.bim` file for
every SNP.  If they differ, the frequency is flipped (`p → 1 − p`) before
being passed to the GRM calculation.

## Example

```bash
# generate test data with plink
plink --dummy 200 5000 --make-bed --out test_data

# create a pair file
echo "per0 per1
per0 per2
per1 per2" > my_pairs.txt

# compute GRM for those pairs using a plink2 frequency file
plink2 --bfile test_data --freq --out test_data
./grm_pairs test_data my_pairs.txt output.txt

# or using a plink1.9 frequency file
plink --bfile test_data --freq --out test_data
./grm_pairs test_data my_pairs.txt output.txt

# validate against plink ground truth
plink --bfile test_data --make-rel square --out truth
```

## Running the test

A self-contained Python test verifies allele-flip correctness using a small
synthetic dataset:

```bash
python3 test/test_allele_flip.py
```

The test creates a 4-sample / 3-SNP dataset and runs `grm_pairs` four times —
once with each combination of `.frq`/`.afreq` format and allele orientation
(no flip needed / flip needed for all SNPs) — and asserts that all four runs
produce identical GRM values matching an analytical reference.
