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

Allele frequencies are computed directly from the `.bed` file in a first
pass, guaranteeing consistency with the genotype encoding.

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

## Example

```bash
# generate test data with plink
plink --dummy 200 5000 --make-bed --out test_data

# create a pair file
echo "per0 per1
per0 per2
per1 per2" > my_pairs.txt

# compute GRM for those pairs
./grm_pairs test_data my_pairs.txt output.txt

# validate against plink ground truth
plink --bfile test_data --make-rel square --out truth
```
