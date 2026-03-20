# phenoCovBins - GRM Bin Analysis

A C implementation for analyzing phenotype-GRM covariance with binning, equivalent to the Fortran `phenoCovBins` program. Includes a preprocessing step to sort phenotype files.

## Compilation

```bash
make
```

This creates two executables:
- `sort_pheno` - Preprocesses and sorts phenotype files
- `pheno_cov_bins` - Main analysis program

## Workflow

### Step 1: Sort Phenotype File

Before running the analysis, sort your phenotype file to match the GRM individual ordering:

```bash
./sort_pheno fileNameGrmPrefix fileNamePhenoInput fileNamePhenoOutput
```

#### Arguments

- `fileNameGrmPrefix`: Prefix for GRM files (program expects `fileNameGrmPrefix.grm.id`)
- `fileNamePhenoInput`: Input phenotype file
- `fileNamePhenoOutput`: Output sorted phenotype file

#### Input Format

The input phenotype file should be two-column format with individual ID and phenotype value:
```
individual_id phenotype_value
```

Missing values can be coded as:
- `NA`
- Empty values
- Any non-numeric string

#### Output Features

The output file will have:
- Same individual ID format as the GRM ID file
- Phenotype values in GRM individual order
- Missing values coded as `-999` for:
  - Individuals in GRM but not in input phenotype file
  - Individuals with non-numeric phenotype values in input
  - Individuals with negative IID values in GRM

### Step 2: Run GRM Bin Analysis

After preparing the phenotype file, run the main analysis:

```bash
./pheno_cov_bins nIndiv fileNamePhenoPrefix fileNameGrmPrefix
```

#### Arguments

- `nIndiv`: Number of individuals in the dataset
- `fileNamePhenoPrefix`: Prefix for the **sorted** phenotype file (program expects `fileNamePhenoPrefix.txt`)
- `fileNameGrmPrefix`: Prefix for GRM files (program expects `fileNameGrmPrefix.grm.bin` and `fileNameGrmPrefix.grm.id`)

#### Input Files

- **Phenotype file** (`*.txt`): Output from `sort_pheno` step
- **GRM files** (`.grm.bin` and `.grm.id`): PLINK2 format binary GRM files
  - Generated with: `plink2 --bfile grm_input --make-grm-bin --out grm_full`

#### Output

- Text file: `phenCovariance_fileNamePhenoPrefix.txt`
  - Columns: bin_index, lower_bound, upper_bound, avg_grm_value, avg_covariance, count

## Binning Strategy

- **-0.3 to 0.02**: Bins at 0.001 intervals (fine granularity for near-zero relationships)
- **Above 0.02 to 1.05**: Bins at 0.005 intervals (coarser granularity for higher relationships)

## Notes

- After sort_pheno runs, your phenotype file will be in GRM order with proper missing data encoding
- The pheno_cov_bins program expects the phenotype file to already be sorted and formatted correctly
- Phenotype standardization (z-score normalization) is done internally by pheno_cov_bins
- Individuals with negative IID values in the GRM ID file are automatically marked as missing (-999)
  regardless of whether they appear in the input phenotype file
