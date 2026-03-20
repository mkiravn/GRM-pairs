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
./pheno_cov_bins nIndiv fileNamePhenoPrefix fileNameGrmPrefix [--id-file ID_FILE] [--bin-file BIN_FILE]
```

#### Arguments

- `nIndiv`: Number of individuals in the dataset/shard
- `fileNamePhenoPrefix`: Prefix for the **sorted** phenotype file (program expects `fileNamePhenoPrefix.txt`)
- `fileNameGrmPrefix`: Prefix for GRM files (program expects `fileNameGrmPrefix.grm.bin` and `fileNameGrmPrefix.grm.id`)
- `--id-file ID_FILE`: Optional separate ID file (for shards with different individual sets)
- `--bin-file BIN_FILE`: Optional separate binary file (for shards)

#### Examples

**Full GRM:**
```bash
./pheno_cov_bins 500000 pheno_sorted grm_full
```

**GRM Shard:**
```bash
./pheno_cov_bins 5000 pheno_sorted --id-file grm_full.grm.id.1 --bin-file grm_full.grm.bin.1
```

## DNAnexus Usage

For running on DNAnexus with GRM shards:

```bash
pheno_file="results/phenotypes/p51_i0/p51_i0_residualized_cov1234.txt"
grm="grm_full.grm.bin.1"
id="grm_full.grm.id"

dx run app-swiss-army-knife \
  -iimage="ubuntu:20.04" \
  -icmd="
    set -euo pipefail
    
    # Get number of individuals in this shard
    N=\$(wc -l < ${id})
    echo \"Number of individuals in shard: \$N\"
    
    # Clone and build tools
    git clone https://github.com/mkiravn/GRM-pairs.git
    cd GRM-pairs/grm_bin && make && cd ../..
    
    # Sort phenotype (using global ID file)
    ./GRM-pairs/grm_bin/sort_pheno grm_full ${pheno_file} pheno_sorted.txt
    
    # Run analysis on shard
    ./GRM-pairs/grm_bin/pheno_cov_bins \$N pheno_sorted --id-file ${id} --bin-file ${grm}
    
    # Clean up
    rm -rf GRM-pairs/
  " \
  --destination "/" \
  --name "grm_bin_shard_analysis" \
  --instance-type "mem1_ssd1_v2_x4" \
  --yes
```

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
