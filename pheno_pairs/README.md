# Pheno Pairs

A C program to calculate cross products of phenotypic values for pairs of individuals.

## Compilation

```bash
make
```

or

```bash
gcc -Wall -Wextra -O2 -std=c99 -o pheno_pairs main.c
```

## Usage

```bash
./pheno_pairs <pair_file> <pheno_file> <out_file>
```

### Input Files

1. **Pair file**: Tab or space-delimited file with pairs of individual IDs
   ```
   IID1    IID2
   sample1 sample2
   sample1 sample3
   sample2 sample4
   ```

2. **Phenotype file**: Tab or space-delimited file with individual IDs and multiple phenotype values
   - First column: Individual IDs (any header name)  
   - Remaining columns: Phenotype values (any header names)
   ```
   IID     height  weight  age
   sample1 1.75    70.5    25
   sample2 1.82    85.0    30
   sample3 1.68    60.0    22
   ```
   
   Or even:
   ```
   ID      bmi     score   valid
   bob     23.5    85.2    1
   alice   21.3    78.0    0
   charlie 25.1    92.5    1
   ```

### Output

Tab-delimited file with cross products for each phenotype column:
```
IID1    IID2    crossproduct_height crossproduct_weight crossproduct_age
sample1 sample2 3.185              5992.5             750
sample1 sample3 2.94               4230               550
sample2 sample3 3.0576             5100               660
```

## Features

- **Multi-phenotype support**: Handles multiple phenotype columns and calculates cross products for each
- **Dynamic column naming**: Uses column headers to name output columns (e.g., "crossproduct_height")
- **Simple and robust**: Assumes first column is always IID, remaining columns are phenotype values
- **Flexible input**: Works with any column headers or no headers at all
- **Missing data handling**: Outputs "NA" for cross products when phenotype values are missing or non-numeric
- **Reports summary statistics**: Shows number of pairs processed and success/failure counts
- **Memory-efficient**: Processes large files efficiently

## Example

```bash
# Create test files
echo -e "IID1\tIID2\nsample1\tsample2\nsample1\tsample3" > pairs.txt
echo -e "IID\tpheno\nsample1\t1.5\nsample2\t-2.0\nsample3\t0.5" > pheno.txt

# Run the program
./pheno_pairs pairs.txt pheno.txt output.txt

# Check results
cat output.txt
```

This will produce:
```
IID1    IID2    crossproduct_measurement
sample1 sample2 -3
sample1 sample3 0.75
```
