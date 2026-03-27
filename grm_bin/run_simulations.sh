#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: run_simulations.sh <input_bfile> <out_dir> [n_individuals] [n_snps] [n_reps] [h2] [seed]

Environment variables:
  PLINK     Path to plink/plink2 command used for subsetting and QC (default: plink)
  GCTA      Path to gcta64 used to write .grm.bin/.grm.id (default: gcta64)
  MAF       Minor-allele frequency threshold for QC subset (default: 0.05)
  GENO      Missingness threshold for QC subset (default: 0.01)
  JK_BLOCKS Number of blocks for grm_bin_cov_jackknife (default: 50)

Example:
  ./run_simulations.sh /data/1000g/chr20 ./sim_1kg 800 20000 20 0.4 123
EOF
  exit 1
fi

INPUT_BFILE="$1"
OUT_DIR="$2"
N_INDIVIDUALS="${3:-800}"
N_SNPS="${4:-20000}"
N_REPS="${5:-20}"
H2="${6:-0.4}"
SEED="${7:-1}"

PLINK="${PLINK:-plink}"
GCTA="${GCTA:-gcta64}"
MAF="${MAF:-0.05}"
GENO="${GENO:-0.01}"
JK_BLOCKS="${JK_BLOCKS:-50}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "${OUT_DIR}"

echo "Building grm_bin binaries"
make -C "${SCRIPT_DIR}"

echo "Sampling individuals from ${INPUT_BFILE}.fam"
awk '{print $1, $2}' "${INPUT_BFILE}.fam" | shuf --random-source=<(yes "${SEED}") | head -n "${N_INDIVIDUALS}" > "${OUT_DIR}/keep.txt"

echo "PLINK QC + LD pruning"
"${PLINK}" \
  --bfile "${INPUT_BFILE}" \
  --keep "${OUT_DIR}/keep.txt" \
  --maf "${MAF}" \
  --geno "${GENO}" \
  --indep-pairwise 200 50 0.2 \
  --out "${OUT_DIR}/prune"

"${PLINK}" \
  --bfile "${INPUT_BFILE}" \
  --keep "${OUT_DIR}/keep.txt" \
  --extract "${OUT_DIR}/prune.prune.in" \
  --make-bed \
  --out "${OUT_DIR}/subset_ld"

if [[ "${N_SNPS}" -gt 0 ]]; then
  echo "Limiting to ${N_SNPS} SNPs for simulation speed"
  awk 'NR<='"${N_SNPS}"' {print $2}' "${OUT_DIR}/subset_ld.bim" > "${OUT_DIR}/extract.snps"
  "${PLINK}" \
    --bfile "${OUT_DIR}/subset_ld" \
    --extract "${OUT_DIR}/extract.snps" \
    --make-bed \
    --out "${OUT_DIR}/simulation_panel"
else
  cp "${OUT_DIR}/subset_ld.bed" "${OUT_DIR}/simulation_panel.bed"
  cp "${OUT_DIR}/subset_ld.bim" "${OUT_DIR}/simulation_panel.bim"
  cp "${OUT_DIR}/subset_ld.fam" "${OUT_DIR}/simulation_panel.fam"
fi

echo "Building GRM in .grm.bin format"
"${GCTA}" \
  --bfile "${OUT_DIR}/simulation_panel" \
  --make-grm-bin \
  --out "${OUT_DIR}/simulation_panel"

echo "Running ${N_REPS} phenotype replicates"
for rep in $(seq 1 "${N_REPS}"); do
  rep_dir="${OUT_DIR}/rep_${rep}"
  rep_seed=$((SEED + rep - 1))
  mkdir -p "${rep_dir}"

  Rscript "${SCRIPT_DIR}/simulate_pheno_from_grm.R" \
    "${OUT_DIR}/simulation_panel" \
    "${rep_dir}/pheno.txt" \
    --h2 "${H2}" \
    --seed "${rep_seed}" \
    --replicate "${rep}" \
    --truth-out "${rep_dir}/truth.tsv" \
    --meta-out "${rep_dir}/meta.tsv"

  "${SCRIPT_DIR}/grm_bin_cov" \
    "${rep_dir}/pheno.txt" \
    "${OUT_DIR}/simulation_panel" \
    "${rep_dir}/binned.tsv"

  "${SCRIPT_DIR}/grm_bin_cov_jackknife" \
    "${rep_dir}/pheno.txt" \
    "${OUT_DIR}/simulation_panel" \
    "${rep_dir}/jackknife.tsv" \
    "${JK_BLOCKS}"
done

Rscript "${SCRIPT_DIR}/summarize_simulations.R" \
  "${OUT_DIR}" \
  "${OUT_DIR}/simulation_summary.tsv"

cat <<EOF
Finished.
Panel prefix: ${OUT_DIR}/simulation_panel
Summary:      ${OUT_DIR}/simulation_summary.tsv
EOF
