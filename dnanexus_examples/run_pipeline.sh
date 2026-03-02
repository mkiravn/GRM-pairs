#!/usr/bin/env bash
# run_pipeline.sh -- Submit the full GRM-pairs analysis pipeline on DNAnexus.
#
# Usage:
#   bash run_pipeline.sh [OPTIONS]
#
# Options:
#   --project   <name>   DNAnexus project name (default: current project)
#   --dest      <path>   Output folder in the project  (default: /results)
#   --bfile     <name>   PLINK bfile prefix, e.g. my_data (default: my_data)
#   --pairs     <path>   DNAnexus path to the pair-list file (default: /data/pairs.txt)
#   --pheno     <path>   DNAnexus path to the phenotype file (default: /data/pheno/phenotypes.txt)
#   --covar     <path>   DNAnexus path to the covariate file (default: /data/pheno/covariates.txt)
#   --covar-idx <str>    Covariate column indices, 1-based  (default: "1,2")
#   --intervals <str>    HE-regression GRM intervals        (default: "0-0.05,0.05-0.35,0.35-0.75")
#   --bin-width <float>  grm_pheno_join bin width           (default: 0.001)
#   --seed      <int>    Bootstrap random seed              (default: 42)
#   --instance  <type>   DNAnexus instance type             (default: mem1_ssd1_v2_x4)
#
# Prerequisites:
#   - dx CLI installed and authenticated  (dx login / dx select <project>)
#   - Compiled binaries already uploaded to /tools/ in the project
#     (see README.md for build + upload instructions)

set -euo pipefail

# ── Defaults ─────────────────────────────────────────────────────────────────
DEST="/results"
BFILE="my_data"
PAIRS_DX="/data/pairs.txt"
PHENO_DX="/data/pheno/phenotypes.txt"
COVAR_DX="/data/pheno/covariates.txt"
COVAR_IDX="1,2"
INTERVALS="0-0.05,0.05-0.35,0.35-0.75"
BIN_WIDTH="0.001"
SEED="42"
INSTANCE="mem1_ssd1_v2_x4"

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dest)      DEST="$2";       shift 2 ;;
        --bfile)     BFILE="$2";      shift 2 ;;
        --pairs)     PAIRS_DX="$2";   shift 2 ;;
        --pheno)     PHENO_DX="$2";   shift 2 ;;
        --covar)     COVAR_DX="$2";   shift 2 ;;
        --covar-idx) COVAR_IDX="$2";  shift 2 ;;
        --intervals) INTERVALS="$2";  shift 2 ;;
        --bin-width) BIN_WIDTH="$2";  shift 2 ;;
        --seed)      SEED="$2";       shift 2 ;;
        --instance)  INSTANCE="$2";   shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

PLINK_DIR="/data/plink"

echo "=== GRM-pairs DNAnexus pipeline ==="
echo "  PLINK bfile : ${PLINK_DIR}/${BFILE}.{bed,bim,fam,afreq}"
echo "  Pair list   : ${PAIRS_DX}"
echo "  Phenotypes  : ${PHENO_DX}"
echo "  Covariates  : ${COVAR_DX} (indices: ${COVAR_IDX})"
echo "  HE intervals: ${INTERVALS}"
echo "  Bin width   : ${BIN_WIDTH}"
echo "  Output dir  : ${DEST}"
echo "  Instance    : ${INSTANCE}"
echo ""

# ── Step 1: Residualise phenotypes ───────────────────────────────────────────
echo "[1/5] Submitting regress_y..."
JOB1=$(dx run swiss-army-knife \
    -iin="/tools/regress_y" \
    -iin="${PHENO_DX}" \
    -iin="${COVAR_DX}" \
    -icmd="
        chmod +x regress_y
        ./regress_y \
            \$(basename '${PHENO_DX}') \
            \$(basename '${COVAR_DX}') \
            '${COVAR_IDX}' \
            residualized.txt
    " \
    --instance-type "${INSTANCE}" \
    --destination "${DEST}" \
    --name "1_regress_y" \
    -y --brief)
echo "  Job ID: ${JOB1}"

# ── Step 2: Compute pairwise GRM ─────────────────────────────────────────────
echo "[2/5] Submitting grm_pairs..."
JOB2=$(dx run swiss-army-knife \
    -iin="/tools/grm_pairs" \
    -iin="${PLINK_DIR}/${BFILE}.bed" \
    -iin="${PLINK_DIR}/${BFILE}.bim" \
    -iin="${PLINK_DIR}/${BFILE}.fam" \
    -iin="${PLINK_DIR}/${BFILE}.afreq" \
    -iin="${PAIRS_DX}" \
    -icmd="
        chmod +x grm_pairs
        ./grm_pairs '${BFILE}' \$(basename '${PAIRS_DX}') grm_output.txt
    " \
    --instance-type "${INSTANCE}" \
    --destination "${DEST}" \
    --name "2_grm_pairs" \
    -y --brief)
echo "  Job ID: ${JOB2}"

# ── Step 3: Compute phenotype cross-products ──────────────────────────────────
echo "[3/5] Submitting pheno_pairs (depends on step 1)..."
JOB3=$(dx run swiss-army-knife \
    -iin="/tools/pheno_pairs" \
    -iin="${PAIRS_DX}" \
    -iin="${DEST}/residualized.txt" \
    -icmd="
        chmod +x pheno_pairs
        ./pheno_pairs \$(basename '${PAIRS_DX}') residualized.txt pheno_cross.txt
    " \
    --depends-on "${JOB1}" \
    --instance-type "${INSTANCE}" \
    --destination "${DEST}" \
    --name "3_pheno_pairs" \
    -y --brief)
echo "  Job ID: ${JOB3}"

# ── Step 4: Join GRM + phenotype cross-products and bin ───────────────────────
echo "[4/5] Submitting grm_pheno_join (depends on steps 2 and 3)..."
JOB4=$(dx run swiss-army-knife \
    -iin="/tools/grm_pheno_join" \
    -iin="${DEST}/grm_output.txt" \
    -iin="${DEST}/pheno_cross.txt" \
    -icmd="
        chmod +x grm_pheno_join
        ./grm_pheno_join grm_output.txt pheno_cross.txt joined.txt '${BIN_WIDTH}'
    " \
    --depends-on "${JOB2}" \
    --depends-on "${JOB3}" \
    --instance-type "${INSTANCE}" \
    --destination "${DEST}" \
    --name "4_grm_pheno_join" \
    -y --brief)
echo "  Job ID: ${JOB4}"

# ── Step 5: HE regression ────────────────────────────────────────────────────
echo "[5/5] Submitting he_regression (depends on step 4)..."
JOB5=$(dx run swiss-army-knife \
    -iin="/tools/he_regression" \
    -iin="${DEST}/joined.txt" \
    -icmd="
        chmod +x he_regression
        ./he_regression \
            -f joined.txt \
            -i '${INTERVALS}' \
            -b 1000 \
            -s '${SEED}' \
            -o he_results.txt
    " \
    --depends-on "${JOB4}" \
    --instance-type "${INSTANCE}" \
    --destination "${DEST}" \
    --name "5_he_regression" \
    -y --brief)
echo "  Job ID: ${JOB5}"

echo ""
echo "All jobs submitted successfully."
echo ""
echo "Monitor progress:"
echo "  dx watch ${JOB5}"
echo "  dx find jobs --state '*' --brief"
echo ""
echo "Download results when complete:"
echo "  mkdir -p ./grm_results"
echo "  dx download ${DEST}/grm_output.txt  -o ./grm_results/"
echo "  dx download ${DEST}/pheno_cross.txt -o ./grm_results/"
echo "  dx download ${DEST}/joined.txt      -o ./grm_results/"
echo "  dx download ${DEST}/he_results.txt  -o ./grm_results/"
