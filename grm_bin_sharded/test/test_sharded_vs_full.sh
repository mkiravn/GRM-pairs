#!/usr/bin/env bash
# End-to-end local test: does accumulate+merge over shards reproduce
# blocked_grm_bin_means's single-pass output exactly, on real 1000 Genomes
# genotypes and a simulated additive phenotype?
#
# No working plink/gcta on this machine, so a dense GRM is manufactured from
# the already-validated grm_pairs sparse pairwise tool (see build_fake_grm.py)
# and the phenotype is simulated directly from the .bed file (see
# simulate_pheno.py, standing in for `plink2 --score`).
#
# Usage: bash test_sharded_vs_full.sh [n_individuals] [n_shards] [nblocks]
set -euo pipefail

N_INDIV="${1:-80}"
N_SHARDS="${2:-4}"
NBLOCKS="${3:-10}"
SEED=1
H2=0.5
N_SNPS=5000

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARDED_DIR="$(dirname "${SCRIPT_DIR}")"
REPO_DIR="$(dirname "${SHARDED_DIR}")"

BFILE="${REPO_DIR}/1kg/1000g_chr1"
GRM_PAIRS_BIN="${REPO_DIR}/grm_pairs/grm_pairs"
BLOCKED_BIN="${REPO_DIR}/full_grm_bin/blocked_grm_bin_means"
SHARD_TOOL="${SHARDED_DIR}/grm_shard_tool"
BINS_FILE="${REPO_DIR}/full_grm_bin/bins.txt"

WORK_DIR="$(mktemp -d)"
trap 'rm -rf "${WORK_DIR}"' EXIT
echo "[INFO] Work dir: ${WORK_DIR}"

for bin in "${GRM_PAIRS_BIN}" "${BLOCKED_BIN}" "${SHARD_TOOL}"; do
    if [ ! -x "${bin}" ]; then
        echo "[ERROR] Missing binary: ${bin} (build it first)" >&2
        exit 1
    fi
done

# --- pick a subset of individuals, in a fixed order ---------------------
IDS_FILE="${WORK_DIR}/ids.txt"
awk -v n="${N_INDIV}" 'NR <= n {print $2}' "${BFILE}.fam" > "${IDS_FILE}"
echo "[INFO] Using ${N_INDIV} individuals"

# --- all-pairs list for grm_pairs ----------------------------------------
PAIRS_FILE="${WORK_DIR}/pairs.txt"
awk '
    { ids[NR] = $1 }
    END {
        for (i = 1; i <= NR; i++)
            for (j = i + 1; j <= NR; j++)
                print ids[i], ids[j]
    }
' "${IDS_FILE}" > "${PAIRS_FILE}"

# --- manufacture a dense GRM ----------------------------------------------
GRM_PAIRS_OUT="${WORK_DIR}/grm_pairs_out.tsv"
"${GRM_PAIRS_BIN}" "${BFILE}" "${PAIRS_FILE}" "${GRM_PAIRS_OUT}"

GRM_PREFIX="${WORK_DIR}/fake"
python3 "${SCRIPT_DIR}/build_fake_grm.py" "${GRM_PAIRS_OUT}" "${IDS_FILE}" "${BFILE}" "${GRM_PREFIX}"

# --- simulate an additive phenotype ---------------------------------------
PHENO_FILE="${WORK_DIR}/pheno.txt"
python3 "${SCRIPT_DIR}/simulate_pheno.py" "${BFILE}" "${IDS_FILE}" "${N_SNPS}" "${H2}" "${SEED}" "${PHENO_FILE}"

# --- reference: single-pass blocked_grm_bin_means -------------------------
REF_PREFIX="${WORK_DIR}/reference"
"${BLOCKED_BIN}" \
    --grm-prefix "${GRM_PREFIX}" \
    --pheno "${PHENO_FILE}" \
    --bins "${BINS_FILE}" \
    --out-prefix "${REF_PREFIX}" \
    --nblocks "${NBLOCKS}" \
    --seed "${SEED}"

# --- sharded path: split, accumulate per shard, merge ---------------------
SHARD_PREFIX="${WORK_DIR}/shard"
python3 "${SCRIPT_DIR}/split_grm_shards.py" "${GRM_PREFIX}.grm.bin" "${N_INDIV}" "${N_SHARDS}" "${SHARD_PREFIX}"

ACC_LIST="${WORK_DIR}/acc_list.txt"
> "${ACC_LIST}"
for k in $(seq 1 "${N_SHARDS}"); do
    ACC_OUT="${SHARD_PREFIX}.${k}.acc.tsv"
    "${SHARD_TOOL}" accumulate \
        --grm-id "${GRM_PREFIX}.grm.id" \
        --shard "${SHARD_PREFIX}.${k}" \
        --parallel "${k}" "${N_SHARDS}" \
        --pheno "${PHENO_FILE}" \
        --bins "${BINS_FILE}" \
        --nblocks "${NBLOCKS}" --seed "${SEED}" \
        --out "${ACC_OUT}"
    echo "${ACC_OUT}" >> "${ACC_LIST}"
done

SHARDED_PREFIX="${WORK_DIR}/sharded"
"${SHARD_TOOL}" merge \
    --acc-list "${ACC_LIST}" \
    --bins "${BINS_FILE}" \
    --nblocks "${NBLOCKS}" \
    --out-prefix "${SHARDED_PREFIX}"

# --- non-vacuous check: make sure real pairs actually got counted, not just
# a silently-empty comparison (e.g. from a FID/IID mismatch dropping every
# phenotype) ------------------------------------------------------------------
TOTAL_N=$(awk -F'\t' 'NR>1{sum+=$7} END{print sum+0}' "${REF_PREFIX}.full.tsv")
EXPECTED_PAIRS=$(wc -l < "${PAIRS_FILE}")
echo "[INFO] Reference full.tsv covers ${TOTAL_N} pairs (of ${EXPECTED_PAIRS} generated)"
if [ "${TOTAL_N}" -eq 0 ]; then
    echo "[FAIL] Reference output has zero pairs in every bin -- phenotype alignment is broken" >&2
    exit 1
fi

# --- compare ----------------------------------------------------------------
if diff -q "${REF_PREFIX}.full.tsv" "${SHARDED_PREFIX}.full.tsv" > /dev/null && \
   diff -q "${REF_PREFIX}.jk.tsv" "${SHARDED_PREFIX}.jk.tsv" > /dev/null; then
    echo "[PASS] sharded accumulate+merge output matches single-pass blocked_grm_bin_means exactly"
else
    echo "[FAIL] outputs differ"
    echo "--- full.tsv diff ---"
    diff "${REF_PREFIX}.full.tsv" "${SHARDED_PREFIX}.full.tsv" | head -n 20
    echo "--- jk.tsv diff ---"
    diff "${REF_PREFIX}.jk.tsv" "${SHARDED_PREFIX}.jk.tsv" | head -n 20
    exit 1
fi
