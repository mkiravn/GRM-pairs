#!/usr/bin/env bash
# Preliminary local test: does splitting a pair list into chunks and running
# grm_pairs per chunk (as separate "jobs" would run on AoU) give the same
# result as one grm_pairs run over the whole pair list?
#
# Uses the 1000 Genomes chr1 data already checked into ../../1kg/, so this
# runs entirely locally with no AoU access.
#
# Usage: bash test_chunking_consistency.sh [n_chunks]
set -euo pipefail

N_CHUNKS="${1:-4}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AOU_DIR="$(dirname "${SCRIPT_DIR}")"
REPO_DIR="$(dirname "${AOU_DIR}")"

BFILE="${REPO_DIR}/1kg/1000g_chr1"
GRM_PAIRS_BIN="${REPO_DIR}/grm_pairs/grm_pairs"

WORK_DIR="$(mktemp -d)"
trap 'rm -rf "${WORK_DIR}"' EXIT

echo "[INFO] Work dir: ${WORK_DIR}"

if [ ! -x "${GRM_PAIRS_BIN}" ]; then
    echo "[INFO] Building grm_pairs binary"
    (cd "${REPO_DIR}/grm_pairs" && make)
fi

# --- build a modest all-vs-all pair list from the first 60 samples ----------
PAIRS_FILE="${WORK_DIR}/pairs.txt"
awk 'NR <= 60 {print $2}' "${BFILE}.fam" > "${WORK_DIR}/ids.txt"
awk '
    { ids[NR] = $1 }
    END {
        for (i = 1; i <= NR; i++)
            for (j = i + 1; j <= NR; j++)
                print ids[i], ids[j]
    }
' "${WORK_DIR}/ids.txt" > "${PAIRS_FILE}"
echo "[INFO] Built $(wc -l < "${PAIRS_FILE}") pairs from 60 samples"

# --- single-shot run ---------------------------------------------------------
SINGLE_OUT="${WORK_DIR}/single_shot.txt"
"${GRM_PAIRS_BIN}" "${BFILE}" "${PAIRS_FILE}" "${SINGLE_OUT}"

# --- chunked run --------------------------------------------------------------
python3 "${AOU_DIR}/chunking/split_pairs.py" "${PAIRS_FILE}" "${N_CHUNKS}" "${WORK_DIR}/chunk"

CHUNK_OUTPUTS=()
for CHUNK_FILE in "${WORK_DIR}"/chunk_*.txt; do
    CHUNK_OUT="${CHUNK_FILE%.txt}.grm.txt"
    bash "${AOU_DIR}/plink/run_chunk.sh" "${BFILE}" "${CHUNK_FILE}" "${CHUNK_OUT}" "${GRM_PAIRS_BIN}"
    CHUNK_OUTPUTS+=("${CHUNK_OUT}")
done

MERGED_OUT="${WORK_DIR}/merged.txt"
bash "${AOU_DIR}/plink/merge_chunks.sh" "${MERGED_OUT}" "${CHUNK_OUTPUTS[@]}"

# --- compare ------------------------------------------------------------------
SINGLE_SORTED="${WORK_DIR}/single_sorted.txt"
MERGED_SORTED="${WORK_DIR}/merged_sorted.txt"
tail -n +2 "${SINGLE_OUT}" | sort > "${SINGLE_SORTED}"
tail -n +2 "${MERGED_OUT}" | sort > "${MERGED_SORTED}"

if diff -q "${SINGLE_SORTED}" "${MERGED_SORTED}" > /dev/null; then
    echo "[PASS] chunked+merged output matches single-shot output ($(wc -l < "${SINGLE_SORTED}") rows)"
else
    echo "[FAIL] chunked+merged output differs from single-shot output"
    diff "${SINGLE_SORTED}" "${MERGED_SORTED}" | head -n 20
    exit 1
fi
