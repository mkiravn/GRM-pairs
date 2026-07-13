#!/usr/bin/env bash
# dsub_submit_example.sh -- submit one grm_pairs chunk job to Google Cloud
# Batch from a laptop, analogous to ../../dnanexus_examples/run_pipeline.sh
# but for the All of Us / Google Cloud setup instead of DNAnexus.
#
# This is illustrative, not wired up: bucket paths, project, and region are
# placeholders. Requires `dsub` (https://github.com/DataBiosphere/dsub) and
# `gcloud` auth against your AoU workspace project.
#
# Usage: bash dsub_submit_example.sh
set -euo pipefail

PROJECT="my-aou-workspace-project"
REGION="us-central1"
BUCKET="gs://my-aou-workspace-bucket"

TOOLS_DIR="${BUCKET}/tools"
DATA_DIR="${BUCKET}/data/plink"     # expects my_data.bed/.bim/.fam/.afreq here
CHUNK_DIR="${BUCKET}/data/pair_chunks"
RESULTS_DIR="${BUCKET}/results"
LOG_DIR="${BUCKET}/logs"

# ── Upload compiled binary and inputs (one-time) ────────────────────────────
# gsutil cp ../../grm_pairs/grm_pairs "${TOOLS_DIR}/"
# gsutil cp my_data.bed my_data.bim my_data.fam my_data.afreq "${DATA_DIR}/"
# python3 ../chunking/split_pairs.py pairs.txt 20 pairs_chunk
# gsutil cp pairs_chunk_*.txt "${CHUNK_DIR}/"

# ── Submit one job per pair chunk ────────────────────────────────────────────
# For a handful of chunks a simple loop is fine. For many chunks, dsub's
# --tasks <TSV> flag (one row per chunk) submits them as a single task array
# instead -- worth switching to once the chunk count grows.
for CHUNK in pairs_chunk_0001 pairs_chunk_0002; do
    dsub \
        --provider google-batch \
        --project "${PROJECT}" \
        --regions "${REGION}" \
        --logging "${LOG_DIR}" \
        --machine-type n2-standard-4 \
        --disk-size 50 \
        --input-recursive BFILE_DIR="${DATA_DIR}" \
        --input PAIR_CHUNK="${CHUNK_DIR}/${CHUNK}.txt" \
        --input GRM_PAIRS_BIN="${TOOLS_DIR}/grm_pairs" \
        --output OUT_FILE="${RESULTS_DIR}/${CHUNK}.grm.txt" \
        --command '
            chmod +x "${GRM_PAIRS_BIN}"
            BFILE_PREFIX="${BFILE_DIR}/my_data"
            "${GRM_PAIRS_BIN}" "${BFILE_PREFIX}" "${PAIR_CHUNK}" "${OUT_FILE}"
        ' \
        --name "grm_pairs_${CHUNK}" \
        --wait=false
done

echo "Submitted jobs. Check status with: dstat --provider google-batch --project ${PROJECT} --status '*'"
echo "Once all chunks finish, merge with:"
echo "  gsutil cp ${RESULTS_DIR}/pairs_chunk_*.grm.txt ."
echo "  bash ../plink/merge_chunks.sh merged.grm.txt pairs_chunk_*.grm.txt"
