#!/usr/bin/env bash
# Run grm_pairs on one pair chunk. This is the unit of work dispatched to a
# single batch job.
#
# Usage: run_chunk.sh <bfile_prefix> <pair_chunk_file> <out_file> [grm_pairs_bin]
set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <bfile_prefix> <pair_chunk_file> <out_file> [grm_pairs_bin]" >&2
    exit 1
fi

BFILE="$1"
PAIR_CHUNK="$2"
OUT_FILE="$3"
GRM_PAIRS_BIN="${4:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../grm_pairs" && pwd)/grm_pairs}"

"$GRM_PAIRS_BIN" "$BFILE" "$PAIR_CHUNK" "$OUT_FILE"
