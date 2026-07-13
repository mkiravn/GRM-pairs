#!/usr/bin/env bash
# Concatenate per-chunk grm_pairs outputs into one file: one header, rows in
# chunk order.
#
# Usage: merge_chunks.sh <out_file> <chunk_output_1> [chunk_output_2 ...]
set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <out_file> <chunk_output_1> [chunk_output_2 ...]" >&2
    exit 1
fi

OUT_FILE="$1"
shift

head -n 1 "$1" > "$OUT_FILE"
for chunk_file in "$@"; do
    tail -n +2 "$chunk_file" >> "$OUT_FILE"
done
