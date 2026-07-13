#!/usr/bin/env python3
"""Split a whitespace-delimited pair file into N contiguous chunk files.

Usage:
    split_pairs.py <pair_file> <n_chunks> <out_prefix>

Writes <out_prefix>_0001.txt .. <out_prefix>_NNNN.txt, each with a
contiguous slice of the input lines (blank lines and lines starting with
'#' are dropped before splitting). Chunk sizes differ by at most one line.
"""
import sys


def main() -> None:
    if len(sys.argv) != 4:
        sys.exit(__doc__)

    pair_file, n_chunks_str, out_prefix = sys.argv[1:4]
    n_chunks = int(n_chunks_str)
    if n_chunks < 1:
        sys.exit("n_chunks must be >= 1")

    with open(pair_file) as f:
        lines = [line for line in f if line.strip() and not line.startswith("#")]

    n = len(lines)
    if n == 0:
        sys.exit(f"No pairs found in {pair_file}")
    if n_chunks > n:
        sys.exit(f"n_chunks ({n_chunks}) exceeds number of pairs ({n})")

    base, rem = divmod(n, n_chunks)
    start = 0
    width = len(str(n_chunks))
    for chunk_idx in range(n_chunks):
        size = base + (1 if chunk_idx < rem else 0)
        chunk_lines = lines[start:start + size]
        start += size

        out_path = f"{out_prefix}_{chunk_idx + 1:0{width}d}.txt"
        with open(out_path, "w") as out:
            out.writelines(chunk_lines)
        print(f"wrote {len(chunk_lines)} pairs to {out_path}")


if __name__ == "__main__":
    main()
