#!/usr/bin/env python3
"""Cut a dense .grm.bin into N shards using plink 1.9's actual --parallel
split point (verified against plink_common.c's parallel_bounds()/
triangle_divide() and plink_calc.c's calc_rel()), not an approximation --
so this reproduces real plink shard boundaries, and grm_shard_tool's row-
range recovery should match them at offset 0 with no search-window fallback
needed.

Usage: split_grm_shards.py <grm.bin> <n_ids> <n_shards> <out_prefix>
Writes <out_prefix>.1 .. <out_prefix>.<n_shards>.
"""
import math
import sys


def cumulative_entries(i: int) -> int:
    """Floats through row i (exclusive), diagonal included -- the file layout."""
    return i * (i + 1) // 2


def triangle_divide_off_diag(target: int) -> int:
    """plink's triangle_divide() with modif=-1: smallest v with v*(v-1) >= target."""
    if target == 0:
        return 1
    v = int(math.sqrt(target)) + 2
    while v > 1 and (v - 1) * (v - 2) >= target:
        v -= 1
    while v * (v - 1) < target:
        v += 1
    return v


def plink_parallel_row_start(n_ids: int, parallel_idx: int, n_shards: int) -> int:
    """plink's calc_rel() min_sample for shard `parallel_idx` (0-based) of
    `n_shards`, off-diagonal-balanced -- then the min_sample==1 -> 0 fix."""
    ct_tot = n_ids * (n_ids - 1)
    target = (ct_tot * parallel_idx) // n_shards
    v = triangle_divide_off_diag(target)
    return 0 if v == 1 else v


def main() -> None:
    if len(sys.argv) != 5:
        sys.exit(__doc__)

    grm_bin, n_ids_str, n_shards_str, out_prefix = sys.argv[1:5]
    n_ids = int(n_ids_str)
    n_shards = int(n_shards_str)

    total = cumulative_entries(n_ids)
    with open(grm_bin, "rb") as f:
        data = f.read()
    assert len(data) == total * 4, "grm.bin size does not match n_ids"

    row_bounds = [plink_parallel_row_start(n_ids, k, n_shards) for k in range(n_shards)]
    row_bounds.append(n_ids)

    for k in range(1, n_shards + 1):
        row_start, row_end = row_bounds[k - 1], row_bounds[k]
        start_entry = cumulative_entries(row_start)
        end_entry = cumulative_entries(row_end)
        out_path = f"{out_prefix}.{k}"
        with open(out_path, "wb") as f:
            f.write(data[start_entry * 4 : end_entry * 4])
        print(f"shard {k}/{n_shards}: rows [{row_start}, {row_end}) -> {out_path} "
              f"({end_entry - start_entry} floats)")


if __name__ == "__main__":
    main()
