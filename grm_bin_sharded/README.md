# grm_bin_sharded

Bins a GRM by relatedness and computes a delete-block jackknife SE (blocks
of **individuals**, never pairs) on the binned phenotype cross-product —
same statistics as `../full_grm_bin/blocked_grm_bin_means.cpp` — but split
across however many row-chunked shards the GRM was computed in (e.g.
`plink --make-grm-bin --parallel k n`), without ever assembling the full
dense `.grm.bin`.

## Why

`full_grm_bin`/`grm_bin` assume the entire `.grm.bin` fits in one file you
stream once. At scale that full N×N matrix may be too large to want to
build or move around at all. But the accumulation itself — a sum, a sum of
squares, and a count per (relatedness bin, jackknife block) — is additive:
you can accumulate those sums separately per shard and just add the small
per-shard accumulator matrices together afterward. This tool does exactly
that: `accumulate` runs on one shard and produces a tiny sums file;
`merge` adds several of those together and produces the same
`.full.tsv`/`.jk.tsv` output `blocked_grm_bin_means` would have produced
from the whole matrix in one pass.

## Building

```bash
make
```

## Usage

### 1. `accumulate` — one job per shard

```bash
./grm_shard_tool accumulate \
    --grm-id mydata.grm.id \
    --shard  mydata.grm.bin.3 \
    --parallel 3 20 \
    --pheno  aligned_pheno.txt \
    --bins   bins.txt \
    --nblocks 50 --seed 1 \
    --out    shard_3.acc.tsv
```

- `--grm-id`: the one shared `.grm.id` for the whole matrix (all N
  individuals, in plink's row/col order).
- `--shard`: this shard's raw `.grm.bin.k` file (just the floats for its
  row range — nothing else in it identifies which rows they are).
- `--parallel K N`: same meaning as plink's own flag — this is shard `K`
  of `N` total. Used to *estimate* which rows this shard covers; the
  estimate is then checked against the shard file's actual size and
  adjusted within a small window if needed (see "Row-range recovery"
  below) — it errors out rather than silently misaligning if it can't
  find a match.
- `--pheno`: phenotype aligned to `.grm.id` order, `FID IID value` per
  line (same convention as `prepare_pheno.R`/`prep_pheno.R`), missing as
  `NA`/`NaN`/`.`.
- `--bins`: bin edges, one `left right` pair per line (same format as
  `full_grm_bin/bins.txt`).
- `--nblocks`/`--seed`: jackknife block count and RNG seed. **Must be
  identical across every shard's `accumulate` call for a given run** —
  the block assignment is a pure function of (N, nblocks, seed), so every
  shard computes the same assignment independently with no coordination.
- `--out`: this shard's accumulator file — a small TSV of raw sums, sums
  of squares, and counts per bin (`scope=full`) and per (block, bin)
  (`scope=drop`). Not a finished answer — just cached partial totals.

### 2. `merge` — once, after all shards are done

```bash
ls shard_*.acc.tsv > acc_list.txt
./grm_shard_tool merge \
    --acc-list acc_list.txt \
    --bins bins.txt \
    --nblocks 50 \
    --out-prefix merged
```

Sums the `full`/`drop` accumulators across every shard's file and writes
`merged.full.tsv` and `merged.jk.tsv` — identical schema to
`full_grm_bin/blocked_grm_bin_means`'s own output.

### 3. `ranges` — human sanity check, no shard file needed

```bash
./grm_shard_tool ranges --n-ids 2504 --parallel 3 20
```

Prints the row range shard 3 of 20 *should* cover, by the same formula
`accumulate` uses internally, for eyeballing against plink's own log
output.

## Row-range recovery

A GRM shard file (`.grm.bin.k`) is just a flat array of floats — nothing in
it says which rows of the full matrix it covers. `accumulate` recovers the
row range from `--parallel k n` using **plink 1.9's own split formula**,
read out of its actual source (`plink_common.c: parallel_bounds()` /
`triangle_divide()`, `plink_calc.c: calc_rel()`), not a guess:

- The split boundary is chosen by balancing **off-diagonal pair counts**
  only (`vv·(vv−1) ≥ target`, i.e. `N·(N−1)` scale) — the diagonal
  (self-relatedness/IBC) is cheap to compute and just gets folded into
  whichever row it lands in when writing.
- Each row in the file still holds `row_index + 1` floats (off-diagonal
  entries, then the diagonal entry) — confirmed from `calc_rel_grm_emitn`
  / the `.grm.bin` write loop in `calc_rel()`.
- plink defines `triangle_divide(0, -1) = 1`, then snaps that back to `0`
  for the first shard's start row (so row 0, which has no off-diagonal
  entries of its own, isn't dropped). `plink_parallel_row_start()`
  reproduces this exactly, including the snap.

`accumulate` still verifies the computed boundary against the shard
file's actual size (rows must divide it evenly, no partial row at the
end) and searches a small window around it if the exact formula doesn't
match — a defensive fallback for a plink version difference or a
different tool (GCTA, plink2) producing the shards, not the primary
mechanism anymore. It errors loudly rather than silently misaligning if
nothing nearby fits.

This has been checked against plink 1.9's real algorithm (source above),
and against locally-cut shards using that exact formula (`test/`,
`split_grm_shards.py`) — not yet against byte-for-byte real plink
`--parallel` output, since no working `plink`/`plink2` build is available
on this machine (the local ones are the wrong architecture / lack AVX2).
Worth a final check against a real plink shard the first time one exists.

## Non-goals (for now)

- No `--exclude-pairs` support (present in `full_grm_bin`, not yet ported
  here).
- Plink-shard path only — no Hail equivalent yet.
