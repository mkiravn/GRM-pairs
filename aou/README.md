# AoU (All of Us) chunked GRM computation

Preliminary work for running the GRM-pairs pipeline against All of Us data,
executed as batch jobs from a laptop against Google Cloud (`gcloud batch` /
`dsub`), instead of the DNAnexus/UKB setup in `../dnanexus_examples/`.

## Why this exists

The existing tools assume one of two things:

- `grm_pairs` (plink path): a single pass over one whole `.bed/.bim/.fam`
  fileset for one pair list, in one job (see `../grm_pairs/calc_rel_pairs.c`).
- `grm_bin` / `full_grm_bin`: a full precomputed `.grm.bin`/`.grm.id` (e.g.
  from GCTA) streamed once.

Neither splits the work across multiple batch jobs. For AoU-scale data we
want to take a master ID/pair list, split it into chunks, and dispatch one
job per chunk (plink genotypes for now; Hail MatrixTable support later),
then merge the per-chunk outputs.

## Status

Plink-path chunking is implemented and tested locally against the 1000
Genomes chr1 data already in `../1kg/`. SNP-chunking and the Hail path are
stubs — decide whether either is needed once we see real AoU job timings.

## Layout

- `chunking/split_pairs.py` — splits a master pair list into N contiguous
  chunk files.
- `plink/run_chunk.sh` — runs the existing `grm_pairs` binary on one pair
  chunk against a plink bfile. This is the whole "chunk" unit of work.
- `plink/merge_chunks.sh` — concatenates per-chunk `grm_pairs` outputs back
  into one file (single header, rows in chunk order).
- `hail/grm_pairs_hail.py` — rough, unverified sketch of the same estimator
  against a Hail MatrixTable, for one ID chunk. Not run against real AoU data
  yet.
- `batch/dsub_submit_example.sh` — example `dsub` invocation for submitting
  one chunk as a Google Cloud Batch job from a laptop, analogous to
  `../dnanexus_examples/run_pipeline.sh`.
- `test/test_chunking_consistency.sh` — local preliminary test: builds a pair
  list from `../1kg/1000g_chr1`, runs `grm_pairs` once on the whole list vs.
  split+chunked+merged, and diffs the two outputs.

## Chunking axis

Chunking by ID/pair sublist (not SNP blocks): each job does one full pass
over the genotype file for its slice of pairs, and outputs need no further
combination beyond concatenation — no partial-sum merge step, unlike
SNP-chunking which would require summing numerators/n_valid across blocks
before dividing. Revisit SNP-chunking if a real AoU pass-per-chunk over the
full genotype file turns out to be the bottleneck (e.g. srWGS-scale variant
counts), since that's a case where re-reading the whole genotype file per ID
chunk becomes wasteful and splitting on SNPs instead (single genotype pass,
partial sums per chunk, merge before dividing) would pay off.

## Genotype format assumptions for AoU

- Plink path: AoU exports plink-converted filesets for array and short-read
  common-variant data; same `.bed/.bim/.fam` + `.afreq`/`.frq` assumptions as
  `../grm_pairs/README.md` should hold.
- Hail path: AoU's srWGS callset ships as Hail MatrixTables/VariantDatasets.
  `hail/grm_pairs_hail.py` is a placeholder for reading genotype dosages
  directly off a MatrixTable subset by sample ID, without an intermediate
  plink export.
