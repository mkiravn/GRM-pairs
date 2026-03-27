#!/usr/bin/env bash

set -euo pipefail

if [[ $# -lt 3 ]]; then
  cat <<'EOF'
Usage: simulate_pheno_from_bfile.sh <bfile> <grm_prefix> <out_pheno> [options]

Options:
  --h2 <x>            Trait heritability on the standardized liability scale (default: 0.4)
  --pi <x>            Fraction of causal SNPs (default: 1.0)
  --seed <x>          RNG seed (default: 1)
  --effects-out <f>   Save sampled SNP effects to this file
  --score-out <f>     Save the plink2 .sscore output prefix here
  --genetic-out <f>   Save the aligned standardized genetic component here

Environment:
  PLINK2              Path to plink2 binary (default: plink2)
EOF
  exit 1
fi

BFILE="$1"
GRM_PREFIX="$2"
OUT_PHENO="$3"
shift 3

H2="0.4"
PI="1.0"
SEED="1"
EFFECTS_OUT=""
SCORE_OUT=""
GENETIC_OUT=""
PLINK2="${PLINK2:-plink2}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --h2)
      H2="$2"
      shift 2
      ;;
    --pi)
      PI="$2"
      shift 2
      ;;
    --seed)
      SEED="$2"
      shift 2
      ;;
    --effects-out)
      EFFECTS_OUT="$2"
      shift 2
      ;;
    --score-out)
      SCORE_OUT="$2"
      shift 2
      ;;
    --genetic-out)
      GENETIC_OUT="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

effects_file="${EFFECTS_OUT:-$tmpdir/effects.tsv}"
score_prefix="${SCORE_OUT:-$tmpdir/score}"

Rscript "$SCRIPT_DIR/sample_effect_sizes_from_bim.R" \
  "${BFILE}.bim" \
  "$effects_file" \
  --h2 "$H2" \
  --pi "$PI" \
  --seed "$SEED"

"$PLINK2" \
  --bfile "$BFILE" \
  --score "$effects_file" 1 2 3 header-read variance-standardize cols=maybefid,maybesid,scoresums \
  --out "$score_prefix"

score_file="${score_prefix}.sscore"

score_args=()
if [[ -n "$GENETIC_OUT" ]]; then
  score_args+=(--genetic-out "$GENETIC_OUT")
fi

Rscript "$SCRIPT_DIR/simulate_pheno_from_score.R" \
  "$score_file" \
  "$GRM_PREFIX" \
  "$OUT_PHENO" \
  --h2 "$H2" \
  --seed "$SEED" \
  "${score_args[@]}"

echo "Wrote score-based phenotype to $OUT_PHENO"
