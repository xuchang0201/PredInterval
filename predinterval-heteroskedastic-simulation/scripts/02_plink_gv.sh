#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)

OUTPUT_ROOT="${OUTPUT_ROOT:-${REPO_ROOT}/outputs/run}"
GENO_PREFIX="${GENO_PREFIX:-${REPO_ROOT}/data/reference/geno_82k_hm3_chr1}"
PLINK_BIN="${PLINK_BIN:-plink2}"
THREADS="${THREADS:-8}"
H2="${H2:-0.05}"
CAUSAL="${CAUSAL:-1}"
SEED_START="${SEED_START:-0}"
SEED_END="${SEED_END:-9}"

WEIGHTS_DIR="${OUTPUT_ROOT}/01-weights"
GV_DIR="${OUTPUT_ROOT}/02-gv"
mkdir -p "${GV_DIR}"

for seed in $(seq "${SEED_START}" "${SEED_END}"); do
  weight_seed=$((seed + 1))
  "${PLINK_BIN}" \
    --bfile "${GENO_PREFIX}" \
    --score "${WEIGHTS_DIR}/sim.h2${H2}.causal${CAUSAL}.seed${weight_seed}.tsv" 2 6 9 header cols=+scoresums \
    --threads "${THREADS}" \
    --out "${GV_DIR}/sim_gv.h2${H2}.causal${CAUSAL}.seed${weight_seed}"
done
