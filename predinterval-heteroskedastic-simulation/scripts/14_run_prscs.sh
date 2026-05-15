#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)

OUTPUT_ROOT="${OUTPUT_ROOT:-${REPO_ROOT}/outputs/run}"
GENO_PREFIX="${GENO_PREFIX:-${REPO_ROOT}/data/reference/geno_82k_hm3_chr1}"
PYTHON_BIN="${PYTHON_BIN:-${REPO_ROOT}/.venv/bin/python}"
PRSCS_PY="${PRSCS_PY:-${REPO_ROOT}/third_party/PRScs/PRScs.py}"
PRSCS_REF_DIR="${PRSCS_REF_DIR:-${REPO_ROOT}/data/reference/prscs_ref/ldblk_1kg_eur}"
H2="${H2:-0.05}"
CAUSAL="${CAUSAL:-1}"
SEED_START="${SEED_START:-0}"
SEED_END="${SEED_END:-9}"
GAMMAS="${GAMMAS:-100}"
N_GWAS="${N_GWAS:-40000}"
CHROM="${CHROM:-1}"
N_ITER="${N_ITER:-1000}"
N_BURNIN="${N_BURNIN:-500}"
THIN="${THIN:-5}"
PRSCS_SEED="${PRSCS_SEED:-}"

if [[ ! -f "${PRSCS_PY}" || ! -d "${PRSCS_REF_DIR}" ]]; then
  echo "Set PRSCS_PY and PRSCS_REF_DIR before running Stage 14." >&2
  exit 1
fi

SST_DIR="${OUTPUT_ROOT}/13-sst"
PRSCS_DIR="${OUTPUT_ROOT}/14-prscs"
mkdir -p "${PRSCS_DIR}"

for gamma in ${GAMMAS}; do
  for seed in $(seq "${SEED_START}" "${SEED_END}"); do
    for cv in 0 1 2 3 4; do
      cmd=(
        "${PYTHON_BIN}" "${PRSCS_PY}"
        --ref_dir="${PRSCS_REF_DIR}"
        --bim_prefix="${GENO_PREFIX}"
        --sst_file="${SST_DIR}/gwas.h2${H2}.causal${CAUSAL}.gamma${gamma}.seed${seed}.cv${cv}.sst"
        --n_gwas="${N_GWAS}"
        --n_iter="${N_ITER}"
        --n_burnin="${N_BURNIN}"
        --thin="${THIN}"
        --chrom="${CHROM}"
        --out_dir="${PRSCS_DIR}/pgs.h2${H2}.causal${CAUSAL}.gamma${gamma}.seed${seed}.cv${cv}"
      )
      if [[ -n "${PRSCS_SEED}" ]]; then
        cmd+=(--seed="${PRSCS_SEED}")
      fi
      "${cmd[@]}"
    done
  done
done
