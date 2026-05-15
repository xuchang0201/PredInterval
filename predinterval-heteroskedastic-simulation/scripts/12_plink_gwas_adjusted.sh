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
GAMMAS="${GAMMAS:-100}"

PHENO_DIR="${OUTPUT_ROOT}/11-cv/pheno"
COVAR_DIR="${OUTPUT_ROOT}/11-cv/covar"
GWAS_DIR="${OUTPUT_ROOT}/12-gwas"
mkdir -p "${GWAS_DIR}"

for gamma in ${GAMMAS}; do
  for seed in $(seq "${SEED_START}" "${SEED_END}"); do
    for cv in 0 1 2 3 4; do
      "${PLINK_BIN}" \
        --bfile "${GENO_PREFIX}" \
        --pheno "${PHENO_DIR}/pheno.h2${H2}.causal${CAUSAL}.gamma${gamma}.seed${seed}.cv${cv}.tsv" \
        --pheno-col-nums 3 \
        --covar "${COVAR_DIR}/covar.h2${H2}.causal${CAUSAL}.gamma${gamma}.seed${seed}.cv${cv}.tsv" \
        --covar-col-nums 3 \
        --glm hide-covar \
        --no-input-missing-phenotype \
        --threads "${THREADS}" \
        --out "${GWAS_DIR}/gwas.h2${H2}.causal${CAUSAL}.gamma${gamma}.seed${seed}.cv${cv}"
    done
  done
done
