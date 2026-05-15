#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/.." && pwd)

UV_BIN="${UV_BIN:-uv}"

"${UV_BIN}" run python "${SCRIPT_DIR}/01_simulate_gv_weights.py"
bash "${SCRIPT_DIR}/02_plink_gv.sh"
"${UV_BIN}" run python "${SCRIPT_DIR}/03_generate_traits.py"
"${UV_BIN}" run python "${SCRIPT_DIR}/04_make_cv_inputs_raw.py"
bash "${SCRIPT_DIR}/05_plink_gwas_raw.sh"
"${UV_BIN}" run python "${SCRIPT_DIR}/06_prepare_prscs_sst_raw.py"
bash "${SCRIPT_DIR}/07_run_prscs_raw.sh"
bash "${SCRIPT_DIR}/08_plink_score_raw.sh"
"${UV_BIN}" run python "${SCRIPT_DIR}/09_assemble_sim_df_pgs_raw.py"
"${UV_BIN}" run python "${SCRIPT_DIR}/10_adjust_phenotype_direct.py"
"${UV_BIN}" run python "${SCRIPT_DIR}/11_make_cv_inputs_adjusted.py"
bash "${SCRIPT_DIR}/12_plink_gwas_adjusted.sh"
"${UV_BIN}" run python "${SCRIPT_DIR}/13_prepare_prscs_sst.py"
bash "${SCRIPT_DIR}/14_run_prscs.sh"
bash "${SCRIPT_DIR}/15_plink_score_adjusted.sh"
"${UV_BIN}" run python "${SCRIPT_DIR}/16_assemble_sim_df_pgs_adjusted.py"
"${UV_BIN}" run python "${SCRIPT_DIR}/17_run_predinterval_comparison.py"
"${UV_BIN}" run python "${SCRIPT_DIR}/18_plot_coverage.py"

echo "Finished reproduction under ${REPO_ROOT}/outputs/run"
