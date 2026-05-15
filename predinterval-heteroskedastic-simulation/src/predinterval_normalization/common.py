from __future__ import annotations

import argparse
import os
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_ROOT = Path(os.environ.get("OUTPUT_ROOT", PROJECT_ROOT / "outputs" / "run"))
DEFAULT_GENO_PREFIX = Path(os.environ.get("GENO_PREFIX", PROJECT_ROOT / "data" / "reference" / "geno_82k_hm3_chr1"))
DEFAULT_PRSCS_PY = Path(os.environ.get("PRSCS_PY", PROJECT_ROOT / "third_party" / "PRScs" / "PRScs.py"))
DEFAULT_PRSCS_REF_DIR = Path(
    os.environ.get("PRSCS_REF_DIR", PROJECT_ROOT / "data" / "reference" / "prscs_ref" / "ldblk_1kg_eur")
)
DEFAULT_THREADS = int(os.environ.get("THREADS", "8"))

DEFAULT_H2 = float(os.environ.get("H2", "0.05"))
DEFAULT_CAUSAL = float(os.environ.get("CAUSAL", "1.0"))
DEFAULT_GAMMA_VALUES = [int(value) for value in os.environ.get("GAMMAS", "100").split()]
DEFAULT_SEED_START = int(os.environ.get("SEED_START", "0"))
DEFAULT_SEED_END = int(os.environ.get("SEED_END", "9"))
DEFAULT_SEED_INDICES = list(range(DEFAULT_SEED_START, DEFAULT_SEED_END + 1))
DEFAULT_N_SIM = len(DEFAULT_SEED_INDICES)

N_FOLDS = 5
FOLD_SIZE = 10_000
TRAIN_SIZE = N_FOLDS * FOLD_SIZE
TEST_SIZE = 10_000
N_OUTPUT_ROWS = TRAIN_SIZE + TEST_SIZE


def format_tag(value: float) -> str:
    return f"{value:g}"


def stage_dir(output_root: Path, stage_name: str) -> Path:
    return output_root / stage_name


def append_suffix(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}{suffix}")


def weight_path(output_root: Path, h2: float, causal: float, seed_index: int) -> Path:
    return stage_dir(output_root, "01-weights") / (
        f"sim.h2{format_tag(h2)}.causal{format_tag(causal)}.seed{seed_index + 1}.tsv"
    )


def gv_prefix(output_root: Path, h2: float, causal: float, seed_index: int) -> Path:
    return stage_dir(output_root, "02-gv") / (
        f"sim_gv.h2{format_tag(h2)}.causal{format_tag(causal)}.seed{seed_index + 1}"
    )


def trait_path(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int) -> Path:
    return stage_dir(output_root, "03-traits") / (
        f"sim_df.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.tsv"
    )


def cv_pheno_path(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int, fold_index: int) -> Path:
    return stage_dir(output_root, "04-cv") / "pheno" / (
        f"pheno.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.cv{fold_index}.tsv"
    )


def cv_covar_path(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int, fold_index: int) -> Path:
    return stage_dir(output_root, "04-cv") / "covar" / (
        f"covar.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.cv{fold_index}.tsv"
    )


def gwas_prefix(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int, fold_index: int) -> Path:
    return stage_dir(output_root, "05-gwas") / (
        f"gwas.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.cv{fold_index}"
    )


def raw_sst_path(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int, fold_index: int) -> Path:
    return stage_dir(output_root, "06-sst") / (
        f"gwas.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.cv{fold_index}.sst"
    )


def raw_prscs_prefix(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int, fold_index: int) -> Path:
    return stage_dir(output_root, "07-prscs") / (
        f"pgs.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.cv{fold_index}"
    )


def raw_score_prefix(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int, fold_index: int) -> Path:
    return stage_dir(output_root, "08-scores") / (
        f"score.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.cv{fold_index}"
    )


def sim_df_pgs_path(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int) -> Path:
    return stage_dir(output_root, "09-sim_df_pgs") / (
        f"sim_df_PGS.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.tsv"
    )


def adjusted_trait_path(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int) -> Path:
    return stage_dir(output_root, "10-adjusted-traits") / (
        f"sim_df.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.tsv"
    )


def rebuild_cv_pheno_path(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int, fold_index: int) -> Path:
    return stage_dir(output_root, "11-cv") / "pheno" / (
        f"pheno.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.cv{fold_index}.tsv"
    )


def rebuild_cv_covar_path(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int, fold_index: int) -> Path:
    return stage_dir(output_root, "11-cv") / "covar" / (
        f"covar.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.cv{fold_index}.tsv"
    )


def rebuild_gwas_prefix(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int, fold_index: int) -> Path:
    return stage_dir(output_root, "12-gwas") / (
        f"gwas.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.cv{fold_index}"
    )


def rebuild_sst_path(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int, fold_index: int) -> Path:
    return stage_dir(output_root, "13-sst") / (
        f"gwas.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.cv{fold_index}.sst"
    )


def rebuild_prscs_prefix(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int, fold_index: int) -> Path:
    return stage_dir(output_root, "14-prscs") / (
        f"pgs.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.cv{fold_index}"
    )


def rebuild_score_prefix(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int, fold_index: int) -> Path:
    return stage_dir(output_root, "15-scores") / (
        f"score.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.cv{fold_index}"
    )


def rebuilt_sim_df_pgs_path(output_root: Path, h2: float, causal: float, gamma_value: int, seed_index: int) -> Path:
    return stage_dir(output_root, "16-sim_df_pgs") / (
        f"sim_df_PGS.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.tsv"
    )


def add_common_generation_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--geno-prefix", type=Path, default=DEFAULT_GENO_PREFIX)
    parser.add_argument("--h2", type=float, default=DEFAULT_H2)
    parser.add_argument("--causal", type=float, default=DEFAULT_CAUSAL)
    parser.add_argument("--seed-indices", type=int, nargs="*", default=DEFAULT_SEED_INDICES)
    parser.add_argument("--gamma-values", type=int, nargs="*", default=DEFAULT_GAMMA_VALUES)
    parser.add_argument("--overwrite", action="store_true")
