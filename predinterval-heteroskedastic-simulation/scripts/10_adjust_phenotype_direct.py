from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
import fcntl
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from predinterval_normalization.common import (
    TRAIN_SIZE,
    add_common_generation_args,
    adjusted_trait_path,
    format_tag,
    sim_df_pgs_path,
    stage_dir,
)

EPS = 1e-8


def _fit_linear_model(design: np.ndarray, response: np.ndarray) -> np.ndarray:
    beta_hat, *_ = np.linalg.lstsq(design, response, rcond=None)
    return beta_hat


def _quant_codes(values: pd.Series) -> pd.Series:
    return pd.qcut(values, q=10, labels=False, duplicates="drop")


@contextmanager
def _locked_file(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _qc_shard_path(output_dir: Path, h2: float, causal: float, gamma_value: int, seed_index: int) -> Path:
    return output_dir / (
        f"adjustment_qc.h2{format_tag(h2)}.causal{format_tag(causal)}.gamma{gamma_value}.seed{seed_index}.tsv"
    )


def _persist_qc_rows(output_dir: Path, qc_rows: list[dict[str, float | int | str]]) -> None:
    if not qc_rows:
        return

    lock_path = output_dir / ".adjustment_qc.lock"
    with _locked_file(lock_path):
        for row in qc_rows:
            shard_path = _qc_shard_path(
                output_dir,
                float(row["h2"]),
                float(row["causal"]),
                int(row["gamma"]),
                int(row["seed"]),
            )
            pd.DataFrame([row]).to_csv(shard_path, sep="\t", index=False)

        shard_paths = sorted(output_dir.glob("adjustment_qc.h2*.causal*.gamma*.seed*.tsv"))
        merged = pd.concat((pd.read_csv(path, sep="\t") for path in shard_paths), ignore_index=True)
        merged = merged.sort_values(["h2", "causal", "gamma", "seed"]).reset_index(drop=True)
        merged.to_csv(output_dir / "adjustment_qc.tsv", sep="\t", index=False)


def _transform_one(
    output_root: str,
    h2: float,
    causal: float,
    gamma_value: int,
    seed_index: int,
    overwrite: bool,
) -> dict[str, float | int | str]:
    input_path = sim_df_pgs_path(Path(output_root), h2, causal, gamma_value, seed_index)
    output_path = adjusted_trait_path(Path(output_root), h2, causal, gamma_value, seed_index)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep="\t")
    if len(df) < TRAIN_SIZE:
        raise ValueError(f"Expected at least {TRAIN_SIZE} rows in {input_path}, found {len(df)}")

    fit_df = df.iloc[:TRAIN_SIZE].copy()
    fit_quant = fit_df["quant"].to_numpy(dtype=float)
    fit_y = fit_df["y"].to_numpy(dtype=float)
    fit_design = np.column_stack([np.ones_like(fit_quant), fit_quant, fit_quant**2])

    mu_beta = _fit_linear_model(fit_design, fit_y)
    fit_mean_hat = fit_design @ mu_beta
    resid_mu = fit_y - fit_mean_hat

    log_resid_sq = np.log(np.maximum(np.square(resid_mu), EPS))
    var_beta = _fit_linear_model(fit_design, log_resid_sq)

    quant = df["quant"].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=float)
    design = np.column_stack([np.ones_like(quant), quant, quant**2])

    mean_hat = design @ mu_beta
    log_var_hat = design @ var_beta
    sd_hat = np.maximum(np.exp(0.5 * log_var_hat), np.sqrt(EPS))
    y_new = (y - mean_hat) / sd_hat

    transformed = df[["FID", "IID", "quant", "GV"]].copy()
    transformed["y"] = y_new
    if overwrite or not output_path.exists():
        transformed.to_csv(output_path, sep="\t", index=False)

    raw_q = _quant_codes(df["quant"])
    adj_q = _quant_codes(transformed["quant"])
    return {
        "h2": h2,
        "causal": causal,
        "gamma": gamma_value,
        "seed": seed_index,
        "input_file": input_path.name,
        "output_file": output_path.name,
        "transform": "fit_mu: lm(y~quant+quant^2); fit_var: lm(log(resid_mu^2)~quant+quant^2); y_new=(y-mu_hat)/exp(0.5*logvar_hat)",
        "fit_rows": TRAIN_SIZE,
        "apply_rows": int(len(df)),
        "raw_y_var_all": float(df["y"].var(ddof=1)),
        "raw_y_var_q0": float(df.loc[raw_q == 0, "y"].var(ddof=1)),
        "raw_y_var_q9": float(df.loc[raw_q == 9, "y"].var(ddof=1)),
        "adj_y_var_all": float(transformed["y"].var(ddof=1)),
        "adj_y_var_q0": float(transformed.loc[adj_q == 0, "y"].var(ddof=1)),
        "adj_y_var_q9": float(transformed.loc[adj_q == 9, "y"].var(ddof=1)),
        "adj_y_mean_all": float(transformed["y"].mean()),
        "adj_y_sd_all": float(transformed["y"].std(ddof=1)),
        "train_adj_y_mean": float(pd.Series(y_new[:TRAIN_SIZE]).mean()),
        "train_adj_y_sd": float(pd.Series(y_new[:TRAIN_SIZE]).std(ddof=1)),
        "test_adj_y_mean": float(pd.Series(y_new[TRAIN_SIZE:]).mean()),
        "test_adj_y_sd": float(pd.Series(y_new[TRAIN_SIZE:]).std(ddof=1)),
        "mu_beta0": float(mu_beta[0]),
        "mu_beta1_quant": float(mu_beta[1]),
        "mu_beta2_quant2": float(mu_beta[2]),
        "var_beta0": float(var_beta[0]),
        "var_beta1_quant": float(var_beta[1]),
        "var_beta2_quant2": float(var_beta[2]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Adjust stage-09 phenotypes by fitting direct-y mean/log-variance models with quant and quant^2 on the first 50000 rows only, then apply the learned transform to all 60000."
    )
    parser.add_argument("--jobs", type=int, default=1)
    add_common_generation_args(parser)
    args = parser.parse_args()

    output_dir = stage_dir(args.output_root, "10-adjusted-traits")
    output_dir.mkdir(parents=True, exist_ok=True)
    tasks = [
        (str(args.output_root), args.h2, args.causal, gamma_value, seed_index, args.overwrite)
        for gamma_value in args.gamma_values
        for seed_index in args.seed_indices
    ]
    qc_rows: list[dict[str, float | int | str]] = []

    if args.jobs == 1:
        for task in tqdm(tasks, desc="adjust phenotypes"):
            qc_rows.append(_transform_one(*task))
    else:
        max_workers = min(args.jobs, len(tasks), os.cpu_count() or args.jobs)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_transform_one, *task) for task in tasks]
            for future in tqdm(as_completed(futures), total=len(futures), desc="adjust phenotypes"):
                qc_rows.append(future.result())

    _persist_qc_rows(output_dir, qc_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
