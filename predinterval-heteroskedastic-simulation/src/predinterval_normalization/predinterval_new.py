from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import GammaRegressor, LinearRegression
from tqdm import tqdm

from predinterval_normalization.quant_groups import (
    assign_extreme_quant_groups,
    enumerate_gamma_values,
    gamma_tag,
    select_nonredundant_columns,
)
from predinterval_normalization.predinterval_old import (
    FOLD_SIZE,
    N_FOLDS,
    TEST_SIZE,
    TRAIN_SIZE,
    _coverage,
    _delta_r2_rows,
)

EPS = 1e-8


def _add_covariates(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["quant_sq"] = out["quant"] ** 2
    return out


def _assign_honest_training_scores(df_train: pd.DataFrame) -> pd.DataFrame:
    out = df_train.copy()
    fold_ids = np.repeat(np.arange(N_FOLDS), FOLD_SIZE)
    pgs_matrix = out[[f"PGS_{fold}" for fold in range(N_FOLDS)]].to_numpy()
    out["fold"] = fold_ids
    out["SCORESUM"] = pgs_matrix[np.arange(len(out)), fold_ids]
    return out


def _feature_frame(df: pd.DataFrame) -> pd.DataFrame:
    columns = select_nonredundant_columns(df, ["SCORESUM", "quant", "quant_sq"])
    return df[columns]


def _conformal_scale(nonconformity: np.ndarray, conf_level: float) -> float:
    n = len(nonconformity)
    rank = int(np.ceil(conf_level * (n + 1)))
    rank = min(max(rank, 1), n)
    return float(np.sort(nonconformity)[rank - 1])


def _build_covariate_predictions(
    df_train: pd.DataFrame,
    df_test: pd.DataFrame,
    alphas: tuple[float, ...],
) -> dict[float, pd.DataFrame]:
    train = _assign_honest_training_scores(_add_covariates(df_train))
    test = _add_covariates(df_test)

    nonconformity = np.empty(len(train), dtype=float)
    mu_hat_test = np.empty((len(test), N_FOLDS), dtype=float)
    sigma_hat_test = np.empty((len(test), N_FOLDS), dtype=float)

    for fold in range(N_FOLDS):
        train_mask = train["fold"] != fold
        held_mask = ~train_mask
        train_fold = train.loc[train_mask].copy()
        held_fold = train.loc[held_mask].copy()

        x_train = _feature_frame(train_fold)
        mean_model = LinearRegression().fit(x_train, train_fold["y"])
        train_abs_resid = np.maximum(np.abs(train_fold["y"] - mean_model.predict(x_train)), EPS)
        var_model = GammaRegressor(alpha=1e-6, max_iter=1000).fit(x_train, train_abs_resid)

        x_held = _feature_frame(held_fold)
        mu_held = mean_model.predict(x_held)
        sigma_held = np.maximum(var_model.predict(x_held), EPS)
        nonconformity[held_mask.to_numpy()] = np.abs(held_fold["y"] - mu_held) / sigma_held

        test_fold = test.copy()
        test_fold["SCORESUM"] = test_fold[f"PGS_{fold}"]
        x_test = _feature_frame(test_fold)
        mu_hat_test[:, fold] = mean_model.predict(x_test)
        sigma_hat_test[:, fold] = np.maximum(var_model.predict(x_test), EPS)

    mu_hat = mu_hat_test.mean(axis=1)
    sigma_hat = sigma_hat_test.mean(axis=1)

    predictions: dict[float, pd.DataFrame] = {}
    for alpha in alphas:
        conformal_scale = _conformal_scale(nonconformity, conf_level=alpha)
        pred = test.copy()
        pred["mean"] = mu_hat
        pred["lower"] = mu_hat - sigma_hat * conformal_scale
        pred["upper"] = mu_hat + sigma_hat * conformal_scale
        z_score = stats.norm.ppf(0.5 + alpha / 2)
        pred["std"] = (pred["upper"] - pred["lower"]) / (2 * z_score)
        predictions[alpha] = pred
    return predictions


def _run_one_wcontext(
    input_dir: str,
    h2: float,
    gamma_index: int,
    gamma_value: float,
    seed: int,
) -> tuple[list[list[float | int | str]], list[list[float | int | str]], list[list[float | int | str]], list[float]]:
    df = pd.read_csv(
        Path(input_dir) / f"sim_df_PGS.h2{gamma_tag(h2)}.causal1.gamma{gamma_tag(gamma_value)}.seed{seed}.tsv",
        sep="\t",
    )
    df["quant_q"] = assign_extreme_quant_groups(df["quant"])
    delta_r2_values = _delta_r2_rows(df)

    df_train = df.iloc[:TRAIN_SIZE].copy()
    df_test = df.iloc[TRAIN_SIZE : TRAIN_SIZE + TEST_SIZE].copy()
    groups = {
        "PredInterval_marginal": df_test,
        "PredInterval_quant_q_0": df_test[df_test["quant_q"] == 0],
        "PredInterval_quant_q_9": df_test[df_test["quant_q"] == 9],
    }

    coverage_95_rows: list[list[float | int | str]] = []
    coverage_50_rows: list[list[float | int | str]] = []
    length_rows: list[list[float | int | str]] = []
    predictions = _build_covariate_predictions(df_train=df_train, df_test=df_test, alphas=(0.95, 0.5))
    for alpha in (0.95, 0.5):
        pred = predictions[alpha]
        pred_groups = {name: pred.loc[current.index] for name, current in groups.items()}
        current_rows = coverage_95_rows if alpha == 0.95 else coverage_50_rows
        for method, current in pred_groups.items():
            current_rows.append(
                [
                    gamma_index,
                    method,
                    seed,
                    _coverage(current["y"], current["lower"], current["upper"]),
                ]
            )
            if alpha == 0.95:
                length_rows.append([gamma_index, method, seed, float(current["std"].mean())])

    return coverage_95_rows, coverage_50_rows, length_rows, delta_r2_values


def run_wcontext(
    input_dir: Path,
    output_dir: Path,
    n_sim: int,
    jobs: int,
    gamma_values: list[float] | None = None,
    h2: float = 0.5,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    coverage_95_rows: list[list[float | int | str]] = []
    coverage_50_rows: list[list[float | int | str]] = []
    length_rows: list[list[float | int | str]] = []
    gamma_pairs = enumerate_gamma_values(gamma_values)
    delta_r2_map: dict[int, list[float]] = {gamma_index: [] for gamma_index, _ in gamma_pairs}

    tasks = [(str(input_dir), h2, gamma_index, gamma_value, seed) for gamma_index, gamma_value in gamma_pairs for seed in range(n_sim)]
    desc = "PredInterval wcontext covariates"
    if jobs == 1:
        iterator = (_run_one_wcontext(*task) for task in tasks)
        progress = tqdm(iterator, total=len(tasks), desc=desc)
        for task, result in zip(tasks, progress, strict=True):
            rows95, rows50, current_lengths, delta_values = result
            coverage_95_rows.extend(rows95)
            coverage_50_rows.extend(rows50)
            length_rows.extend(current_lengths)
            delta_r2_map[task[2]].extend(delta_values)
    else:
        max_workers = min(jobs, len(tasks), os.cpu_count() or jobs)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {executor.submit(_run_one_wcontext, *task): task for task in tasks}
            for future in tqdm(as_completed(future_to_task), total=len(future_to_task), desc=desc):
                gamma_index = future_to_task[future][2]
                rows95, rows50, current_lengths, delta_values = future.result()
                coverage_95_rows.extend(rows95)
                coverage_50_rows.extend(rows50)
                length_rows.extend(current_lengths)
                delta_r2_map[gamma_index].extend(delta_values)

    coverage_95_rows.sort(key=lambda row: (row[0], row[1], row[2]))
    coverage_50_rows.sort(key=lambda row: (row[0], row[1], row[2]))
    length_rows.sort(key=lambda row: (row[0], row[1], row[2]))
    delta_r2_rows = [[gamma_index, abs(float(np.mean(values))) * 100] for gamma_index, values in sorted(delta_r2_map.items())]

    pd.DataFrame(coverage_95_rows).to_csv(
        output_dir / "f2_alpha95_predinterval.tsv", sep="\t", index=False, header=False
    )
    pd.DataFrame(coverage_50_rows).to_csv(
        output_dir / "f2_alpha50_predinterval.tsv", sep="\t", index=False, header=False
    )
    pd.DataFrame(length_rows).to_csv(output_dir / "f2_length_predinterval.tsv", sep="\t", index=False, header=False)
    pd.DataFrame(delta_r2_rows).to_csv(output_dir / "f2_deltaR2.tsv", sep="\t", index=False, header=False)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the covariate-aware PredInterval variant on wcontext sim_df_PGS inputs."
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-sim", type=int, default=30)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--gamma-values", type=float, nargs="*", default=None)
    parser.add_argument("--h2", type=float, default=0.5)
    args = parser.parse_args()

    run_wcontext(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        n_sim=args.n_sim,
        jobs=args.jobs,
        gamma_values=args.gamma_values,
        h2=args.h2,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
