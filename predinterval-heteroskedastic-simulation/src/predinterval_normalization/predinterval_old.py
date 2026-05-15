from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from tqdm import tqdm

from predinterval_normalization.quant_groups import assign_extreme_quant_groups, enumerate_gamma_values, gamma_tag

WO_H2S = [0.2, 0.5, 0.8]
WO_CAUSAL_PCTS = [0.001, 0.01, 0.1, 1.0]
N_FOLDS = 5
FOLD_SIZE = 10_000
TRAIN_SIZE = N_FOLDS * FOLD_SIZE
TEST_SIZE = 10_000


def _coverage(y: pd.Series, lower: pd.Series, upper: pd.Series) -> float:
    return float(y.between(lower, upper).mean())


def _r2(x: pd.Series, y: pd.Series) -> float:
    return float(stats.pearsonr(x, y).statistic**2)


def _resolve_wocontext_path(input_dir: Path, h2: float, causal_percent: float, seed: int) -> Path:
    candidates = [
        input_dir / f"sim_df_PGS.h2{gamma_tag(h2)}.causal{gamma_tag(causal_percent)}.seed{seed}.tsv",
        input_dir / f"sim_df_PGS.train5000.h2{gamma_tag(h2)}.causal{gamma_tag(causal_percent)}.seed{seed}.tsv",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(
        f"Could not find released wocontext TSV for h2={h2}, causal={causal_percent}, seed={seed} under {input_dir}"
    )


def _build_predinterval_predictions(
    df_train: pd.DataFrame,
    df_test: pd.DataFrame,
    alphas: tuple[float, ...],
    row_chunk_size: int = 200,
) -> dict[float, pd.DataFrame]:
    ricv: list[np.ndarray] = []
    for fold in range(N_FOLDS):
        current = df_train.iloc[fold * FOLD_SIZE : (fold + 1) * FOLD_SIZE]
        ricv.append(np.abs(current["y"] - current[f"PGS_{fold}"]).to_numpy())

    pgs = df_test[[f"PGS_{fold}" for fold in range(N_FOLDS)]].to_numpy()
    n_test = len(df_test)
    upper_qs = {alpha: np.empty(n_test, dtype=float) for alpha in alphas}
    lower_qs = {alpha: np.empty(n_test, dtype=float) for alpha in alphas}
    lower_alphas = tuple(1 - alpha for alpha in alphas)
    for start in range(0, n_test, row_chunk_size):
        end = min(start + row_chunk_size, n_test)
        chunk = pgs[start:end]

        upper_matrix = np.hstack([chunk[:, [fold]] + ricv[fold][None, :] for fold in range(N_FOLDS)])
        upper_values = np.quantile(upper_matrix, alphas, axis=1)
        for idx, alpha in enumerate(alphas):
            upper_qs[alpha][start:end] = upper_values[idx]

        lower_matrix = np.hstack([chunk[:, [fold]] - ricv[fold][None, :] for fold in range(N_FOLDS)])
        lower_values = np.quantile(lower_matrix, lower_alphas, axis=1)
        for idx, alpha in enumerate(alphas):
            lower_qs[alpha][start:end] = lower_values[idx]

    predictions = {}
    mean = pgs.mean(axis=1)
    for alpha in alphas:
        z_score = stats.norm.ppf(0.5 + alpha / 2)
        std = (upper_qs[alpha] - lower_qs[alpha]) / z_score / 2
        pred = df_test.copy()
        pred["mean"] = mean
        pred["std"] = std
        pred["lower"] = lower_qs[alpha]
        pred["upper"] = upper_qs[alpha]
        predictions[alpha] = pred
    return predictions


def _run_one_wocontext(input_dir: str, h2: float, causal_percent: float, seed: int, alpha: float) -> dict[str, float | int | str]:
    path = _resolve_wocontext_path(Path(input_dir), h2, causal_percent, seed)
    df = pd.read_csv(path, sep="\t")
    df_train = df.iloc[:TRAIN_SIZE].copy()
    df_test = df.iloc[TRAIN_SIZE : TRAIN_SIZE + TEST_SIZE].copy()
    pred = _build_predinterval_predictions(df_train=df_train, df_test=df_test, alphas=(alpha,))[alpha]
    return {
        "adjust": "PredInterval",
        "n": causal_percent,
        "seed": seed,
        "col": "marginal",
        "coverage": _coverage(pred["y"], pred["lower"], pred["upper"]),
        "r2": _r2(pred["mean"], pred["y"]),
        "length": float(pred["std"].mean()),
    }


def run_wocontext(input_dir: Path, output_dir: Path, n_sim: int, alpha: float, jobs: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for h2 in WO_H2S:
        rows: list[dict[str, float | int | str]] = []
        tasks = [(str(input_dir), h2, causal_percent, seed, alpha) for causal_percent in WO_CAUSAL_PCTS for seed in range(n_sim)]
        desc = f"PredInterval wocontext h2={h2}"
        if jobs == 1:
            for task in tqdm(tasks, desc=desc):
                rows.append(_run_one_wocontext(*task))
        else:
            max_workers = min(jobs, len(tasks), os.cpu_count() or jobs)
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [executor.submit(_run_one_wocontext, *task) for task in tasks]
                for future in tqdm(as_completed(futures), total=len(futures), desc=desc):
                    rows.append(future.result())
        rows.sort(key=lambda row: (row["n"], row["seed"]))
        pd.DataFrame(rows).to_csv(output_dir / f"f1_predinterval.{h2}.tsv", sep="\t", index=False)


def _delta_r2_rows(df: pd.DataFrame) -> list[float]:
    df = df.copy()
    df["quant_q"] = assign_extreme_quant_groups(df["quant"])
    df = df.iloc[-TEST_SIZE:]
    top = df[df["quant_q"] == 9]
    bottom = df[df["quant_q"] == 0]
    values = []
    for fold in range(N_FOLDS):
        overall = stats.spearmanr(df[f"PGS_{fold}"], df["y"]).statistic**2
        top_r2 = stats.spearmanr(top[f"PGS_{fold}"], top["y"]).statistic**2
        bottom_r2 = stats.spearmanr(bottom[f"PGS_{fold}"], bottom["y"]).statistic**2
        values.append((top_r2 - bottom_r2) / overall)
    return values


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
    predictions = _build_predinterval_predictions(df_train=df_train, df_test=df_test, alphas=(0.95, 0.5))
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
    desc = "PredInterval wcontext"
    if jobs == 1:
        iterator = ( _run_one_wcontext(*task) for task in tasks )
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
        description="Reproduce PredInterval Figure 1/2 simulation results from released sim_df_PGS intermediates."
    )
    parser.add_argument("--dataset", choices=["wocontext", "wcontext", "all"], default="all")
    parser.add_argument("--wocontext-input-dir", type=Path, default=Path("extracted/df_woContext"))
    parser.add_argument("--wcontext-input-dir", type=Path, default=Path("extracted/df_wContext"))
    parser.add_argument("--output-root", type=Path, default=Path("outputs/predinterval"))
    parser.add_argument("--n-sim", type=int, default=30)
    parser.add_argument("--wocontext-alpha", type=float, default=0.95)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--gamma-values", type=float, nargs="*", default=None)
    args = parser.parse_args()

    if args.dataset in {"wocontext", "all"}:
        run_wocontext(
            input_dir=args.wocontext_input_dir,
            output_dir=args.output_root / "wocontext",
            n_sim=args.n_sim,
            alpha=args.wocontext_alpha,
            jobs=args.jobs,
        )
    if args.dataset in {"wcontext", "all"}:
        run_wcontext(
            input_dir=args.wcontext_input_dir,
            output_dir=args.output_root / "wcontext",
            n_sim=args.n_sim,
            jobs=args.jobs,
            gamma_values=args.gamma_values,
            h2=0.5,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
