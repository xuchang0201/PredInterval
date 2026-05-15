from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from predinterval_normalization.common import DEFAULT_GAMMA_VALUES, DEFAULT_H2, DEFAULT_OUTPUT_ROOT
from predinterval_normalization.predinterval_new import run_wcontext as run_predinterval_new
from predinterval_normalization.predinterval_old import run_wcontext as run_predinterval_old


def _read_metric_table(path: Path, metric: str, state: str, method_label: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, names=["gamma_index", "method_group", "seed", "value"])
    df["metric"] = metric
    df["state"] = state
    df["method"] = method_label
    df["group"] = df["method_group"].str.replace(r"^PredInterval_", "", regex=True)
    return df[["state", "method", "group", "gamma_index", "seed", "metric", "value"]]


def _collect_output(root: Path, state: str, method_label: str) -> pd.DataFrame:
    specs = [
        ("f2_alpha95_predinterval.tsv", "coverage95"),
        ("f2_alpha50_predinterval.tsv", "coverage50"),
        ("f2_length_predinterval.tsv", "length"),
    ]
    frames = [_read_metric_table(root / filename, metric, state, method_label) for filename, metric in specs]
    return pd.concat(frames, ignore_index=True)


def _write_summaries(output_root: Path, long_df: pd.DataFrame, gamma_values: list[float]) -> None:
    gamma_map = {index: value for index, value in enumerate(gamma_values)}
    long_df = long_df.copy()
    long_df["gamma_value"] = long_df["gamma_index"].map(gamma_map)
    long_df.to_csv(output_root / "coverage_long.tsv", sep="\t", index=False)

    summary = (
        long_df.groupby(["state", "method", "group", "metric", "gamma_index", "gamma_value"])["value"]
        .agg(["mean", "std"])
        .reset_index()
        .sort_values(["state", "metric", "method", "group"])
        .round(6)
    )
    summary.to_csv(output_root / "method_summary.tsv", sep="\t", index=False)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the original and covariate-aware PredInterval methods on raw and normalized sim_df_PGS files."
    )
    parser.add_argument("--raw-input-dir", type=Path, default=DEFAULT_OUTPUT_ROOT / "09-sim_df_pgs")
    parser.add_argument("--adjusted-input-dir", type=Path, default=DEFAULT_OUTPUT_ROOT / "16-sim_df_pgs")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT / "17-predinterval-comparison")
    parser.add_argument("--n-sim", type=int, default=10)
    parser.add_argument("--gamma-values", type=float, nargs="*", default=[float(value) for value in DEFAULT_GAMMA_VALUES])
    parser.add_argument("--h2", type=float, default=DEFAULT_H2)
    parser.add_argument("--old-jobs", type=int, default=1)
    parser.add_argument("--new-jobs", type=int, default=1)
    args = parser.parse_args()

    raw_old_root = args.output_root / "raw" / "predinterval_old"
    raw_new_root = args.output_root / "raw" / "predinterval_new"
    adjusted_old_root = args.output_root / "adjusted" / "predinterval_old"
    adjusted_new_root = args.output_root / "adjusted" / "predinterval_new"

    run_predinterval_old(
        input_dir=args.raw_input_dir,
        output_dir=raw_old_root,
        n_sim=args.n_sim,
        jobs=args.old_jobs,
        gamma_values=args.gamma_values,
        h2=args.h2,
    )
    run_predinterval_new(
        input_dir=args.raw_input_dir,
        output_dir=raw_new_root,
        n_sim=args.n_sim,
        jobs=args.new_jobs,
        gamma_values=args.gamma_values,
        h2=args.h2,
    )
    run_predinterval_old(
        input_dir=args.adjusted_input_dir,
        output_dir=adjusted_old_root,
        n_sim=args.n_sim,
        jobs=args.old_jobs,
        gamma_values=args.gamma_values,
        h2=args.h2,
    )
    run_predinterval_new(
        input_dir=args.adjusted_input_dir,
        output_dir=adjusted_new_root,
        n_sim=args.n_sim,
        jobs=args.new_jobs,
        gamma_values=args.gamma_values,
        h2=args.h2,
    )

    long_df = pd.concat(
        [
            _collect_output(raw_old_root, "Raw phenotype", "PGS Only"),
            _collect_output(raw_new_root, "Raw phenotype", "PGS + Covariates"),
            _collect_output(adjusted_old_root, "Normalized phenotype", "PGS Only"),
            _collect_output(adjusted_new_root, "Normalized phenotype", "PGS + Covariates"),
        ],
        ignore_index=True,
    )
    args.output_root.mkdir(parents=True, exist_ok=True)
    _write_summaries(args.output_root, long_df, args.gamma_values)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
