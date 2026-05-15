from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from predinterval_normalization.common import (
    add_common_generation_args,
    append_suffix,
    gv_prefix,
    stage_dir,
    trait_path,
)


def _standardize(values: pd.Series) -> pd.Series:
    return (values - values.mean()) / values.std(ddof=1)


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate the original continuous wcontext traits that match the released extracted setup.")
    add_common_generation_args(parser)
    args = parser.parse_args()

    trait_dir = stage_dir(args.output_root, "03-traits")
    trait_dir.mkdir(parents=True, exist_ok=True)
    qc_rows: list[dict[str, float | int]] = []

    fam = pd.read_csv(args.geno_prefix.with_suffix(".fam"), header=None, sep="\t")
    fam = fam[[0, 1]].copy()
    fam.columns = ["FID", "IID"]

    for gamma_value in tqdm(args.gamma_values, desc="gamma"):
        gamma_scalar = gamma_value / 100.0
        rng = np.random.RandomState(42)
        quant = pd.Series(rng.normal(0, 1, size=len(fam)))
        quant_q = pd.qcut(quant, q=10).cat.codes

        for seed_index in args.seed_indices:
            output_path = trait_path(args.output_root, args.h2, args.causal, gamma_value, seed_index)
            if output_path.exists() and not args.overwrite:
                continue

            gv = pd.read_csv(
                append_suffix(gv_prefix(args.output_root, args.h2, args.causal, seed_index), ".sscore"),
                sep="\t",
            )
            df = fam.copy()
            df["quant"] = quant.to_numpy()
            df["GV"] = _standardize(gv["SCORE1_SUM"])

            scale = np.sqrt(np.exp(np.log(1 / args.h2 - 1) + df["quant"] * gamma_scalar))
            y = pd.Series(rng.normal(loc=df["GV"], scale=scale))
            y = _standardize(y)

            df["y"] = y
            df = df.sample(frac=1, random_state=seed_index).iloc[:60000].copy()
            df.to_csv(output_path, sep="\t", index=False)

            qc_rows.append(
                {
                    "gamma": gamma_value,
                    "gamma_scalar": gamma_scalar,
                    "seed": seed_index,
                    "raw_y_var_all": float(y.var(ddof=1)),
                    "raw_y_var_q0": float(y.loc[quant_q == 0].var(ddof=1)),
                    "raw_y_var_q9": float(y.loc[quant_q == 9].var(ddof=1)),
                }
            )

    if qc_rows:
        pd.DataFrame(qc_rows).sort_values(["gamma", "seed"]).to_csv(trait_dir / "trait_qc.tsv", sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
