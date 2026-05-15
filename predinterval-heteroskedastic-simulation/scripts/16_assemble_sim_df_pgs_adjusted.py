from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from predinterval_normalization.common import (
    N_FOLDS,
    add_common_generation_args,
    adjusted_trait_path,
    append_suffix,
    rebuild_score_prefix,
    rebuilt_sim_df_pgs_path,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Merge stage-15 PRS-CS PLINK scores back onto stage-10 adjusted phenotypes to create rebuilt sim_df_PGS files.")
    add_common_generation_args(parser)
    args = parser.parse_args()

    output_dir = args.output_root / "16-sim_df_pgs"
    output_dir.mkdir(parents=True, exist_ok=True)

    tasks = [(gamma_value, seed_index) for gamma_value in args.gamma_values for seed_index in args.seed_indices]
    for gamma_value, seed_index in tqdm(tasks, desc="assemble rebuilt sim_df_PGS"):
        output_path = rebuilt_sim_df_pgs_path(args.output_root, args.h2, args.causal, gamma_value, seed_index)
        if output_path.exists() and not args.overwrite:
            continue

        df = pd.read_csv(adjusted_trait_path(args.output_root, args.h2, args.causal, gamma_value, seed_index), sep="\t")
        for fold_index in range(N_FOLDS):
            pgs = pd.read_csv(
                append_suffix(
                    rebuild_score_prefix(args.output_root, args.h2, args.causal, gamma_value, seed_index, fold_index),
                    ".sscore",
                ),
                sep="\t",
            )
            df = df.merge(
                pgs[["IID", "SCORE1_SUM"]].rename(columns={"SCORE1_SUM": f"PGS_{fold_index}"}),
                on="IID",
                how="left",
            )

        missing_cols = [col for col in [f"PGS_{fold_index}" for fold_index in range(N_FOLDS)] if df[col].isna().any()]
        if missing_cols:
            raise ValueError(
                f"Missing PGS values for gamma={gamma_value}, seed={seed_index}: {', '.join(missing_cols)}"
            )
        df.to_csv(output_path, sep="\t", index=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
