from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd
from tqdm import tqdm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from predinterval_normalization.common import (
    FOLD_SIZE,
    N_FOLDS,
    TRAIN_SIZE,
    add_common_generation_args,
    cv_covar_path,
    cv_pheno_path,
    trait_path,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Create leave-one-fold-out GWAS phenotype and covariate tables from simulated continuous wcontext traits.")
    add_common_generation_args(parser)
    args = parser.parse_args()

    (args.output_root / "04-cv" / "pheno").mkdir(parents=True, exist_ok=True)
    (args.output_root / "04-cv" / "covar").mkdir(parents=True, exist_ok=True)

    tasks = [(gamma_value, seed_index) for gamma_value in args.gamma_values for seed_index in args.seed_indices]
    for gamma_value, seed_index in tqdm(tasks, desc="cv splits"):
        df = pd.read_csv(trait_path(args.output_root, args.h2, args.causal, gamma_value, seed_index), sep="\t")
        df = df.iloc[:TRAIN_SIZE].copy()
        covar = df[["FID", "IID", "quant"]]
        pheno = df[["FID", "IID", "y"]]

        for fold_index in range(N_FOLDS):
            fold_slice = slice(fold_index * FOLD_SIZE, (fold_index + 1) * FOLD_SIZE)
            held_out = df.iloc[fold_slice].IID
            current_covar = covar[~covar.IID.isin(held_out)]
            current_pheno = pheno[~pheno.IID.isin(held_out)]

            pheno_path = cv_pheno_path(args.output_root, args.h2, args.causal, gamma_value, seed_index, fold_index)
            covar_path = cv_covar_path(args.output_root, args.h2, args.causal, gamma_value, seed_index, fold_index)
            if pheno_path.exists() and covar_path.exists() and not args.overwrite:
                continue
            current_pheno.to_csv(pheno_path, sep="\t", index=False)
            current_covar.to_csv(covar_path, sep="\t", index=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
