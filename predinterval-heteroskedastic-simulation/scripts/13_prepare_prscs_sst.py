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
    append_suffix,
    rebuild_gwas_prefix,
    rebuild_sst_path,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert stage-12 PLINK GWAS outputs into PRS-CS summary-stat tables.")
    add_common_generation_args(parser)
    args = parser.parse_args()

    sst_dir = args.output_root / "13-sst"
    sst_dir.mkdir(parents=True, exist_ok=True)

    tasks = [
        (gamma_value, seed_index, fold_index)
        for gamma_value in args.gamma_values
        for seed_index in args.seed_indices
        for fold_index in range(N_FOLDS)
    ]
    for gamma_value, seed_index, fold_index in tqdm(tasks, desc="prepare PRS-CS sst"):
        output_path = rebuild_sst_path(args.output_root, args.h2, args.causal, gamma_value, seed_index, fold_index)
        if output_path.exists() and not args.overwrite:
            continue

        glm_path = append_suffix(
            rebuild_gwas_prefix(args.output_root, args.h2, args.causal, gamma_value, seed_index, fold_index),
            ".y.glm.linear",
        )
        sst = pd.read_csv(glm_path, sep="\t")
        sst = sst.drop(columns=["A1"], errors="ignore")
        sst = sst.rename(columns={"ID": "SNP", "ALT": "A1", "REF": "A2"})
        sst[["SNP", "A1", "A2", "BETA", "SE"]].to_csv(output_path, sep="\t", index=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
