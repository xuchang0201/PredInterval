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
    stage_dir,
    weight_path,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Simulate causal SNP weights for the original continuous wcontext pipeline.")
    add_common_generation_args(parser)
    args = parser.parse_args()

    weights_dir = stage_dir(args.output_root, "01-weights")
    weights_dir.mkdir(parents=True, exist_ok=True)

    snplist = pd.read_csv(args.geno_prefix.with_suffix(".bim"), sep="\t", header=None)
    freq = pd.read_csv(args.geno_prefix.with_suffix(".afreq"), sep="\t")
    n_snp = len(snplist)
    n_causal_snps = int(args.causal * n_snp)
    if n_causal_snps < 1:
        raise ValueError("--causal produced zero causal SNPs")

    for seed_index in tqdm(args.seed_indices, desc="simulate causal weights"):
        output_path = weight_path(args.output_root, args.h2, args.causal, seed_index)
        if output_path.exists() and not args.overwrite:
            continue

        causal_snp = snplist.sample(n_causal_snps, replace=False, random_state=seed_index).rename(columns={1: "ID"})
        np.random.seed(seed_index)
        causal_snp["beta"] = np.random.normal(0, np.sqrt(args.h2 / n_causal_snps), n_causal_snps)
        causal_snp = causal_snp.merge(freq[["ID", "ALT_FREQS"]], on="ID", how="inner")
        causal_snp["beta_freqAdj"] = causal_snp["beta"] / np.sqrt(
            2 * causal_snp["ALT_FREQS"] * (1 - causal_snp["ALT_FREQS"])
        )
        causal_snp.to_csv(output_path, sep="\t", index=False)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
