# PredInterval Normalization Pipeline

This repository reproduces the `h2 = 0.05` simulation pipeline used to evaluate PredInterval before and after phenotype normalization under heteroskedastic outcomes.

It includes:

- simulation of SNP effects, genetic values, and phenotypes;
- raw PGS construction with `GWAS -> PRS-CS -> PLINK score`;
- direct phenotype normalization for the `h2 = 0.05` setting;
- rebuilt PGS construction after normalization;
- two PredInterval implementations:
  - `PGS Only`
  - `PGS + Covariates`
- contextual coverage summaries and coverage plots.

## External Requirements

- Python 3.12+
- `uv`
- `plink2`
- PRS-CS reference files

Install Python dependencies with:

```bash
uv sync
```

## Input Files

Place the genotype files and PRS-CS reference files under `data/reference/`:

```text
data/reference/
  geno_82k_hm3_chr1.bed
  geno_82k_hm3_chr1.bim
  geno_82k_hm3_chr1.fam
  geno_82k_hm3_chr1.afreq
  prscs_ref/
    ldblk_1kg_eur/
      ...
```

The genotype files used in this pipeline come from https://doi.org/10.5281/zenodo.18840380, which was released with the CalPred paper, "CalPred yields calibrated intervals for polygenic risk prediction".

## Scope

This release corresponds to the `h2 = 0.05` experiment with direct normalization of `y` using the contextual variable. It does not reproduce the separate `h2 = 0.5` residual-based normalization analysis, which uses a different stage-10 procedure.

## Default Setting

The default run corresponds to:

- `h2 = 0.05`
- `causal = 1.0`
- `gamma = 100`
- `seed = 0..9`

## Full Run

Run the complete pipeline with:

```bash
bash scripts/run_full_pipeline.sh
```

By default, outputs are written to `outputs/run/`.

## Pipeline Stages

Raw phenotype pipeline:

1. `scripts/01_simulate_gv_weights.py`
2. `scripts/02_plink_gv.sh`
3. `scripts/03_generate_traits.py`
4. `scripts/04_make_cv_inputs_raw.py`
5. `scripts/05_plink_gwas_raw.sh`
6. `scripts/06_prepare_prscs_sst_raw.py`
7. `scripts/07_run_prscs_raw.sh`
8. `scripts/08_plink_score_raw.sh`
9. `scripts/09_assemble_sim_df_pgs_raw.py`

Normalization and rebuilt pipeline:

1. `scripts/10_adjust_phenotype_direct.py`
2. `scripts/11_make_cv_inputs_adjusted.py`
3. `scripts/12_plink_gwas_adjusted.sh`
4. `scripts/13_prepare_prscs_sst.py`
5. `scripts/14_run_prscs.sh`
6. `scripts/15_plink_score_adjusted.sh`
7. `scripts/16_assemble_sim_df_pgs_adjusted.py`

PredInterval and coverage outputs:

1. `scripts/17_run_predinterval_comparison.py`
2. `scripts/18_plot_coverage.py`

## Manual Execution

If you want to run stages one by one:

```bash
uv run python scripts/01_simulate_gv_weights.py
bash scripts/02_plink_gv.sh
uv run python scripts/03_generate_traits.py
uv run python scripts/04_make_cv_inputs_raw.py
bash scripts/05_plink_gwas_raw.sh
uv run python scripts/06_prepare_prscs_sst_raw.py
bash scripts/07_run_prscs_raw.sh
bash scripts/08_plink_score_raw.sh
uv run python scripts/09_assemble_sim_df_pgs_raw.py
uv run python scripts/10_adjust_phenotype_direct.py
uv run python scripts/11_make_cv_inputs_adjusted.py
bash scripts/12_plink_gwas_adjusted.sh
uv run python scripts/13_prepare_prscs_sst.py
bash scripts/14_run_prscs.sh
bash scripts/15_plink_score_adjusted.sh
uv run python scripts/16_assemble_sim_df_pgs_adjusted.py
uv run python scripts/17_run_predinterval_comparison.py
uv run python scripts/18_plot_coverage.py
```

## Main Outputs

- `outputs/run/09-sim_df_pgs/`
  Raw `sim_df_PGS` files.
- `outputs/run/10-adjusted-traits/`
  Normalized phenotypes.
- `outputs/run/16-sim_df_pgs/`
  Rebuilt `sim_df_PGS` files after normalization.
- `outputs/run/17-predinterval-comparison/coverage_long.tsv`
  Per-seed coverage and interval-length results.
- `outputs/run/17-predinterval-comparison/method_summary.tsv`
  Mean and standard deviation across seeds.
- `outputs/run/18-figures/coverage95.pdf`
- `outputs/run/18-figures/coverage50.pdf`

## Environment Variables

These optional variables can be used to override defaults:

- `OUTPUT_ROOT`
- `GENO_PREFIX`
- `PLINK_BIN`
- `PYTHON_BIN`
- `PRSCS_PY`
- `PRSCS_REF_DIR`
- `CAUSAL`
- `GAMMAS`
- `SEED_START`
- `SEED_END`
- `THREADS`
