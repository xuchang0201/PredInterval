from __future__ import annotations

import numpy as np
import pandas as pd


def assign_extreme_quant_groups(values: pd.Series) -> pd.Series:
    """Assign labels compatible with q0/q9 summaries even when quant is discrete.

    For continuous quant with at least 10 unique values, this matches the original
    decile-based `pd.qcut(..., q=10)` behavior. For discrete inputs with fewer than
    10 unique values, the minimum quant value is labeled `0` and the maximum quant
    value is labeled `9`, which preserves the intended "bottom" vs "top" subgroup
    comparison used throughout the Figure 2 summaries.
    """

    if values.nunique(dropna=False) >= 10:
        return pd.qcut(values, q=10).cat.codes

    out = pd.Series(-1, index=values.index, dtype=int)
    min_value = values.min()
    max_value = values.max()
    out.loc[values == min_value] = 0
    out.loc[values == max_value] = 9
    if min_value == max_value:
        out.loc[:] = 0
    return out


def select_nonredundant_columns(df: pd.DataFrame, columns: list[str]) -> list[str]:
    """Keep covariates in order while dropping constant or duplicate columns.

    This is mainly needed for binary-quant experiments where `quant**2` (or
    `quant_sq`) can be identical to `quant`, which otherwise makes the variance
    design matrix singular for CalPred.
    """

    kept: list[str] = []
    for column in columns:
        series = df[column]
        if series.nunique(dropna=False) <= 1:
            continue
        if any(series.equals(df[existing]) for existing in kept):
            continue
        kept.append(column)
    return kept


def gamma_tag(value: float) -> str:
    return f"{value:g}"


def default_wcontext_gamma_values() -> list[float]:
    return [float(value) for value in range(0, 101, 5)]


def enumerate_gamma_values(gamma_values: list[float] | None) -> list[tuple[int, float]]:
    values = default_wcontext_gamma_values() if gamma_values is None else [float(value) for value in gamma_values]
    return list(enumerate(values))
