from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "outputs" / ".mplconfig"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(ROOT / "src"))

from predinterval_normalization.common import DEFAULT_OUTPUT_ROOT

MM_TO_INCH = 1.0 / 25.4
PANEL_WIDTH_MM = 45.0
FIG_HEIGHT_MM = 60.0
TITLE_FONTSIZE = 6.5
AXIS_LABEL_FONTSIZE = 6.0
TICK_FONTSIZE = 5.5
LEGEND_FONTSIZE = 5.0

METHOD_ORDER = ["PGS Only", "PGS + Covariates"]
GROUP_ORDER = ["quant_q_0", "marginal", "quant_q_9"]
GROUP_LABELS = {
    "quant_q_0": "Bottom decile",
    "marginal": "Overall",
    "quant_q_9": "Top decile",
}
STATE_PANEL_TITLES = {
    "Raw phenotype": "w/o normalization",
    "Normalized phenotype": "w/ normalization",
}
GROUP_PALETTE = {
    "PGS Only": {
        "quant_q_0": "#A9D2F8",
        "marginal": "#4A90E2",
        "quant_q_9": "#165AA7",
    },
    "PGS + Covariates": {
        "quant_q_0": "#D8B9F5",
        "marginal": "#9B6AE3",
        "quant_q_9": "#5B2C90",
    },
}


def _save(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def _plot_metric(df: pd.DataFrame, metric: str, title: str, reference_line: float, output_stem: Path) -> None:
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(2 * PANEL_WIDTH_MM * MM_TO_INCH, FIG_HEIGHT_MM * MM_TO_INCH),
        sharey=False,
    )
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.20, top=0.76, wspace=0.20)
    rng = np.random.default_rng(0)
    x_base = np.arange(len(METHOD_ORDER), dtype=float)
    offsets = {"quant_q_0": -0.22, "marginal": 0.0, "quant_q_9": 0.22}
    legend_handles: list[plt.Line2D] = []

    current = df[df["metric"] == metric].copy()
    for panel_idx, (ax, state) in enumerate(zip(axes, ["Raw phenotype", "Normalized phenotype"], strict=True)):
        panel = current[current["state"] == state]
        values = panel["value"].to_numpy(dtype=float)
        lower = min(values.min(), reference_line) - 0.03
        upper = max(values.max(), reference_line) + 0.03
        for method_idx, method in enumerate(METHOD_ORDER):
            subset = panel[panel["method"] == method]
            for group in GROUP_ORDER:
                group_values = subset.loc[subset["group"] == group, "value"].to_numpy(dtype=float)
                x_center = x_base[method_idx] + offsets[group]
                jitter = rng.uniform(-0.03, 0.03, size=group_values.shape[0])
                color = GROUP_PALETTE[method][group]
                ax.scatter(
                    np.full(group_values.shape[0], x_center) + jitter,
                    group_values,
                    s=10,
                    color=color,
                    alpha=0.28,
                    linewidths=0,
                    zorder=2,
                )
                ax.errorbar(
                    [x_center],
                    [float(group_values.mean())],
                    yerr=[float(group_values.std(ddof=1))],
                    fmt="o",
                    markersize=4.2,
                    color=color,
                    ecolor=color,
                    elinewidth=1.0,
                    capsize=2.0,
                    zorder=4,
                )

                if panel_idx == 0:
                    legend_handles.append(
                        plt.Line2D(
                            [0],
                            [0],
                            marker="o",
                            linestyle="None",
                            markersize=4.5,
                            markerfacecolor=color,
                            markeredgecolor=color,
                            label=f"{method} {GROUP_LABELS[group]}",
                        )
                    )
        ax.axhline(reference_line, color="black", linestyle="--", linewidth=0.9, zorder=1)
        ax.set_xticks(x_base)
        ax.set_xticklabels(METHOD_ORDER, fontsize=TICK_FONTSIZE)
        ax.set_title(STATE_PANEL_TITLES[state], fontsize=AXIS_LABEL_FONTSIZE)
        ax.set_ylim(lower, upper)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda value, _: f"{value:.0%}"))
        ax.grid(axis="y", linestyle=":", linewidth=0.5, color="#CCCCCC")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if panel_idx == 0:
            ax.set_ylabel("Coverage", fontsize=AXIS_LABEL_FONTSIZE)
        ax.tick_params(axis="both", labelsize=TICK_FONTSIZE, pad=1.5)

    unique_handles = []
    seen = set()
    for handle in legend_handles:
        if handle.get_label() in seen:
            continue
        seen.add(handle.get_label())
        unique_handles.append(handle)

    fig.legend(
        handles=unique_handles,
        frameon=False,
        fontsize=LEGEND_FONTSIZE,
        loc="upper center",
        ncol=3,
        bbox_to_anchor=(0.5, 0.93),
        handletextpad=0.5,
        columnspacing=1.1,
    )
    fig.suptitle(title, fontsize=TITLE_FONTSIZE, y=0.995)
    _save(fig, output_stem)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot contextual coverage for the raw vs normalized PredInterval comparison.")
    parser.add_argument(
        "--input-table",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT / "17-predinterval-comparison" / "coverage_long.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT / "18-figures",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input_table, sep="\t")
    _plot_metric(df, "coverage95", "95% contextual coverage", 0.95, args.output_dir / "coverage95")
    _plot_metric(df, "coverage50", "50% contextual coverage", 0.50, args.output_dir / "coverage50")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
