#!/usr/bin/env python3
"""Generate Part 1 report figures for Option 1 from saved CSV metrics.

Usage:
    uv run python reporting/generate_option1_part1_figures.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OPTION1_DIR = ROOT / "master_output_dir" / "option1"
METRICS_CSV = OPTION1_DIR / "option1_metrics.csv"
OUT_DIR = OPTION1_DIR / "python_report" / "figures_part1"

PART1_METHODS = ["FDM", "FVM-S1", "FVM-S2"]
COLORS = {
    "FDM": "#0072B2",
    "FVM-S1": "#E69F00",
    "FVM-S2": "#7A3E9D",
}
MARKERS = {
    "FDM": "o",
    "FVM-S1": "D",
    "FVM-S2": "^",
}


def configure_style() -> None:
    plt.style.use("default")
    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 300,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linestyle": "--",
            "font.size": 11,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "legend.fontsize": 10,
        }
    )


def load_metrics() -> pd.DataFrame:
    if not METRICS_CSV.exists():
        raise FileNotFoundError(f"Missing metrics CSV: {METRICS_CSV}")
    return pd.read_csv(METRICS_CSV)


def save_fig(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_metric(df: pd.DataFrame, y_col: str, ylabel: str, title: str, out_name: str, yscale: str = "linear") -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.8))

    for method in PART1_METHODS:
        subset = df[df["method"] == method].sort_values("states_per_electrode")
        if subset.empty:
            continue
        ax.plot(
            subset["states_per_electrode"],
            subset[y_col],
            color=COLORS[method],
            marker=MARKERS[method],
            linewidth=2.2,
            markersize=6.5,
            label=method,
        )

    if yscale == "log":
        ax.set_yscale("log")

    ax.set_xlabel("Radial nodes / control volumes per electrode")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False)
    save_fig(fig, OUT_DIR / out_name)


def main() -> None:
    configure_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    metrics = load_metrics()

    plot_metric(
        metrics,
        y_col="rmse_vs_exp_V",
        ylabel="RMSE vs experiment [V]",
        title="Grid Convergence: Voltage Error vs Spatial Resolution",
        out_name="grid_convergence.png",
        yscale="log",
    )
    plot_metric(
        metrics,
        y_col="runtime_s",
        ylabel="Runtime [s]",
        title="Runtime Comparison Across Spatial Resolution",
        out_name="runtime_comparison.png",
    )
    plot_metric(
        metrics,
        y_col="max_rel_inventory_drift",
        ylabel="Max relative inventory drift [-]",
        title="Lithium Conservation Across Spatial Resolution",
        out_name="li_conservation.png",
        yscale="log",
    )

    print(f"Part 1 figures exported to: {OUT_DIR}")


if __name__ == "__main__":
    main()
