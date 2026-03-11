#!/usr/bin/env python3
"""Generate publication-style figures and a Markdown report from SPM CSV outputs.

Usage:
    uv run python reporting/generate_option1_report.py
    uv run python reporting/generate_option1_report.py --input-dir master_output_dir/option1

The script is intentionally data-driven:
1. It discovers result folders inside ``master_output_dir``.
2. It regenerates figures from any CSV summaries that are present.
3. It writes a Markdown report focused on observability and method tradeoffs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


METHOD_ORDER = ["FDM", "FVM-S0", "FVM-S1", "FVM-S2", "PADE2", "PADE3"]
METHOD_COLORS = {
    "FDM": "#0072B2",
    "FVM-S0": "#D55E00",
    "FVM-S1": "#E69F00",
    "FVM-S2": "#7B3294",
    "PADE2": "#56B4E9",
    "PADE3": "#A50F15",
}
METHOD_MARKERS = {
    "FDM": "o",
    "FVM-S0": "s",
    "FVM-S1": "D",
    "FVM-S2": "^",
    "PADE2": "v",
    "PADE3": "P",
}


@dataclass
class ExperimentData:
    name: str
    directory: Path
    metrics: pd.DataFrame | None
    obs_summary: pd.DataFrame | None
    obs_sweep: pd.DataFrame | None
    voltage_curves: pd.DataFrame | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("master_output_dir"),
        help="Path to a single result directory or the master output directory.",
    )
    parser.add_argument(
        "--report-subdir",
        default="python_report",
        help="Subdirectory created inside each experiment folder.",
    )
    return parser.parse_args()


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


def read_csv_if_exists(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    df = pd.read_csv(path)
    return df.replace([np.inf, -np.inf], np.nan)


def discover_experiments(input_dir: Path) -> list[ExperimentData]:
    input_dir = input_dir.resolve()

    if (input_dir / "option1_metrics.csv").exists() or any(input_dir.glob("*_metrics.csv")):
        directories = [input_dir]
    else:
        directories = sorted(
            {
                path.parent.resolve()
                for path in input_dir.rglob("*_metrics.csv")
                if path.is_file()
            }
        )

    experiments: list[ExperimentData] = []
    for directory in directories:
        metrics_file = next(iter(sorted(directory.glob("*_metrics.csv"))), None)
        obs_summary_file = next(
            iter(sorted(directory.glob("*_observability_summary.csv"))), None
        )
        obs_sweep_file = next(
            iter(sorted(directory.glob("*_observability_sweep_summary.csv"))), None
        )
        voltage_file = directory / "voltage_vs_time_curves.csv"

        experiments.append(
            ExperimentData(
                name=directory.name,
                directory=directory,
                metrics=read_csv_if_exists(metrics_file) if metrics_file else None,
                obs_summary=read_csv_if_exists(obs_summary_file) if obs_summary_file else None,
                obs_sweep=read_csv_if_exists(obs_sweep_file) if obs_sweep_file else None,
                voltage_curves=read_csv_if_exists(voltage_file),
            )
        )
    return experiments


def ordered_methods(df: pd.DataFrame | None) -> list[str]:
    if df is None or "method" not in df.columns:
        return METHOD_ORDER
    present = list(dict.fromkeys(df["method"].dropna().astype(str)))
    return [m for m in METHOD_ORDER if m in present] + [m for m in present if m not in METHOD_ORDER]


def style_for_method(method: str) -> dict[str, object]:
    return {
        "color": METHOD_COLORS.get(method, "#4D4D4D"),
        "marker": METHOD_MARKERS.get(method, "o"),
        "linewidth": 2.2,
        "markersize": 6.5,
    }


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def metric_line_plot(
    df: pd.DataFrame,
    y_column: str,
    title: str,
    ylabel: str,
    output_path: Path,
    yscale: str = "linear",
) -> bool:
    if df is None or y_column not in df.columns or "total_states" not in df.columns:
        return False

    methods = ordered_methods(df)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))

    plotted = False
    for method in methods:
        subset = df[df["method"] == method].sort_values("total_states")
        if subset.empty:
            continue
        ax.plot(subset["total_states"], subset[y_column], label=method, **style_for_method(method))
        plotted = True

    if not plotted:
        plt.close(fig)
        return False

    ax.set_xscale("log")
    if yscale == "log":
        ax.set_yscale("log")
    ax.set_title(title)
    ax.set_xlabel("Total states")
    ax.set_ylabel(ylabel)
    ax.legend(ncols=2 if len(methods) > 4 else 1, frameon=False)
    save_figure(fig, output_path)
    return True


def plot_voltage_curves(df: pd.DataFrame | None, output_path: Path) -> bool:
    if df is None or "time_s" not in df.columns:
        return False

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    x = df["time_s"]

    if "V_exp" in df.columns:
        ax.plot(x, df["V_exp"], color="black", linewidth=2.8, label="Experimental")

    simulation_columns = [col for col in df.columns if col not in {"time_s", "V_exp"}]
    for method in [m for m in METHOD_ORDER if m in simulation_columns] + [
        c for c in simulation_columns if c not in METHOD_ORDER
    ]:
        ax.plot(x, df[method], label=method, **style_for_method(method))

    ax.set_title("Voltage response against experiment")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Voltage [V]")
    ax.legend(ncols=3 if len(simulation_columns) >= 5 else 2, frameon=False)
    save_figure(fig, output_path)
    return True


def plot_accuracy_runtime(df: pd.DataFrame | None, output_path: Path) -> bool:
    if (
        df is None
        or "runtime_s" not in df.columns
        or "rmse_vs_dense_fdm_V" not in df.columns
        or "method" not in df.columns
    ):
        return False

    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    for method in ordered_methods(df):
        subset = df[df["method"] == method]
        if subset.empty:
            continue
        ax.scatter(
            subset["runtime_s"],
            subset["rmse_vs_dense_fdm_V"],
            s=72,
            label=method,
            color=METHOD_COLORS.get(method, "#4D4D4D"),
            alpha=0.88,
        )

    ax.set_title("Accuracy-runtime tradeoff")
    ax.set_xlabel("Runtime [s]")
    ax.set_ylabel("RMSE vs dense FDM [V]")
    ax.set_yscale("log")
    ax.legend(ncols=2, frameon=False)
    save_figure(fig, output_path)
    return True


def plot_observability_sweep(df: pd.DataFrame | None, output_path: Path) -> bool:
    required = {"method", "total_states", "min_rank_fraction"}
    if df is None or not required.issubset(df.columns):
        return False

    fig, axes = plt.subplots(3, 1, figsize=(7.2, 10.0), sharex=True)
    methods = ordered_methods(df)

    panels: list[tuple[str, str, str]] = [
        ("min_rank_fraction", "Minimum observable fraction", "linear"),
    ]

    if "min_sigma_min" in df.columns:
        panels.append(("min_sigma_min", r"Minimum $\sigma_{\min}(\mathcal{O})$", "log"))

    if "finite_condition_fraction" in df.columns:
        panels.append(("finite_condition_fraction", "Finite condition fraction", "linear"))
    elif "median_finite_condition_number" in df.columns:
        panels.append(
            ("median_finite_condition_number", r"Median finite cond$(\mathcal{O})$", "log")
        )

    while len(panels) < 3:
        panels.append(("min_rank_fraction", "Minimum observable fraction", "linear"))

    for ax, (column, ylabel, yscale) in zip(axes, panels[:3]):
        for method in methods:
            subset = df[df["method"] == method].sort_values("total_states")
            if subset.empty or column not in subset.columns:
                continue
            y = subset[column].copy()
            if yscale == "log":
                y = y.clip(lower=np.finfo(float).tiny)
                ax.set_yscale("log")
            ax.plot(subset["total_states"], y, label=method, **style_for_method(method))
        ax.set_xscale("log")
        ax.set_ylabel(ylabel)

    axes[0].set_title("Observability metrics across model order")
    axes[-1].set_xlabel("Total states")
    axes[0].legend(ncols=2, frameon=False)
    save_figure(fig, output_path)
    return True


def plot_observability_summary(df: pd.DataFrame | None, output_path: Path) -> bool:
    if df is None or "method" not in df.columns or "total_states" not in df.columns:
        return False

    summary = df.copy()
    if "min_rank" in summary.columns:
        summary["rank_fraction"] = summary["min_rank"] / summary["total_states"]
    elif "full_rank_fraction" in summary.columns:
        summary["rank_fraction"] = summary["full_rank_fraction"]
    else:
        return False

    summary = summary.sort_values(
        by="method", key=lambda s: s.map({m: i for i, m in enumerate(METHOD_ORDER)}).fillna(999)
    )

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8))
    colors = [METHOD_COLORS.get(m, "#4D4D4D") for m in summary["method"]]

    axes[0].bar(summary["method"], summary["rank_fraction"], color=colors)
    axes[0].set_ylim(0.0, min(1.05, max(1.0, summary["rank_fraction"].max() * 1.1)))
    axes[0].set_ylabel("Minimum rank / total states")
    axes[0].set_title("Representative observability rank")
    axes[0].tick_params(axis="x", rotation=25)

    if "min_sigma_min" in summary.columns:
        sigma = summary["min_sigma_min"].clip(lower=np.finfo(float).tiny)
        axes[1].bar(summary["method"], sigma, color=colors)
        axes[1].set_yscale("log")
        axes[1].set_ylabel(r"Minimum $\sigma_{\min}(\mathcal{O})$")
        axes[1].set_title("Smallest singular value")
        axes[1].tick_params(axis="x", rotation=25)
    else:
        axes[1].axis("off")

    save_figure(fig, output_path)
    return True


def best_row(df: pd.DataFrame, column: str, smallest: bool = True) -> pd.Series | None:
    if df is None or column not in df.columns:
        return None
    clean = df[np.isfinite(df[column])]
    if clean.empty:
        return None
    idx = clean[column].idxmin() if smallest else clean[column].idxmax()
    return clean.loc[idx]


def best_reduced_order_row(df: pd.DataFrame, max_states: int, column: str) -> pd.Series | None:
    if df is None or column not in df.columns or "total_states" not in df.columns:
        return None
    subset = df[df["total_states"] <= max_states]
    if subset.empty:
        return None
    return best_row(subset, column, smallest=True)


def representative_observability_summary(obs_summary: pd.DataFrame | None) -> list[str]:
    if obs_summary is None or not {"method", "min_rank", "total_states"}.issubset(obs_summary.columns):
        return []

    lines = []
    ranked = obs_summary.copy()
    ranked["rank_fraction"] = ranked["min_rank"] / ranked["total_states"]
    ranked = ranked.sort_values("rank_fraction", ascending=False)

    tol = 1e-12
    best_fraction = ranked["rank_fraction"].max()
    worst_fraction = ranked["rank_fraction"].min()
    best_rows = ranked[np.abs(ranked["rank_fraction"] - best_fraction) <= tol]
    worst_rows = ranked[np.abs(ranked["rank_fraction"] - worst_fraction) <= tol]

    best = best_rows.iloc[0]
    worst = worst_rows.iloc[0]
    best_methods = ", ".join(f"`{m}`" for m in best_rows["method"])
    worst_methods = ", ".join(f"`{m}`" for m in worst_rows["method"])

    lines.append(
        f"At the representative order, {best_methods} achieve the strongest observability with "
        f"minimum rank fraction {best['rank_fraction']:.3f} ({int(best['min_rank'])}/{int(best['total_states'])})."
    )
    lines.append(
        f"{worst_methods} are weakest with minimum rank fraction "
        f"{worst['rank_fraction']:.3f} ({int(worst['min_rank'])}/{int(worst['total_states'])})."
    )
    return lines


def observability_sweep_summary(obs_sweep: pd.DataFrame | None) -> list[str]:
    if obs_sweep is None or not {"method", "total_states", "min_rank_fraction"}.issubset(obs_sweep.columns):
        return []

    lines = []
    for method in ordered_methods(obs_sweep):
        subset = obs_sweep[obs_sweep["method"] == method].sort_values("total_states")
        if subset.empty:
            continue
        first = subset.iloc[0]
        last = subset.iloc[-1]
        lines.append(
            f"`{method}` decreases from a minimum observable fraction of "
            f"{first['min_rank_fraction']:.3f} at {int(first['total_states'])} states to "
            f"{last['min_rank_fraction']:.3f} at {int(last['total_states'])} states."
        )
    if "finite_condition_fraction" in obs_sweep.columns:
        finite = obs_sweep["finite_condition_fraction"].fillna(0.0)
        if finite.max() == 0.0:
            lines.append(
                "No sampled case produced a finite observability condition number, which means the "
                "condition number is not a useful discriminator here and rank-based metrics should "
                "carry the interpretation."
            )
    return lines


def build_report_text(experiment: ExperimentData, figure_paths: dict[str, Path]) -> str:
    metrics = experiment.metrics
    obs_summary = experiment.obs_summary
    obs_sweep = experiment.obs_sweep

    best_exp = best_row(metrics, "rmse_vs_exp_V", smallest=True) if metrics is not None else None
    best_dense = None
    if metrics is not None and "rmse_vs_dense_fdm_V" in metrics.columns:
        dense_candidates = metrics[metrics["rmse_vs_dense_fdm_V"] > 0]
        if not dense_candidates.empty:
            best_dense = best_row(dense_candidates, "rmse_vs_dense_fdm_V", smallest=True)
    best_runtime = best_row(metrics, "runtime_s", smallest=True) if metrics is not None else None
    best_reduced = (
        best_reduced_order_row(metrics, max_states=30, column="rmse_vs_exp_V")
        if metrics is not None
        else None
    )

    obs_lines = representative_observability_summary(obs_summary)
    obs_lines.extend(observability_sweep_summary(obs_sweep))

    key_findings: list[str] = []
    if best_exp is not None:
        key_findings.append(
            f"The smallest voltage RMSE against experiment is achieved by `{best_exp['method']}` "
            f"at {int(best_exp['total_states'])} states with RMSE = {best_exp['rmse_vs_exp_V']:.4e} V."
        )
    if best_dense is not None:
        key_findings.append(
            f"The closest model to the dense FDM reference is `{best_dense['method']}` at "
            f"{int(best_dense['total_states'])} states with RMSE = {best_dense['rmse_vs_dense_fdm_V']:.4e} V."
        )
    if best_reduced is not None:
        key_findings.append(
            f"Among reduced-order models with at most 30 states, `{best_reduced['method']}` gives the "
            f"best match to experiment with RMSE = {best_reduced['rmse_vs_exp_V']:.4e} V."
        )
    if best_runtime is not None:
        key_findings.append(
            f"The fastest simulated configuration is `{best_runtime['method']}` at "
            f"{int(best_runtime['total_states'])} states with runtime = {best_runtime['runtime_s']:.4f} s, "
            "but speed must be interpreted together with error and observability."
        )

    figure_lines = []
    for label, path in figure_paths.items():
        figure_lines.append(f"- {label}: `{path.as_posix()}`")

    observability_bullets = "\n".join(f"- {line}" for line in obs_lines) if obs_lines else "- Observability CSV summaries were not available."
    findings_bullets = "\n".join(f"- {line}" for line in key_findings) if key_findings else "- No quantitative metrics were available."

    return f"""# Final Project Report

## SPM Discretization Comparison with Observability Analysis

### Introduction
This report analyzes the `Option 1` final-project task: comparing discretization methods for the Single Particle Model (SPM) of a lithium-ion cell. The main goal is not only to compare accuracy and runtime, but to understand how the choice of spatial discretization changes the **ability to infer internal states from terminal voltage**, i.e. observability. That question matters directly for estimator design because a battery management system relies on voltage and current to reconstruct internal lithium concentration states that cannot be measured directly.

### Background
The SPM models solid-phase lithium diffusion in each electrode particle and maps the resulting surface concentrations to terminal voltage through open-circuit potential (OCP) and overpotential terms. After spatial discretization, the model has the linear state dynamics

`x_dot = A x + B I`

with a nonlinear voltage output `V = h(x, I)`. Around a trajectory sample `k`, the MATLAB workflow linearizes the voltage with respect to the states, forming a row vector `C_k = dh/dx`. Observability is then quantified through the classical observability matrix

`O_k = [C_k; C_k A; C_k A^2; ...; C_k A^(n-1)]`

where `n` is the total number of states. The project data summarize observability using three complementary diagnostics:

1. **Rank of `O_k`**: measures how many state directions are theoretically observable.
2. **Smallest singular value `sigma_min(O_k)`**: measures how close the matrix is to losing rank numerically.
3. **Condition number `cond(O_k)`**: measures how sensitive state reconstruction is to noise and numerical perturbations.

The MATLAB driver also records the OCP slope difference `|dU_p/dtheta_p - dU_n/dtheta_n|` as a physics-based proxy for output sensitivity. A larger slope difference helps, but the final observability outcome still depends on the full dynamic pair `(A, C_k)`.

### Methodology
The comparison includes FDM, three FVM surface reconstructions (`FVM-S0`, `FVM-S1`, `FVM-S2`), and two Padé reduced-order models (`PADE2`, `PADE3`). The figures and statistics in this report were generated from CSV exports stored in `{experiment.directory.as_posix()}`:

{chr(10).join(figure_lines)}

The Python reporting pipeline regenerates figures directly from updated CSV outputs, so the same report structure can be rerun after additional experiments without rewriting plotting code manually.

### Results
#### Overall Numerical Performance
{findings_bullets}

#### Observability Results
{observability_bullets}

The representative-order observability summary is especially informative. At comparable model orders, the FVM variants retain a larger observable fraction than FDM, while the Padé models have the smallest observable fraction despite their low state counts. This means reduced state dimension alone is not enough; a reduced model must preserve voltage-sensitive state directions if it is to remain useful for estimation.

### Discussion
The key tradeoff is that increasing model order improves voltage fidelity up to a point, but it also expands the number of internal concentration states faster than the terminal-voltage measurement can distinguish them. This appears clearly in the sweep data: rank grows sublinearly compared with total states, so the observable fraction falls as the discretization is refined.

For this dataset, the FVM family offers the strongest compromise. The FVM schemes achieve experimental RMSE comparable to dense FDM at much smaller runtime, and at the representative order they preserve a larger observable fraction than FDM. Among the FVM options, `FVM-S2` is particularly attractive because it is both accurate and fast while staying competitive in the observability summary.

The Padé models show the opposite behavior. They are compact and fast, but the voltage error is orders of magnitude larger than the grid-based methods and the minimum observable fraction is only about one-half. Their larger OCP slope proxy does not translate into good overall observability or accuracy, confirming that scalar sensitivity alone is not a sufficient design criterion.

Another important result is that the sampled observability matrices never achieve a finite condition number in the exported summaries. That indicates severe numerical ill-conditioning across all tested methods and orders. In practice, this means rank fraction and smallest singular value are more reliable diagnostics than condition number for comparing these cases.

### Limitations
- The exported CSVs do not include the full SOC-resolved observability trajectories, so this report focuses on the summarized observability metrics that were saved from MATLAB.
- The analysis is based on the HPPC experiment and should be repeated for UDDS or other drive cycles before making broad claims about estimator performance.
- Observability is evaluated from a local linearization of the nonlinear voltage map, so it should be interpreted as a practical local metric rather than a complete nonlinear observability proof.

### Conclusions
The results support three main conclusions. First, **observability degrades with increasing state dimension for every discretization method tested**. Second, **FVM schemes provide the best balance of voltage accuracy, runtime, and retained observability**, making them the strongest candidates for reduced-order estimation-oriented SPM models. Third, **very compact Padé models are computationally attractive but too inaccurate and weakly observable to recommend for this final-project objective**.

For estimator development, the most defensible next step is to build EKF or AEKF designs on top of the better-performing FVM reduced-order models, especially the `FVM-S1` and `FVM-S2` variants at moderate state counts.
"""


def write_report(report_text: str, output_path: Path) -> None:
    output_path.write_text(report_text, encoding="utf-8")


def generate_for_experiment(experiment: ExperimentData, report_subdir: str) -> None:
    report_dir = ensure_dir(experiment.directory / report_subdir)
    figure_dir = ensure_dir(report_dir / "figures")

    figure_paths: dict[str, Path] = {}

    figure_specs = [
        (
            "Voltage validation",
            figure_dir / "voltage_vs_time.png",
            lambda output: plot_voltage_curves(experiment.voltage_curves, output),
        ),
        (
            "RMSE vs experiment",
            figure_dir / "rmse_vs_experiment.png",
            lambda output: metric_line_plot(
                experiment.metrics,
                "rmse_vs_exp_V",
                "Model accuracy against experiment",
                "RMSE vs experiment [V]",
                output,
                yscale="log",
            ),
        ),
        (
            "RMSE vs dense FDM",
            figure_dir / "rmse_vs_dense_fdm.png",
            lambda output: metric_line_plot(
                experiment.metrics,
                "rmse_vs_dense_fdm_V",
                "Convergence to dense FDM reference",
                "RMSE vs dense FDM [V]",
                output,
                yscale="log",
            ),
        ),
        (
            "Runtime scaling",
            figure_dir / "runtime_scaling.png",
            lambda output: metric_line_plot(
                experiment.metrics,
                "runtime_s",
                "Runtime scaling with model order",
                "Runtime [s]",
                output,
            ),
        ),
        (
            "Lithium conservation",
            figure_dir / "lithium_conservation.png",
            lambda output: metric_line_plot(
                experiment.metrics,
                "max_rel_inventory_drift",
                "Lithium inventory drift",
                "Max relative inventory drift [-]",
                output,
                yscale="log",
            ),
        ),
        (
            "Accuracy-runtime tradeoff",
            figure_dir / "accuracy_vs_runtime.png",
            lambda output: plot_accuracy_runtime(experiment.metrics, output),
        ),
        (
            "Observability sweep",
            figure_dir / "observability_summary_vs_states.png",
            lambda output: plot_observability_sweep(experiment.obs_sweep, output),
        ),
        (
            "Representative observability",
            figure_dir / "representative_observability.png",
            lambda output: plot_observability_summary(experiment.obs_summary, output),
        ),
    ]

    for label, output_path, builder in figure_specs:
        created = builder(output_path)
        if created:
            figure_paths[label] = output_path.relative_to(report_dir)

    report_text = build_report_text(experiment, figure_paths)
    write_report(report_text, report_dir / "final_project_report.md")


def main() -> None:
    args = parse_args()
    configure_style()

    experiments = discover_experiments(args.input_dir)
    if not experiments:
        raise SystemExit(f"No experiment folders with '*_metrics.csv' found under {args.input_dir}")

    for experiment in experiments:
        generate_for_experiment(experiment, report_subdir=args.report_subdir)
        print(f"Generated report assets for {experiment.name}: {experiment.directory / args.report_subdir}")


if __name__ == "__main__":
    main()
