from __future__ import annotations

from pathlib import Path
import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "master_output_dir" / "option1"
FIG_DIR = OUT_DIR / "slide_deck_figures_python"
MAT_PATH = OUT_DIR / "option1_metrics_and_results.mat"
METRICS_CSV = OUT_DIR / "option1_metrics.csv"
SWEEP_CSV = OUT_DIR / "option1_observability_sweep_summary.csv"


THETA_N_0 = 0.027935145151956345
THETA_P_0 = 0.9066422666349141
THETA_N_100 = 0.9219907588763882
THETA_P_100 = 0.267963982283029


COLORS = {
    "FDM": "#0072B2",
    "FVM-S2": "#7A3E9D",
    "PADE2": "#4DBBD5",
    "PADE3": "#8B1E3F",
}

MARKERS = {
    "FDM": "o",
    "FVM-S2": "^",
    "PADE2": "v",
    "PADE3": ">",
}

FOCUS_METHODS = ["FDM", "FVM-S2", "PADE2", "PADE3"]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 450,
            "font.family": "DejaVu Sans",
            "font.size": 16,
            "axes.titlesize": 20,
            "axes.titleweight": "bold",
            "axes.labelsize": 18,
            "axes.labelweight": "bold",
            "axes.linewidth": 1.4,
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
            "legend.fontsize": 14,
            "legend.frameon": False,
            "grid.color": "#D0D0D0",
            "grid.alpha": 0.45,
            "grid.linewidth": 0.8,
        }
    )


def clean_axes(ax: plt.Axes) -> None:
    ax.grid(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=6, width=1.2)


def u_n(theta: np.ndarray) -> np.ndarray:
    return (
        1.9793 * np.exp(-39.3631 * theta)
        + 0.2482
        - 0.0909 * np.tanh(29.8538 * (theta - 0.1234))
        - 0.04478 * np.tanh(14.9159 * (theta - 0.2769))
        - 0.0205 * np.tanh(30.4444 * (theta - 0.6103))
    )


def u_p(theta: np.ndarray) -> np.ndarray:
    return (
        -0.8090 * theta
        + 4.4875
        - 0.0428 * np.tanh(18.5138 * (theta - 0.5542))
        - 17.7326 * np.tanh(15.7890 * (theta - 0.3117))
        + 17.5842 * np.tanh(15.9308 * (theta - 0.3120))
    )


def theta_from_soc(soc: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    theta_n = THETA_N_0 + soc * (THETA_N_100 - THETA_N_0)
    theta_p = THETA_P_0 + soc * (THETA_P_100 - THETA_P_0)
    return theta_n, theta_p


def load_observability_results():
    mat = loadmat(MAT_PATH, squeeze_me=True, struct_as_record=False)
    obs = mat["observability_results"]
    obs_list = np.atleast_1d(obs)
    by_method = {}
    for item in obs_list:
        by_method[str(item.method)] = item
    return by_method


def save_single_figure(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_slide1_figures() -> None:
    soc = np.linspace(0.0, 1.0, 500)
    theta_n, theta_p = theta_from_soc(soc)
    un = u_n(theta_n)
    up = u_p(theta_p)
    voc = up - un

    dun = np.gradient(un, soc)
    dup = np.gradient(up, soc)
    slope_metric = np.abs(dup - dun)

    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    clean_axes(ax)
    ax.plot(soc, voc, color="black", lw=3.2, label=r"Cell OCV = $U_p - U_n$")
    ax.plot(soc, up, color="#B22222", lw=2.8, label=r"Cathode OCP, $U_p$")
    ax.plot(soc, un, color="#1F4E79", lw=2.8, label=r"Anode OCP, $U_n$")
    ax.set_xlabel("SOC [-]")
    ax.set_ylabel("Voltage [V]")
    ax.set_title("OCP Curves Across SOC")
    ax.legend(loc="best")
    save_single_figure(fig, FIG_DIR / "slide1_fig1_ocp_curves.png")

    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    clean_axes(ax)
    ax.plot(soc, np.abs(dup), color="#B22222", lw=2.8, label=r"$|dU_p/dSOC|$")
    ax.plot(soc, np.abs(dun), color="#1F4E79", lw=2.8, label=r"$|dU_n/dSOC|$")
    ax.plot(soc, slope_metric, color="black", lw=3.2, label="Net OCP slope metric")
    ax.set_xlabel("SOC [-]")
    ax.set_ylabel("Slope magnitude")
    ax.set_title("Voltage Sensitivity vs SOC")
    ax.legend(loc="best")
    save_single_figure(fig, FIG_DIR / "slide1_fig2_ocp_slope_metric.png")


def make_slide2_figures(obs_by_method) -> None:
    soc_all = []
    for method in FOCUS_METHODS:
        obs = obs_by_method[method]
        soc_all.append(np.asarray(obs.soc_ref).reshape(-1))
    soc_all = np.concatenate(soc_all)
    soc_min = float(np.min(soc_all))
    soc_max = float(np.max(soc_all))
    soc_pad = max(0.005, 0.08 * (soc_max - soc_min))

    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for method in FOCUS_METHODS:
        obs = obs_by_method[method]
        soc = np.asarray(obs.soc_ref).reshape(-1)
        rank = np.asarray(obs.rank).reshape(-1)
        ax.plot(
            soc,
            rank,
            lw=2.8,
            color=COLORS[method],
            label=method,
        )
    clean_axes(ax)
    ax.set_xlim(min(1.0, soc_max + soc_pad), max(0.0, soc_min - soc_pad))
    ax.invert_xaxis()
    ax.set_xlabel("SOC [-]")
    ax.set_ylabel("Observability rank")
    ax.set_title("Rank of O vs SOC")
    ax.legend(loc="best")
    save_single_figure(fig, FIG_DIR / "slide2_fig1_rank_vs_soc.png")

    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for method in FOCUS_METHODS:
        obs = obs_by_method[method]
        soc = np.asarray(obs.soc_ref).reshape(-1)
        sigma_min = np.asarray(obs.sigma_min).reshape(-1)
        sigma_plot = np.clip(sigma_min, 1e-30, None)
        ax.semilogy(
            soc,
            sigma_plot,
            lw=2.8,
            color=COLORS[method],
            label=method,
        )
    clean_axes(ax)
    ax.set_xlim(min(1.0, soc_max + soc_pad), max(0.0, soc_min - soc_pad))
    ax.invert_xaxis()
    ax.set_xlabel("SOC [-]")
    ax.set_ylabel(r"Smallest singular value, $\sigma_{\min}(O)$")
    ax.set_title("Weakest Observable Direction vs SOC")
    ax.legend(loc="best")
    save_single_figure(fig, FIG_DIR / "slide2_fig2_sigma_min_vs_soc.png")


def make_slide3_figures(sweep_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for method in FOCUS_METHODS:
        df = sweep_df[sweep_df["method"] == method].sort_values("total_states")
        ax.plot(
            df["total_states"],
            df["min_rank_fraction"],
            lw=2.6,
            marker=MARKERS[method],
            ms=7.5,
            color=COLORS[method],
            label=method,
        )
    clean_axes(ax)
    ax.set_xscale("log")
    ax.set_xlabel("Total states")
    ax.set_ylabel("Minimum rank / n")
    ax.set_title("Normalized Rank vs Model Order")
    ax.legend(loc="best")
    save_single_figure(fig, FIG_DIR / "slide3_fig1_rank_fraction_vs_states.png")

    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for method in FOCUS_METHODS:
        df = sweep_df[sweep_df["method"] == method].sort_values("total_states")
        sigma_plot = np.clip(df["min_sigma_min"].to_numpy(), 1e-30, None)
        ax.semilogy(
            df["total_states"],
            sigma_plot,
            lw=2.6,
            marker=MARKERS[method],
            ms=7.5,
            color=COLORS[method],
            label=method,
        )
    clean_axes(ax)
    ax.set_xscale("log")
    ax.set_xlabel("Total states")
    ax.set_ylabel(r"Minimum $\sigma_{\min}(O)$")
    ax.set_title("Weakest Observable Direction vs Model Order")
    ax.legend(loc="best")
    save_single_figure(fig, FIG_DIR / "slide3_fig2_sigma_min_vs_states.png")


def make_slide4_figures(metrics_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for method in FOCUS_METHODS:
        df = metrics_df[metrics_df["method"] == method].sort_values("total_states")
        ax.plot(
            df["total_states"],
            1e3 * df["rmse_vs_exp_V"],
            lw=2.6,
            marker=MARKERS[method],
            ms=7.5,
            color=COLORS[method],
            label=method,
        )
    clean_axes(ax)
    ax.set_xscale("log")
    ax.set_xlabel("Total states")
    ax.set_ylabel("Voltage RMSE vs experiment [mV]")
    ax.set_title("Accuracy vs Model Order")
    ax.legend(loc="best")
    save_single_figure(fig, FIG_DIR / "slide4_fig1_rmse_vs_states.png")

    fig, ax = plt.subplots(figsize=(8.5, 5.8))
    for method in FOCUS_METHODS:
        df = metrics_df[metrics_df["method"] == method].sort_values("obs_min_rank_fraction")
        ax.plot(
            df["obs_min_rank_fraction"],
            1e3 * df["rmse_vs_exp_V"],
            lw=2.6,
            marker=MARKERS[method],
            ms=7.5,
            color=COLORS[method],
            label=method,
        )
    clean_axes(ax)
    ax.set_xlabel("Minimum rank / n")
    ax.set_ylabel("Voltage RMSE vs experiment [mV]")
    ax.set_title("Accuracy vs Observability Strength")
    ax.legend(loc="best")
    save_single_figure(fig, FIG_DIR / "slide4_fig2_rmse_vs_rank_fraction.png")


def write_slide_content() -> None:
    text = """# Observability Slide Deck Content

## Slide 1: Why Observability Changes With SOC
- Terminal voltage sensitivity is SOC-dependent because both OCP and kinetic terms vary with stoichiometry.
- Flat OCP regions reduce the information content of voltage measurements.
- This motivates trajectory-dependent observability analysis rather than a single fixed linear test.

Figures:
- `slide1_fig1_ocp_curves.png`
- `slide1_fig2_ocp_slope_metric.png`

Speaker notes:
- Introduce the link between battery chemistry and estimation.
- Emphasize that voltage carries more state information where the OCP curve is steep.

## Slide 2: Local Observability Along the HPPC Trajectory
- At each time step, the local output Jacobian changes, so the observability matrix changes with SOC.
- The rank stays limited for all methods, showing that voltage alone cannot fully reconstruct all diffusion states.
- The smallest singular value collapses in weakly informative SOC regions, indicating fragile reconstruction.

Figures:
- `slide2_fig1_rank_vs_soc.png`
- `slide2_fig2_sigma_min_vs_soc.png`

Speaker notes:
- Explain that rank is the structural test, while the smallest singular value measures practical strength.
- Highlight the same weak-observability regions across methods.

## Slide 3: Observability vs Model Order
- Increasing model order improves spatial fidelity but does not improve reconstructability from a single voltage output.
- The normalized observability rank generally decreases as the number of states increases.
- FDM and FVM-S2 follow similar trends because they discretize the same diffusion physics.

Figures:
- `slide3_fig1_rank_fraction_vs_states.png`
- `slide3_fig2_sigma_min_vs_states.png`

Speaker notes:
- Stress the tradeoff: more states means a richer model, but also a harder estimation problem.
- Use this slide to justify why high-order physics models are not automatically better for observers.

## Slide 4: Practical Tradeoff With Accuracy
- FDM and FVM-S2 maintain low voltage RMSE while Padé is less accurate, even after fixing the implementation.
- Padé remains computationally attractive because of its low order, but not because it improves observability.
- Main takeaway: OCP slope limits practical observability, and discretization choice mainly shifts the accuracy-cost tradeoff.

Figures:
- `slide4_fig1_rmse_vs_states.png`
- `slide4_fig2_rmse_vs_rank_fraction.png`

Speaker notes:
- Mention the Padé corrections briefly: the model is now physically consistent, but still less accurate.
- End with the central message: reduced order helps computation, but it does not remove the fundamental voltage observability limitation.
"""
    (FIG_DIR / "slide_contents.md").write_text(text)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    setup_style()

    metrics_df = pd.read_csv(METRICS_CSV)
    sweep_df = pd.read_csv(SWEEP_CSV)
    obs_by_method = load_observability_results()

    make_slide1_figures()
    make_slide2_figures(obs_by_method)
    make_slide3_figures(sweep_df)
    make_slide4_figures(metrics_df)
    write_slide_content()

    print(f"Python slide figures exported to: {FIG_DIR}")


if __name__ == "__main__":
    main()
