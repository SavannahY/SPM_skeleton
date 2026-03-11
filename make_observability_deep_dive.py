from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parent
OUT_DIR = ROOT / "master_output_dir" / "option1"
FIG_DIR = OUT_DIR / "observability_deep_dive"
MAT_PATH = OUT_DIR / "option1_metrics_and_results.mat"
DEEP_CSV = OUT_DIR / "observability_deep_dive_timeseries.csv"
METRICS_CSV = OUT_DIR / "option1_metrics.csv"
SWEEP_CSV = OUT_DIR / "option1_observability_sweep_summary.csv"

COLORS = {
    "FDM": "#0072B2",
    "FVM-S0": "#D55E00",
    "FVM-S1": "#E69F00",
    "FVM-S2": "#7A3E9D",
    "Nonlinear SPM-Padé 2": "#4DBBD5",
    "Nonlinear SPM-Padé 3": "#8B1E3F",
    "Local Linear Padé-ECM": "#009E73",
}

MARKERS = {
    "FDM": "o",
    "FVM-S0": "s",
    "FVM-S1": "d",
    "FVM-S2": "^",
    "Nonlinear SPM-Padé 2": "v",
    "Nonlinear SPM-Padé 3": ">",
    "Local Linear Padé-ECM": "p",
}

SPM_METHODS = ["FDM", "FVM-S0", "FVM-S1", "FVM-S2"]
PADE_METHODS = ["Local Linear Padé-ECM", "Nonlinear SPM-Padé 2", "Nonlinear SPM-Padé 3"]
SUMMARY_METHODS = ["FDM", "FVM-S2", "Local Linear Padé-ECM", "Nonlinear SPM-Padé 2", "Nonlinear SPM-Padé 3"]

THETA_N_0 = 0.027935145151956345
THETA_P_0 = 0.9066422666349141
THETA_N_100 = 0.9219907588763882
THETA_P_100 = 0.267963982283029


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 500,
            "font.family": "DejaVu Sans",
            "font.size": 15,
            "axes.titlesize": 19,
            "axes.titleweight": "bold",
            "axes.titlepad": 18,
            "axes.labelsize": 17,
            "axes.labelweight": "bold",
            "axes.linewidth": 1.4,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 12,
            "legend.frameon": False,
            "grid.color": "#D6D6D6",
            "grid.alpha": 0.45,
            "grid.linewidth": 0.8,
        }
    )


def clean_axes(ax: plt.Axes) -> None:
    ax.grid(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=5, width=1.2)


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
    obs = np.atleast_1d(mat["observability_results"])
    return {str(item.method): item for item in obs}

def load_deep_observability_results() -> pd.DataFrame:
    return pd.read_csv(DEEP_CSV)


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_ocp_slope_full_soc(obs_by_method) -> None:
    soc = np.linspace(0.0, 1.0, 700)
    theta_n, theta_p = theta_from_soc(soc)
    un = u_n(theta_n)
    up = u_p(theta_p)
    slope_metric = np.abs(np.gradient(up, soc) - np.gradient(un, soc))

    obs_sample = obs_by_method["FDM"]
    soc_ref = np.asarray(obs_sample.soc_ref).reshape(-1)
    soc_min = float(soc_ref.min())
    soc_max = float(soc_ref.max())

    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    clean_axes(ax)
    ax.plot(soc, slope_metric, color="black", lw=3.0)
    ax.axvspan(soc_min, soc_max, color="#B0B0B0", alpha=0.18, label="HPPC operating window")
    ax.set_xlabel("SOC [-]")
    ax.set_ylabel(r"$|dU_p/dSOC - dU_n/dSOC|$")
    ax.set_title("OCP Slope Metric Across the Full SOC Range")
    ax.legend(loc="best")
    save_figure(fig, FIG_DIR / "obs_fig1_ocp_slope_full_soc.png")


def make_spm_rank_vs_soc(obs_by_method) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    for method in SPM_METHODS:
        obs = obs_by_method[method]
        soc = np.asarray(obs.soc_ref).reshape(-1)
        rank = np.asarray(obs.rank).reshape(-1)
        ax.plot(soc, rank, color=COLORS[method], lw=2.4, label=method)
    clean_axes(ax)
    ax.invert_xaxis()
    ax.set_xlabel("SOC [-]")
    ax.set_ylabel("Observability rank")
    ax.set_title("SPM Family: Rank of O Along the HPPC Trajectory")
    ax.legend(loc="best", ncol=2)
    save_figure(fig, FIG_DIR / "obs_fig2_spm_rank_vs_soc.png")


def make_spm_sigma_nz_vs_soc(deep_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    for method in SPM_METHODS:
        obs = deep_df[deep_df["method"] == method]
        soc = obs["soc_ref"].to_numpy()
        sigma = np.clip(obs["sigma_min_nonzero"].to_numpy(), 1e-30, None)
        ax.semilogy(soc, sigma, color=COLORS[method], lw=2.4, label=method)
    clean_axes(ax)
    ax.invert_xaxis()
    ax.set_xlabel("SOC [-]")
    ax.set_ylabel(r"Smallest nonzero $\sigma(O)$")
    ax.set_title("SPM/FVM: Weakest Observable State Direction Under HPPC Excitation")
    ax.legend(loc="best", ncol=2)
    save_figure(fig, FIG_DIR / "obs_fig3_spm_sigma_nonzero_vs_soc.png")


def make_spm_rank_fraction_vs_states(sweep_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    for method in SPM_METHODS:
        df = sweep_df[sweep_df["method"] == method].sort_values("total_states")
        ax.plot(
            df["total_states"],
            df["min_rank_fraction"],
            color=COLORS[method],
            lw=2.4,
            marker=MARKERS[method],
            ms=6.5,
            label=method,
        )
    clean_axes(ax)
    ax.set_xscale("log")
    ax.set_xlabel("Total states")
    ax.set_ylabel("Minimum rank / n")
    ax.set_title("SPM Family: Normalized Observability Rank vs Model Order")
    ax.legend(loc="best", ncol=2)
    save_figure(fig, FIG_DIR / "obs_fig4_spm_rank_fraction_vs_states.png")


def make_pade_rank_vs_soc(obs_by_method) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    for method in PADE_METHODS:
        obs = obs_by_method[method]
        soc = np.asarray(obs.soc_ref).reshape(-1)
        rank = np.asarray(obs.rank).reshape(-1)
        ax.plot(soc, rank, color=COLORS[method], lw=2.6, label=method)
    clean_axes(ax)
    ax.invert_xaxis()
    ax.set_xlabel("SOC [-]")
    ax.set_ylabel("Observability rank")
    ax.set_title("Padé Family: Rank of O Along the HPPC Trajectory")
    ax.legend(loc="best")
    save_figure(fig, FIG_DIR / "obs_fig5_pade_rank_vs_soc.png")


def make_pade_sigma_nz_vs_soc(deep_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    for method in PADE_METHODS:
        obs = deep_df[deep_df["method"] == method]
        soc = obs["soc_ref"].to_numpy()
        sigma = np.clip(obs["sigma_min_nonzero"].to_numpy(), 1e-30, None)
        ax.semilogy(soc, sigma, color=COLORS[method], lw=2.6, label=method)
    clean_axes(ax)
    ax.invert_xaxis()
    ax.set_xlabel("SOC [-]")
    ax.set_ylabel(r"Smallest nonzero $\sigma(O)$")
    ax.set_title("Padé Models: Smallest Nonzero Singular Value Under HPPC Excitation")
    ax.legend(loc="best")
    save_figure(fig, FIG_DIR / "obs_fig6_pade_sigma_nonzero_vs_soc.png")


def make_effective_rank_vs_soc(deep_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    for method in SUMMARY_METHODS:
        obs = deep_df[deep_df["method"] == method]
        soc = obs["soc_ref"].to_numpy()
        effective_rank = obs["effective_rank"].to_numpy()
        ax.plot(soc, effective_rank, color=COLORS[method], lw=2.5, label=method)
    clean_axes(ax)
    ax.invert_xaxis()
    ax.set_xlabel("SOC [-]")
    ax.set_ylabel("Effective rank")
    ax.set_title("Effective Rank Along the HPPC Trajectory")
    ax.legend(loc="best", ncol=2, fontsize=11)
    save_figure(fig, FIG_DIR / "obs_fig8_effective_rank_vs_soc.png")


def make_pade_effective_rank_vs_soc(deep_df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    for method in PADE_METHODS:
        obs = deep_df[deep_df["method"] == method]
        soc = obs["soc_ref"].to_numpy()
        effective_rank = obs["effective_rank"].to_numpy()
        ax.plot(soc, effective_rank, color=COLORS[method], lw=2.6, label=method)
    clean_axes(ax)
    ax.invert_xaxis()
    ax.set_xlabel("SOC [-]")
    ax.set_ylabel("Effective rank")
    ax.set_title("Local Linear vs Nonlinear Padé: Effective Rank Comparison")
    ax.legend(loc="best")
    save_figure(fig, FIG_DIR / "obs_fig6_pade_effective_rank_vs_soc.png")


def make_representative_rmse_bar(metrics_df: pd.DataFrame) -> None:
    subset = []
    for method, nr in [("FDM", 101), ("FVM-S2", 101), ("Local Linear Padé-ECM", 12), ("Nonlinear SPM-Padé 2", 12), ("Nonlinear SPM-Padé 3", 12)]:
        row = metrics_df[(metrics_df["method"] == method) & (metrics_df["Nr"] == nr)].iloc[0]
        subset.append(row)
    df = pd.DataFrame(subset)

    fig, ax = plt.subplots(figsize=(9.6, 5.8))
    clean_axes(ax)
    colors = [COLORS[m] for m in df["method"]]
    bars = ax.bar(range(len(df)), 1e3 * df["rmse_vs_exp_V"], color=colors, width=0.72)
    ax.set_xticks(range(len(df)))
    ax.set_xticklabels(
        ["FDM\n(101)", "FVM-S2\n(101)", "Local Linear\nPadé-ECM", "Nonlinear\nSPM-Padé 2", "Nonlinear\nSPM-Padé 3"],
        rotation=0,
    )
    ax.set_ylabel("Voltage RMSE vs experiment [mV]")
    ax.set_title("Representative Voltage Error Comparison")
    for bar, val in zip(bars, 1e3 * df["rmse_vs_exp_V"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.15, f"{val:.2f}", ha="center", va="bottom", fontsize=12)
    save_figure(fig, FIG_DIR / "obs_fig7_representative_rmse_bar.png")


def make_accuracy_vs_observability(metrics_df: pd.DataFrame) -> None:
    subset = metrics_df[metrics_df["method"].isin(SUMMARY_METHODS)].copy()
    subset["rmse_mV"] = 1e3 * subset["rmse_vs_exp_V"]

    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    clean_axes(ax)
    for method in SUMMARY_METHODS:
        df = subset[subset["method"] == method].sort_values("total_states")
        ax.scatter(
            df["obs_min_rank_fraction"],
            df["rmse_mV"],
            s=80,
            color=COLORS[method],
            marker=MARKERS[method],
            label=method,
        )
    ax.set_xlabel("Minimum rank / n")
    ax.set_ylabel("Voltage RMSE vs experiment [mV]")
    ax.set_title("Accuracy–Observability Tradeoff Across Model Families")
    ax.legend(loc="best", fontsize=11)
    save_figure(fig, FIG_DIR / "obs_fig8_accuracy_vs_observability.png")


def write_part3_slide_content() -> None:
    text = """# Part 3: Observability

## Slide 1
### Title
Physical Driver of Observability Variation

### Core equations
\\[
\\dot{x} = A x + B u, \\qquad y = h(x,u)
\\]
\\[
C_k = \\frac{\\partial h}{\\partial x}\\bigg|_{x_k,u_k}
\\]

### Recommended figures
- `obs_fig1_ocp_slope_full_soc.png`
- `obs_fig8_accuracy_vs_observability.png`

### Key conclusion
The dominant source of observability variation is the electrochemical output sensitivity, especially the OCP slope, rather than the numerical discretization alone.

### English script
In this first observability slide, the key result is that observability is primarily output-limited. As the OCP slope changes across SOC, the voltage measurement becomes more or less sensitive to internal concentration dynamics. This gives a physical explanation for why the state reconstruction problem becomes stronger in some regions and weaker in others.

## Slide 2
### Title
SPM/FVM Observability Along the HPPC Trajectory

### Core equations
\\[
\\mathcal{O}_k =
\\begin{bmatrix}
C_k \\\\
C_kA \\\\
\\vdots \\\\
C_kA^{n-1}
\\end{bmatrix}
\\]
\\[
\\mathrm{rank}(\\mathcal{O}_k), \\quad \\sigma_{\\min}^{+}(\\mathcal{O}_k)
\\]

### Recommended figures
- `obs_fig2_spm_rank_vs_soc.png`
- `obs_fig3_spm_sigma_nonzero_vs_soc.png`

### Key conclusion
FDM and the FVM family exhibit very similar observability trends because they approximate the same diffusion physics; increasing state dimension does not create additional voltage information.

### English script
Here we focus on the SPM and FVM family. The rank plot shows that these models remain structurally limited along the trajectory, while the smallest nonzero singular value shows that the weakest observable direction remains extremely small. The important interpretation is that adding more spatial states increases model fidelity, but it does not improve how much information terminal voltage contains about those states.

## Slide 3
### Title
Local Linear vs Nonlinear Padé Observability

### Core equations
\\[
\\dot{x}_d = A_d x_d + B_d I
\\]
\\[
c_{s,\\mathrm{surf}}(t) \\approx c_{s,\\mathrm{lin}} + C_d x_d(t)
\\]

### Recommended figures
- `obs_fig5_pade_rank_vs_soc.png`
- `obs_fig6_pade_effective_rank_vs_soc.png`

### Key conclusion
The Local Linear Padé-ECM and the Nonlinear SPM-Padé are not competing versions of the same model, but reduced-order models designed for different purposes; the former is optimized for local voltage fitting, while the latter preserves the unified nonlinear SPM comparison framework.

### English script
This slide separates the two Padé philosophies. The Local Linear Padé-ECM is a control-oriented reduced model built around a chosen operating point, whereas the Nonlinear SPM-Padé keeps the Padé reduction inside the nonlinear SPM framework. That is why they should not be interpreted as simply a low-order and high-order version of the same method. Their voltage errors and their observability properties must be discussed in the context of different modeling goals.

## Slide 4
### Title
Overall Accuracy–Observability Interpretation

### Core equations
\\[
\\kappa(\\mathcal{O}_k) = \\frac{\\sigma_{\\max}(\\mathcal{O}_k)}{\\sigma_{\\min}(\\mathcal{O}_k)}
\\]
\\[
r_{\\mathrm{eff}} = \\exp\\left(-\\sum_i p_i \\log p_i\\right), \\qquad
p_i = \\frac{\\sigma_i}{\\sum_j \\sigma_j}
\\]

### Recommended figures
- `obs_fig7_representative_rmse_bar.png`
- `obs_fig8_effective_rank_vs_soc.png`

### Key conclusion
Voltage accuracy and state reconstructability are related but not equivalent; even when FDM, FVM-S2, and the Local Linear Padé-ECM achieve similar voltage RMSE, their observability remains fundamentally constrained by the electrochemical output map.

### English script
The final message is that a low voltage RMSE does not automatically imply strong state observability. In the whole-trajectory HPPC comparison, FDM, FVM-S2, and the Local Linear Padé-ECM all achieve about 12 mV voltage error, while the Nonlinear SPM-Padé remains closer to 23 mV. However, all models still show weak practical observability in parts of the trajectory. So the most defensible conclusion is that discretization changes the accuracy–computation tradeoff, but the dominant observability bottleneck is still the voltage physics itself.
"""
    (FIG_DIR / "part3_slide_contents.md").write_text(text)


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    setup_style()
    metrics_df = pd.read_csv(METRICS_CSV)
    sweep_df = pd.read_csv(SWEEP_CSV)
    obs_by_method = load_observability_results()
    deep_df = load_deep_observability_results()

    make_ocp_slope_full_soc(obs_by_method)
    make_spm_rank_vs_soc(obs_by_method)
    make_spm_sigma_nz_vs_soc(deep_df)
    make_spm_rank_fraction_vs_states(sweep_df)
    make_pade_rank_vs_soc(obs_by_method)
    make_pade_sigma_nz_vs_soc(deep_df)
    make_representative_rmse_bar(metrics_df)
    make_accuracy_vs_observability(metrics_df)
    make_effective_rank_vs_soc(deep_df)
    make_pade_effective_rank_vs_soc(deep_df)
    write_part3_slide_content()

    print(f"Observability deep-dive figures exported to: {FIG_DIR}")


if __name__ == "__main__":
    main()
