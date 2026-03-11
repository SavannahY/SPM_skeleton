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

Validated by rerunning:
- `export_observability_deep_dive_metrics.m`
- `make_observability_deep_dive.py`

Speaker split for the final presentation:
- `Yunhan`: Part 1, FDM vs FVM accuracy/runtime comparison
- `Siyao`: Nonlinear Padé modeling principle and voltage results
- `You`: Part 3, observability calculation, observability results, and cross-family interpretation

## Slide 1
### Title
How Observability Is Computed for FDM and FVM

### Why this slide belongs to Part 3
This slide should not repeat Yunhan's accuracy/runtime comparison. Its purpose is to explain how the observability calculation is performed on top of the FDM/FVM simulations.

### Core equations
\\[
\\dot{x} = A x + B u, \\qquad y = h(x,u)
\\]
\\[
C_k = \\frac{\\partial h}{\\partial x}\\bigg|_{x_k,u_k}, \\qquad
\\mathcal{O}_k =
\\begin{bmatrix}
C_k \\\\
C_kA \\\\
\\vdots \\\\
C_kA^{n-1}
\\end{bmatrix}
\\]

### Figure(s)
- `master_output_dir/option1/observability_deep_dive/obs_fig1_ocp_slope_full_soc.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig2_spm_rank_vs_soc.png`

### Result to say
For FDM and FVM, the state matrix `A` comes directly from the discretized diffusion model, while the time-varying output Jacobian `C_k` is updated at every sample from the nonlinear voltage equation. This means observability is trajectory-dependent and changes with SOC because the OCP slope and overpotential sensitivity change with the operating point.

### English script
In Part 3, I start from the observability calculation itself. For each simulated time step, I linearize the nonlinear voltage equation around the current state and build a local output Jacobian, denoted by C sub k. I then combine that output Jacobian with the state matrix A to form the observability matrix. This lets us evaluate rank and singular-value-based metrics along the full HPPC trajectory instead of treating observability as a single fixed number.

## Slide 2
### Title
Important FDM/FVM Observability Results

### Core equations
\\[
\\mathrm{rank}(\\mathcal{O}_k), \\qquad \\sigma_{\\min}^{+}(\\mathcal{O}_k)
\\]
\\[
r_{\\mathrm{eff}} = \\exp\\left(-\\sum_i p_i \\log p_i\\right), \\qquad
p_i = \\frac{\\sigma_i}{\\sum_j \\sigma_j}
\\]

### Figure(s)
- `master_output_dir/option1/observability_deep_dive/obs_fig2_spm_rank_vs_soc.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig3_spm_sigma_nonzero_vs_soc.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig4_spm_rank_fraction_vs_states.png`

### Result to say
- Representative whole-trajectory result at the deep-dive setting:
  - FDM, 24 states: RMSE = `12.22 mV`, rank range = `15 to 19`
  - FVM-S2, 22 states: RMSE = `12.32 mV`, rank range = `16 to 19`
- Across model order, FDM minimum rank fraction drops from `0.75` at 12 states to `0.168` at 202 states.
- FVM-S2 shows the same pattern: more states improve diffusion resolution, but the voltage measurement does not reveal all added states.

### Interpretation
The key result is that FDM and FVM behave similarly in observability because they discretize the same SPM diffusion physics. Their main limitation is not the solver family, but the fact that terminal voltage carries limited information about high-order diffusion states, especially when the OCP slope becomes small.

### English script
The main FDM and FVM result is that the two families are much closer in observability than their implementation details might suggest. At the representative setting, FDM gives about 12.22 millivolts RMSE with rank varying between 15 and 19, while FVM-S2 gives about 12.32 millivolts RMSE with rank varying between 16 and 19. As the model order grows, the normalized observable fraction drops, which tells us that adding internal states does not create additional measurement information.

## Slide 3
### Title
Extension to Padé: Nonlinear and Local Linear Observability

### Why this slide belongs to Part 3
Siyao explains the nonlinear Padé principle and voltage behavior. Your role here is to extend the same observability framework to both Padé formulations and compare how observability is computed and what it reveals.

### Core equations
For the nonlinear Padé models:
\\[
\\dot{x} = A_{\\mathrm{Pad\\acute{e}}} x + B_{\\mathrm{Pad\\acute{e}}} u, \\qquad
y = h(x,u)
\\]

For the Local Linear Padé-ECM:
\\[
\\dot{x}_d = A_d x_d + B_d I, \\qquad
c_{s,\\mathrm{surf}}(t) \\approx c_{s,0} + C_d x_d(t)
\\]

### Figure(s)
- `master_output_dir/option1/observability_deep_dive/obs_fig5_pade_rank_vs_soc.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig6_pade_effective_rank_vs_soc.png`

### Result to say
- Nonlinear SPM-Padé 2: `6` states, RMSE = `23.31 mV`, rank fixed at `5`
- Nonlinear SPM-Padé 3: `8` states, RMSE = `23.31 mV`, rank range = `5 to 6`
- Local Linear Padé-ECM: `8` states, RMSE = `12.15 mV`, rank fixed at `5`
- Effective-rank range:
  - Nonlinear SPM-Padé 2: `1.86 to 2.91`
  - Nonlinear SPM-Padé 3: `1.80 to 2.88`
  - Local Linear Padé-ECM: `1.01 to 2.00`

### Interpretation
The two Padé approaches should not be described as better versus worse versions of the same model. The nonlinear Padé models stay inside the nonlinear SPM comparison framework, while the Local Linear Padé-ECM compresses the dynamics around one operating point. That is why the Local Linear Padé-ECM can match FDM-level voltage RMSE, but its effective-rank plot still shows that the information is concentrated in a very small observable subspace.

### English script
To extend the observability study beyond FDM and FVM, I applied the same local linearization procedure to both Padé formulations. The nonlinear SPM-Padé models use the reduced diffusion states inside the nonlinear SPM framework, whereas the Local Linear Padé-ECM reconstructs surface concentration perturbations around a chosen operating point. The important result is that the Local Linear Padé-ECM reaches about 12.15 millivolts RMSE, but its effective rank stays only around one to two, so compact voltage accuracy does not mean rich state observability.

## Slide 4
### Title
Cross-Family Comparison: Accuracy Does Not Equal Observability

### Core equations
\\[
\\kappa(\\mathcal{O}_k) = \\frac{\\sigma_{\\max}(\\mathcal{O}_k)}{\\sigma_{\\min}(\\mathcal{O}_k)}
\\]
\\[
\\text{trajectory-dependent observability} \\Longrightarrow
\\text{compare rank, } \\sigma_{\\min}^{+}, r_{\\mathrm{eff}}, \\text{ and RMSE together}
\\]

### Figure(s)
- `master_output_dir/option1/observability_deep_dive/obs_fig7_representative_rmse_bar.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig8_accuracy_vs_observability.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig8_effective_rank_vs_soc.png`

### Result to say
- FDM: `12.22 mV`
- FVM-S2: `12.32 mV`
- Local Linear Padé-ECM: `12.15 mV`
- Nonlinear SPM-Padé 2 and 3: both about `23.31 mV`

### Interpretation
The final takeaway is that voltage fit and observability are related but not identical. FDM, FVM-S2, and Local Linear Padé-ECM all deliver about `12 mV` whole-trajectory RMSE, yet none of them makes the full state set strongly observable over the full HPPC trajectory. The dominant bottleneck remains the output physics, especially the OCP-slope-driven loss of sensitivity in some SOC regions.

### English script
The last comparison brings all model families onto one slide. In terms of voltage RMSE, FDM, FVM-S2, and the Local Linear Padé-ECM are all near 12 millivolts, while the nonlinear Padé models remain near 23 millivolts. But the observability metrics tell a different story: even accurate compact models can remain weakly observable. So the correct conclusion is that model reduction changes the accuracy and complexity tradeoff, while the main observability bottleneck still comes from the electrochemical output map.

## Validated result files
- `master_output_dir/option1/option1_metrics.csv`
- `master_output_dir/option1/option1_observability_sweep_summary.csv`
- `master_output_dir/option1/observability_deep_dive_timeseries.csv`
- `master_output_dir/option1/observability_deep_dive/obs_fig1_ocp_slope_full_soc.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig2_spm_rank_vs_soc.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig3_spm_sigma_nonzero_vs_soc.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig4_spm_rank_fraction_vs_states.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig5_pade_rank_vs_soc.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig6_pade_effective_rank_vs_soc.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig7_representative_rmse_bar.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig8_accuracy_vs_observability.png`
- `master_output_dir/option1/observability_deep_dive/obs_fig8_effective_rank_vs_soc.png`

## Prompt for another slide-making tool
Copy and paste the prompt below into the slide software together with the files above:

```text
Build a 4-slide research-style presentation section titled "Part 3: Observability" for a lithium-ion battery Single Particle Model project.

Important speaker split:
- Yunhan already covers Part 1: FDM vs FVM comparison in terms of method and overall performance.
- Siyao already covers the principle and voltage results of the Nonlinear Padé model.
- This section must focus only on observability calculation, observability results, and cross-family interpretation.

Use the provided markdown file `part3_slide_contents.md` as the authoritative content source.

Use these figures exactly where indicated in the markdown:
- obs_fig1_ocp_slope_full_soc.png
- obs_fig2_spm_rank_vs_soc.png
- obs_fig3_spm_sigma_nonzero_vs_soc.png
- obs_fig4_spm_rank_fraction_vs_states.png
- obs_fig5_pade_rank_vs_soc.png
- obs_fig6_pade_effective_rank_vs_soc.png
- obs_fig7_representative_rmse_bar.png
- obs_fig8_accuracy_vs_observability.png
- obs_fig8_effective_rank_vs_soc.png

Design requirements:
- Academic, clean, research-grade style
- No decorative icons
- Prioritize figures over text
- Keep equations visible and readable
- Put slide titles clearly above the figure body
- Use concise bullets derived from the markdown, not long paragraphs
- Keep speaker notes aligned with the English script in the markdown

Content requirements:
- Slide 1: explain how observability is computed for FDM and FVM
- Slide 2: show the key FDM/FVM observability results
- Slide 3: extend the observability discussion to Nonlinear SPM-Padé and Local Linear Padé-ECM
- Slide 4: compare FDM, FVM, Nonlinear Padé, and Local Linear Padé together and emphasize that voltage accuracy does not imply strong observability

Do not duplicate Yunhan's method-comparison slide or Siyao's Padé-principle slide. This section is specifically the observability part.
```
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
