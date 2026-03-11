# Final Project Report

## SPM Discretization, Performance, and Observability Analysis

### ENERGY 295 Winter 2026

## 1. Introduction
This project studies how the spatial discretization used in the Single Particle Model (SPM) changes three things that matter for battery modeling and estimation: voltage prediction accuracy, computational cost, and state observability. The motivating question is practical. A battery management system does not directly measure internal lithium concentration states, so even a numerically accurate model is only useful for estimation if those internal states remain observable from current and terminal-voltage data.

The project focuses on the `Option 1` comparison between finite-difference, finite-volume, and Padé reduced-order approximations of the SPM. The available results in `master_output_dir/option1` already contain grid sweeps, voltage traces, and summarized observability metrics. The goal of this report is to interpret those results scientifically, with particular emphasis on how observability is quantified and what the different methods reveal.

## 2. Background
The SPM represents each electrode as a single spherical particle governed by solid-phase diffusion. After discretization in the radial coordinate, the concentration dynamics can be written as

`x_dot = A x + B I`

where `x` is the stacked concentration state vector, `A` is the diffusion operator, and `I` is the applied current. The terminal voltage is nonlinear because it depends on electrode surface stoichiometries through open-circuit potential (OCP) and overpotential terms:

`V = U_p(theta_p) - U_n(theta_n) + eta_p - eta_n - R_0 I`

For state-estimation questions, the important issue is whether the internal concentration states influence voltage strongly enough and distinctly enough to be reconstructed from the measurement history.

## 3. How Observability Is Quantified
The MATLAB workflow in `run_option1_comparison.m` computes observability from a trajectory-dependent linearization of the voltage output. At each sampled time index `k`, the terminal voltage is linearized with respect to the state vector to form

`C_k = partial h / partial x`

where the entries of `C_k` depend on the derivatives of the electrode OCP functions and overpotential functions with respect to stoichiometry. The code builds `C_k` by:

1. Computing the surface stoichiometries `theta_n` and `theta_p`.
2. Numerically differentiating `U_n(theta_n)`, `U_p(theta_p)`, `eta_n(theta_n, I_k)`, and `eta_p(theta_p, I_k)`.
3. Mapping those voltage sensitivities back to the surface-concentration states.

Using the linearized output matrix `C_k` and the discretized state matrix `A`, the classical observability matrix is constructed as

`O_k = [C_k; C_k A; C_k A^2; ...; C_k A^(n-1)]`

with `n` equal to the total number of states. The project data then summarize observability with four numerical indicators:

1. `rank(O_k)`: the number of observable state directions.
2. `rank(O_k) / n`: the observable fraction of the full state vector.
3. `sigma_min(O_k)`: the smallest singular value, which indicates proximity to rank loss.
4. `cond(O_k)`: the condition number, which indicates numerical sensitivity.

The workflow also reports the proxy

`|dU_p/dtheta_p - dU_n/dtheta_n|`

which measures how strongly the difference between cathode and anode OCP slopes contributes to the voltage sensitivity. This scalar quantity helps explain trends, but the full observability result still depends on the complete pair `(A, C_k)`.

## 4. Methods Compared
The study compares six model classes:

- `FDM`: finite difference method.
- `FVM-S0`: finite volume method with zero-order surface reconstruction.
- `FVM-S1`: finite volume method with first-order surface reconstruction.
- `FVM-S2`: finite volume method with second-order surface reconstruction.
- `PADE2`: second-order Padé reduced model.
- `PADE3`: third-order Padé reduced model.

The HPPC dataset is used for the reported comparison. Model order is varied across the sweep, and a representative observability comparison is reported at `Nr = 12`.

## 5. Methodology
The analysis uses the exported CSV files:

- `master_output_dir/option1/option1_metrics.csv`
- `master_output_dir/option1/option1_observability_summary.csv`
- `master_output_dir/option1/option1_observability_sweep_summary.csv`
- `master_output_dir/option1/voltage_vs_time_curves.csv`

To make the workflow reusable, I also wrote a Python report generator at `reporting/generate_option1_report.py`. It scans the output folder, regenerates publication-style figures, and writes a fresh Markdown report each time it is run. This makes the deliverable robust to updated experiments without manually editing plot logic or figure titles.

## 6. Results
### 6.1 Voltage Accuracy
The best overall agreement with experiment is obtained by `FVM-S2` at 22 total states, with

`RMSE vs experiment = 5.3115e-3 V`

This result is slightly better than `FVM-S1` at 22 states (`5.3311e-3 V`) and `FDM` at 24 states (`5.3531e-3 V`). The key observation is that an FVM discretization reaches the same voltage fidelity as FDM with fewer states and lower runtime.

The Padé models perform poorly in this dataset. Both `PADE2` and `PADE3` remain near

`RMSE vs experiment ~ 6.90e-1 V`

across all reported settings, which is roughly two orders of magnitude worse than the grid-based methods. This immediately rules them out for accurate voltage tracking in the present form.

### 6.2 Convergence to the Dense FDM Reference
When compared to the dense `FDM, Nr = 101` reference solution, the FVM models converge steadily as model order increases. The best non-reference result is `FVM-S2` at 200 states with

`RMSE vs dense FDM = 2.4333e-5 V`

This confirms that the FVM family can match the high-resolution FDM solution very closely while preserving a conservative finite-volume structure.

### 6.3 Runtime
The FVM methods are consistently fast at low-to-moderate state count. A particularly strong reduced-order point is `FVM-S2` at 14 states:

- Runtime: `2.9326e-2 s`
- RMSE vs experiment: `5.6600e-3 V`

This is a strong practical operating point because it is both fast and already near the best achievable voltage accuracy in the sweep.

The absolute fastest configurations are the smallest FVM models and some Padé models, but low runtime alone is not sufficient because the Padé models sacrifice too much voltage accuracy.

### 6.4 Lithium Conservation
The finite-volume models show much better conservation behavior than FDM at comparable resolution. At the highest reported order:

- `FDM, 202 states`: max relative inventory drift = `9.8998e-6`
- `FVM-S0, 200 states`: `6.9449e-7`
- `FVM-S1, 200 states`: `6.9449e-7`
- `FVM-S2, 200 states`: `6.9449e-7`

This is consistent with the expected physical advantage of FVM, since the governing diffusion equation is integrated over control volumes in a conservative form.

## 7. Observability Results
### 7.1 Representative Comparison at Nr = 12
The representative observability summary at `Nr = 12` shows a clear ranking:

- `FVM-S0`: minimum rank = `19`, total states = `22`, rank fraction = `0.864`
- `FVM-S1`: minimum rank = `19`, total states = `22`, rank fraction = `0.864`
- `FVM-S2`: minimum rank = `19`, total states = `22`, rank fraction = `0.864`
- `FDM`: minimum rank = `19`, total states = `24`, rank fraction = `0.792`
- `PADE2`: minimum rank = `3`, total states = `6`, rank fraction = `0.500`
- `PADE3`: minimum rank = `4`, total states = `8`, rank fraction = `0.500`

The important interpretation is that the FVM models do not merely match FDM in voltage prediction. At comparable order, they also preserve a larger observable fraction of the internal concentration states. The Padé models are compact, but only about half of their states are observable in this summary, which weakens their value for state estimation.

### 7.2 Observability Across Model Order
The observability sweep shows that **all methods become less observable as the number of states grows**. A few representative trends are:

- `FDM`: rank fraction drops from `0.625` at 8 states to `0.178` at 202 states.
- `FVM-S0`: rank fraction drops from `0.833` at 6 states to `0.185` at 200 states.
- `FVM-S1`: rank fraction drops from `0.833` at 6 states to `0.185` at 200 states.
- `FVM-S2`: rank fraction drops from `0.833` at 6 states to `0.185` at 200 states.
- `PADE2`: remains at `0.500`.
- `PADE3`: remains at `0.500`.

This pattern means the observability rank does increase with model order, but not fast enough to keep pace with the growth of the state dimension. In other words, adding more spatial states does not produce a proportional increase in voltage-identifiable information.

### 7.3 Smallest Singular Value and Conditioning
The smallest singular value for FDM is exactly reported as zero in the sweep summaries, and for the FVM methods it becomes extremely small as the model order increases, eventually reaching values on the order of `1e-40` to `1e-50`. This indicates that the observability matrix is very close to numerical rank deficiency even when the reported rank is nonzero.

An even stronger result is that the finite-condition fraction is zero for every case in the observability sweep. Equivalently, the summarized condition numbers are reported as infinite across the dataset. This means that the observability matrices are numerically ill-conditioned everywhere in the sampled comparisons. Therefore:

- rank fraction is a better comparative metric than condition number here,
- `sigma_min(O)` is more informative than `cond(O)`,
- and any observer built on these models should be designed with numerical regularization and noise sensitivity in mind.

### 7.4 Physical Interpretation
The OCP-slope metric stays near `0.372` for the FDM and FVM families and near `0.809` for the Padé models, yet the Padé models still exhibit poor overall observability and very poor voltage accuracy. This is an important scientific result: a larger scalar voltage-sensitivity proxy does not guarantee better full-state observability. The structure of the state dynamics and the way surface concentration is reconstructed matter just as much.

## 8. Discussion
The strongest overall conclusion is that `FVM-S1` and `FVM-S2` provide the best compromise between accuracy, computational efficiency, mass conservation, and observability. They match or outperform FDM in voltage prediction at lower order, preserve a larger observable fraction at the representative comparison point, and retain the finite-volume conservation advantage.

FDM remains a strong high-fidelity baseline, especially for convergence studies, but it is not the best practical choice when the model is intended for estimation or real-time deployment. As the state count grows, FDM loses observable fraction rapidly while carrying a larger state dimension than the comparable FVM models.

The Padé models are appealing from a runtime standpoint, but for this project they are not competitive as physically meaningful estimation models. Their voltage error is too large and their observable fraction is too small to justify using them over the FVM family.

## 9. Limitations
- The available CSV exports summarize observability over sampled points rather than providing the full time-resolved observability trajectories.
- The current report is based on the HPPC experiment only.
- The observability analysis is local because it uses a linearized output matrix along the trajectory.

These limitations do not invalidate the conclusions, but they do suggest that future work should repeat the same analysis on UDDS and on estimator-in-the-loop simulations.

## 10. Conclusions
This project shows that observability is a necessary companion metric to voltage RMSE and runtime. A method can be accurate yet poorly suited to estimation if the internal states are weakly distinguishable from terminal voltage.

The final conclusions are:

1. Observability degrades as model order increases for every method tested.
2. The FVM family, especially `FVM-S1` and `FVM-S2`, gives the best overall tradeoff between voltage accuracy, computational cost, conservation, and retained observability.
3. The Padé models are fast but not accurate enough and not observable enough to recommend for this project objective.
4. Condition number is not a useful discriminator in this dataset because all sampled cases are numerically ill-conditioned; rank fraction and smallest singular value are more meaningful.

For the next stage of work, the best path is to use a moderate-order FVM model as the basis for EKF or adaptive EKF design.

## 11. Reproducible Figure Generation
To regenerate the report-ready figures from updated CSV results, run:

```bash
uv run python reporting/generate_option1_report.py --input-dir master_output_dir/option1
```

This writes refreshed figures and a generated Markdown report under `master_output_dir/option1/python_report`.
