# Part 3: Observability

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
\[
\dot{x} = A x + B u, \qquad y = h(x,u)
\]
\[
C_k = \frac{\partial h}{\partial x}\bigg|_{x_k,u_k}, \qquad
\mathcal{O}_k =
\begin{bmatrix}
C_k \\
C_kA \\
\vdots \\
C_kA^{n-1}
\end{bmatrix}
\]

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
\[
\mathrm{rank}(\mathcal{O}_k), \qquad \sigma_{\min}^{+}(\mathcal{O}_k)
\]
\[
r_{\mathrm{eff}} = \exp\left(-\sum_i p_i \log p_i\right), \qquad
p_i = \frac{\sigma_i}{\sum_j \sigma_j}
\]

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
\[
\dot{x} = A_{\mathrm{Pad\acute{e}}} x + B_{\mathrm{Pad\acute{e}}} u, \qquad
y = h(x,u)
\]

For the Local Linear Padé-ECM:
\[
\dot{x}_d = A_d x_d + B_d I, \qquad
c_{s,\mathrm{surf}}(t) \approx c_{s,0} + C_d x_d(t)
\]

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
\[
\kappa(\mathcal{O}_k) = \frac{\sigma_{\max}(\mathcal{O}_k)}{\sigma_{\min}(\mathcal{O}_k)}
\]
\[
\text{trajectory-dependent observability} \Longrightarrow
\text{compare rank, } \sigma_{\min}^{+}, r_{\mathrm{eff}}, \text{ and RMSE together}
\]

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
