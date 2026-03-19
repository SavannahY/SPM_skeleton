# ENERGY 295 Final Project Submission

## Project Title
Accuracy, Padé-Based Reduced-Order Approximation, and Trajectory-Dependent Observability of Single Particle Battery Models Under HPPC Excitation

## Group Members
- Jane Yang
- Siyao Wu
- Yunhan Huang

Instructor: Simona Onori  
TA: Sai Thatipamula

## Package Purpose
This submission package is organized as a compact, mirrored version of the project repository so that the main MATLAB and Python workflows can be understood and, in most cases, rerun using the same relative paths as the original project.

It includes:
- final MATLAB and Python code used for the submitted analysis
- the final report source file
- the figure assets referenced by the report
- the saved CSV results used to regenerate the Part I and Part III figures
- the raw data files used by the project workflows

## Folder Structure

```text
final_submission_package/
├── README.md
├── modelParameters.m
├── run_sim.m
├── export_observability_deep_dive_metrics.m
├── make_observability_deep_dive.py
├── pyproject.toml
├── uv.lock
├── cell_model/
├── utils/
├── data/
├── reporting/
│   └── generate_option1_part1_figures.py
└── master_output_dir/
    └── option1/
        ├── option1_metrics.csv
        ├── option1_observability_summary.csv
        ├── option1_observability_sweep_summary.csv
        ├── voltage_vs_time_curves.csv
        ├── observability_deep_dive_timeseries.csv
        ├── observability_deep_dive/
        └── python_report/
            ├── final_project_report_overleaf.tex
            └── figures_part1/
```

## Deliverables Covered
This package is intended to satisfy the course deliverables for:

### Code Submission
- MATLAB and Python code developed for the project
- a README explaining package contents and how to run the code
- scripts that reproduce the final reported Part I and Part III figures

### Written Report
- final report source (`master_output_dir/option1/python_report/final_project_report_overleaf.tex`)
- figure files referenced by the report

## Software Requirements
To run the included workflows, the following are needed:

- MATLAB
- CasADi installed for MATLAB and added to the MATLAB path
- Python 3.10 or later
- `uv` (recommended), or a Python environment with:
  - `matplotlib`
  - `numpy`
  - `pandas`
  - `scipy`

## How to Run the Code and Reproduce the Final Results

### Important setup note for MATLAB / CasADi
The MATLAB scripts in this package assume that CasADi is installed and
added to the MATLAB path. Before running the scripts, open MATLAB in the
root of `final_submission_package/` and check the CasADi path lines in:

- `run_sim.m`
- `export_observability_deep_dive_metrics.m`

If needed, update the CasADi path so it matches your local installation.

### MATLAB-first workflow

#### 1. Single simulation run
To run a direct SPM simulation from MATLAB, use:

```matlab
run_sim
```

This script:
- loads one experimental dataset (currently HPPC in the packaged copy)
- sets the discretization method and state count
- runs the battery simulation through `SPM_sim`
- computes basic validation metrics such as RMSE and conservation error

This is the main MATLAB entry point for a single-case simulation.

#### 2. Part III observability deep-dive results
To regenerate the representative observability results used in Part III,
run:

```matlab
export_observability_deep_dive_metrics
```

This script:
- runs the representative HPPC comparison for
  `FDM`, `FVM-S0`, `FVM-S1`, `FVM-S2`,
  `Nonlinear SPM-Padé 2`, `Nonlinear SPM-Padé 3`,
  and `Local Linear Padé-ECM`
- computes local observability metrics along the trajectory
- writes the deep-dive CSV output to:

```text
master_output_dir/option1/observability_deep_dive_timeseries.csv
```

### Where the Part II Pad\'e code is located
The main Part II Pad\'e-related code is located in the `cell_model/`
folder:

- `cell_model/pade_matrices.m`
  - Pad\'e reduced diffusion matrices for the nonlinear Pad\'e-SPM models
- `cell_model/local_linear_pade_ecm_coeffs.m`
  - coefficient generation for the local linear Pad\'e-ECM formulation
- `cell_model/run_local_linear_pade_ecm_casadi.m`
  - CasADi-based simulation routine for the local linear Pad\'e-ECM
- `cell_model/SPM_sim.m`
  - shared simulation wrapper used by the SPM-family models, including the
    nonlinear Pad\'e variants when selected through solver options

Supporting electrochemical output functions used by the Pad\'e workflows
are also in `cell_model/`, including:

- `U_n.m`
- `U_p.m`
- `eta_anode.m`
- `eta_cathode.m`

### Python post-processing scripts
The Python scripts are used for report-quality figure generation after the
MATLAB or saved-result outputs are available.

#### 3. Part I report figures
The packaged Part I figure generator uses the saved Option 1 CSV summaries:

```bash
uv run python reporting/generate_option1_part1_figures.py
```

This writes:

```text
master_output_dir/option1/python_report/figures_part1/
```

and produces:
- `grid_convergence.png`
- `runtime_comparison.png`
- `li_conservation.png`

#### 4. Part III report figures
After running `export_observability_deep_dive_metrics` in MATLAB, generate
the final observability figures with:

```bash
uv run python make_observability_deep_dive.py
```

This writes:

```text
master_output_dir/option1/observability_deep_dive/
```

The main Part III figures used in the report include:
- `obs_fig1_ocp_slope_full_soc.png`
- `obs_fig2_spm_rank_vs_soc.png`
- `obs_fig3_spm_sigma_nonzero_vs_soc.png`
- `obs_fig4_spm_rank_fraction_vs_states.png`
- `obs_fig6_pade_effective_rank_vs_soc.png`
- `obs_fig8_effective_rank_vs_soc.png`

### Note on Part I full-sweep results
The packaged Part I report figures are regenerated from the saved
`option1_metrics.csv` summaries included in this package. The original
broader sweep workflow that produced those summaries is not packaged here
as a single standalone MATLAB driver script, so the provided Python script
is the intended way to regenerate the final Part I report figures from the
saved results.

## Notes on Final Code State
- The final FDM implementation in `cell_model/fdm_matrices.m` incorporates the TA-provided update.
- The observability deep-dive figures were refreshed after integrating the updated FDM matrices.
- The final observability section of the report was revised to make the methodology explicit, including:
  - the nonlinear output map used before linearization
  - the symbolic form of the local observability matrix
  - the numerical-rank tolerance
  - the entropy-based effective-rank definition

## Report Compilation
The final report source is located at:

```text
master_output_dir/option1/python_report/final_project_report_overleaf.tex
```

To compile the report successfully, keep the report source together with the figure files in:
- `master_output_dir/option1/python_report/figures_part1/`
- `master_output_dir/option1/observability_deep_dive/`

## Attribution
This is a group project completed by Jane Yang, Siyao Wu, and Yunhan Huang for ENERGY 295 under the instruction of Professor Simona Onori, with guidance from TA Sai Thatipamula.
