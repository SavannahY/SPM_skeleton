# ENERGY 295 Winter 2026 Single Particle Model (SPM) for Lithium-Ion Battery Simulation

**Instructor**: Simona Onori
**TA**: Sai Thatipamula
## Overview

This repository contains a dimensional Single Particle Model (SPM) implementation for simulating lithium-ion batteries. The model is calibrated for **LG Chem INR21700** cells (4.85 Ah capacity, NMC cathode, Graphite/Si anode) and includes:

- **Battery simulation** with pre-identified parameters
- **Efficient ODE solving** using CasADi integration with CVODES
- **Validation** against HPPC and drive cycle data
- **Alternative discretization methods** (FVM, CVM) available for comparison

**Contributors**: Anirudh Allam, Trey Weaver, Gabriele Pozzato, Le Xu, Joseph Lucero and Sai Thatipamula
**Date**: January 2026

---

## Table of Contents

1. [Requirements](#requirements)
2. [Project Structure](#project-structure)
3. [Model Description](#model-description)
4. [Quick Start](#quick-start)
5. [Main Scripts](#main-scripts)
6. [Core Functions](#core-functions)
7. [Cell Model Components](#cell-model-components)
8. [Input/Output](#inputoutput)
9. [Parameter Description](#parameter-description)
10. [Citation](#citation)

---

## Student Tasks

**You are required to complete the following file(s):**

| File | Description | Status |
|------|-------------|--------|
| `cell_model/fdm_matrices.m` | Construct the A and B matrices for the Finite Difference Method discretization of spherical diffusion | ❌ **TODO** |

### `fdm_matrices.m`

This function constructs the state-space matrices for the ODE system:

```
dx/dt = A·x + B·I
```

where:
- **A**: State matrix [1/s] containing the discretized diffusion operator
- **B**: Input vector [mol/(m³·A·s)] for the current boundary condition
- **x**: Vector of concentration values at each radial grid point [mol/m³]
- **I**: Applied current [A]

**Function signature:**
```matlab
function [A, B] = fdm_matrices(param, nstates_electrode)
```

**Inputs:**
- `param` - Parameter structure containing diffusivities (`Dsn`, `Dsp`), particle radii (`Rsn`, `Rsp`), grid spacing (`dr_n`, `dr_p`), and other cell parameters
- `nstates_electrode` - Number of spatial nodes per electrode (Nr)

**Outputs:**
- `A` - State matrix of size `[2*Nr, 2*Nr]` (block diagonal for anode and cathode)
- `B` - Input vector of size `[2*Nr, 1]`

**Hints:**
- Refer to `derivations/discretization_methods_derivation.tex` for the mathematical derivation
- Use sparse matrices (`spdiags`) for efficiency
- Remember to apply proper boundary conditions at center (r=0) and surface (r=R)

---

## Requirements

### Software
- **MATLAB** R2020b or later
- **CasADi** 3.7.2+ (included or download from https://web.casadi.org/)
Please download the correct Casadi version from the website above, and
keep the folder in the SPM directory. Then, in run_spm, make sure to
add this folder to your MATLAB path using 
```matlab addpath(genpath('your_casadi_dir_name')) ```.
This is sufficient to run the Single Particle Model. Below, please find 
more information regarding project structure, solver configurations and importing data.
---

## Project Structure

```
SPM_skeleton/
│
├── run_sim.m                     # Main script for battery simulation
├── modelParameters.m             # Parameter definition and initialization
├── plot_ocps.m                   # Plot open circuit potentials
├── README.md                     # This file
│
├── casadi-3.7.2/                 # CasADi library (add to MATLAB path)
│
├── cell_model/                   # Battery electrochemical model
│   ├── SPM_sim.m                 # Main simulation wrapper
│   ├── run_sim_MATLAB.m          # ODE solver (MATLAB integrator)
│   ├── run_sim_CasADi.m          # ODE solver (CasADi integrator)
│   ├── U_n.m                     # Anode open circuit potential
│   ├── U_p.m                     # Cathode open circuit potential
│   ├── eta_anode.m               # Anode overpotential
│   ├── eta_cathode.m             # Cathode overpotential
│   ├── conc_initial_sd.m         # Initial concentration calculator
│   ├── soc_calculation.m         # State of Charge computation
│   ├── volume_average.m          # Volume averaging utility
│   └── fdm_matrices.m            # ⚠️ TODO: Finite Difference Method matrices
│
├── utils/                        # Utility functions
│   ├── load_data.m               # Load and resample experimental data
│   ├── n_states.m                # Calculate number of state variables
│   ├── plot_results.m            # Plotting utilities
│   ├── post_process.m            # Post-process simulation results
│   └── save_results.m            # Save simulation outputs
│
├── data/                         # Input experimental data directory
│   ├── data_HPPC.mat             # HPPC test data
│   ├── data_udds.mat             # UDDS drive cycle data
│   └── data_co20.mat             # CO20 test data
│
├── master_output_dir/            # Simulation results and outputs
│   └── matlab_sim_results_*.mat  # Saved simulation outputs
│
└── derivations/                  # Mathematical derivations
    └── discretization_methods_derivation.tex

```

---

## Model Description

### Single Particle Model (SPM)

The SPM is a reduced-order electrochemical model that represents each electrode (anode and cathode) as a single spherical particle. The model captures:

1. **Solid-phase lithium diffusion** in spherical particles (Fick's second law)
2. **Butler-Volmer kinetics** at electrode-electrolyte interfaces
3. **Open circuit potentials** (OCP) as a function of lithium concentration
4. **Lumped resistance** in the cell

### Governing Equations

**Diffusion in particles:**
```
∂c/∂t = D_eff * (1/r²) * ∂/∂r(r² * ∂c/∂r)
```

**Cell voltage:**
```
V_cell = OCP_p - OCP_n + η_p - η_n - R₀*I
```

Where:
- c: Lithium concentration [mol/m³]
- D_eff: Effective diffusion coefficient [m²/s]
- η: Overpotential (from Butler-Volmer) [V]
- R₀: Lumped resistance [Ω]
- I: Input current [A]
- V_cell: Terminal voltage [V]

### Discretization Methods

The spherical diffusion PDE is converted to a system of ODEs using spatial discretization:

**Primary Method: Finite Difference Method (FDM)**
- Standard grid-based approach with nodes distributed uniformly in radial direction
- Transforms the PDE into: **dx/dt = A·x + B·I**
  - A: State matrix containing diffusion operator [1/s]
  - B: Input vector for current boundary condition [mol/(m³·A·s)]
  - x: Vector of concentration values at each grid point [mol/m³]
- Nr grid points per electrode (typically Nr=101 for high accuracy)
- Direct ODE formulation (no mass matrix required)
- Implemented in `fdm_matrices.m`

All methods convert the diffusion PDE into ODEs that are integrated using CasADi's CVODES solver.

---


## Quick Start

### Battery Simulation

To simulate battery response using pre-identified parameters:

```matlab
% Run simulation
run_sim
```

This will:
- Load model parameters from `modelParameters.m`
- Load experimental data from `./data/` directory
- Simulate battery response to the current profile
- Save detailed outputs (voltage, SOC, concentrations, etc.) to `./master_output_dir/`

**Default settings:**
- Discretization: **FDM (Finite Difference Method)** with Nr=101 grid points
- Integrator: CasADi CVODES
- Initial SOC: 100%
- Temperature: 23°C

---

## Main Scripts

### `run_sim.m`

**Purpose**: Simulate battery with identified parameters

**Configuration:**
```matlab
solver_opts.Nr = [101];         % Number of radial nodes per electrode [-]
solver_opts.method = 'FDM';     % Discretization: 'FDM' (default), 'FVM', or 'CVM'
solver_opts.integrator = 'casadi'; % ODE solver: 'casadi' (default) or 'matlab'
params.SOC_IC = 1.0;            % Initial SOC [-]
params.Q_IC = 4.8768;           % Capacity [A·h]
T_amb = 23;                     % Temperature [°C]
dt = 1.0;                       % Sampling time [s]
```

**Outputs:**
- `matlab_sim_results_*.mat`: Full simulation results

---

## Core Functions

### `modelParameters.m`

Initializes and returns the complete parameter structure with physical constants, geometric parameters, and identified electrochemical parameters.

**Parameter Values:**
```
D_s: 10^-17 to 10^-8 [m²/s]  (solid diffusivity)
ε:   10^-2  to 10^0  [-]       (porosity/volume fraction)
k₀:  10^-10 to 10^0  [m^2.5/(mol^0.5·s)]  (reaction rate)
```

---

## Cell Model Components

### Simulation Functions

#### `SPM_sim.m`
Main simulation wrapper that:
1. Sets up initial conditions
2. Calls the appropriate integrator (CasADi or MATLAB)
3. Returns all simulation outputs including voltage, SOC, and concentration profiles

#### `run_sim_CasADi.m`
Efficient ODE integration using CasADi's symbolic framework with CVODES solver.

**Features:**
- Adaptive BDF time-stepping
- Sparse linear algebra (CSparse)
- GMRES Newton iterations
- Tolerance: reltol=1e-8, abstol=1e-11

### Electrochemistry Functions

#### `U_n.m` and `U_p.m`
Open circuit potential (OCP) functions fitted to experimental data.

**Anode (Graphite):**
```matlab
U_n = 1.9793*exp(-39.3631*θ) + 0.2482 - 0.0909*tanh(...) - ...
```

**Cathode (NMC):**
```matlab
U_p = -0.8090*θ + 4.4875 - 0.0428*tanh(...) - 17.7326*tanh(...) + ...
```

#### `eta_anode.m` and `eta_cathode.m`
Butler-Volmer overpotential calculations:

```matlab
η = 2*V_thermal * asinh(i/(i₀*√(θ*(1-θ))))
```

### Utility Functions

#### `load_data.m`
Loads experimental data from MAT-file, resamples to desired time step, and creates interpolants.

**Inputs:**
- `filename`: Path to data file
- `dt`: Sampling time step [s]
- `Tend`: End time [s]

**Outputs:**
- `t_data`, `I_data`, `V_data`: Resampled time, current, voltage vectors
- `solver_opts`: Structure with interpolants for continuous evaluation

#### `n_states.m`
Calculates the number of state variables per electrode based on discretization method and Nr.

#### `post_process.m`
Extracts and organizes simulation results into a structured format for analysis and visualization.

#### `save_results.m`
Saves simulation results to MAT-file with descriptive filename based on configuration.

#### `plot_results.m`
Provides plotting utilities for visualizing voltage, SOC, concentrations, and overpotentials.

#### `volume_average.m`
Computes volume-averaged concentration in spherical particles:

```matlab
θ_ave = ∫₀¹ 3r²θ(r) dr
```

Keep in mind you will have different implementations based on the spatial discretization scheme chosen.

#### `soc_calculation.m`
Converts surface stoichiometry to State of Charge:

```matlab
SOC_n = (θ_n - θ₀_n)/(θ₁₀₀_n - θ₀_n)
SOC_p = (θ_p - θ₀_p)/(θ₁₀₀_p - θ₀_p)
```

#### `conc_initial_sd.m`
Sets up initial lithium concentration profiles based on desired SOC.

### Discretization Matrices

#### `fdm_matrices.m` (Primary Method) — ⚠️ **TO BE COMPLETED**

> **Student Task**: This file contains only the function signature and comments. You must implement the FDM discretization matrices.

Generates A and B matrices for the FDM discretization:
- **Returns**: `[A, B] = fdm_matrices(param, nstates_electrode)`
- **A matrix**: State matrix [1/s] - contains the diffusion operator
- **B vector**: Input vector [mol/(m³·A·s)] - handles current boundary condition
- **ODE form**: dx/dt = A·x + B·I where x is concentration [mol/m³]
- **Grid**: Uniformly spaced nodes from center (r=0) to surface (r=R)
- **Boundary conditions**: No-flux at center, current-dependent flux at surface

See [Student Tasks](#student-tasks) for detailed implementation guidance.

#### `fvm_matrices.m`
Generates A and B matrices for the FVM discretization:
- **Returns**: `[A, B] = fvm_matrices(param)`
- Control volumes centered at cell interfaces
- Nr-1 states per electrode
- Mass-conservative formulation
- ODE form: dx/dt = A·x + B·I

---

## Input/Output

### Input Data Format

Experimental data files (`.mat`) should contain an `output` structure:
```matlab
output.t  % Time vector [s]
output.I  % Current [A] (positive = discharge, negative in data files)
output.V  % Voltage [V]
```

Note: The `load_data()` utility function handles loading, resampling, and sign correction automatically.

### Output Data Format

Simulation results (`.mat`) contain:
```matlab
all_data.t          % Time [s]
all_data.I          % Current [A]
all_data.V          % Simulated voltage [V]
all_data.V_ref      % Reference voltage [V]
all_data.SOC_ref    % Reference SOC [-]
all_data.Voc        % Open circuit voltage [V]
all_data.OCP_n      % Anode OCP [V]
all_data.OCP_p      % Cathode OCP [V]
all_data.eta_n      % Anode overpotential [V]
all_data.eta_p      % Cathode overpotential [V]
all_data.csn        % Anode concentration profile [mol/m³]
all_data.csp        % Cathode concentration profile [mol/m³]
all_data.csn_ave    % Anode average concentration [mol/m³]
all_data.csp_ave    % Cathode average concentration [mol/m³]
all_data.soc_bulk_n % Anode bulk SOC [-]
all_data.soc_bulk_p % Cathode bulk SOC [-]
all_data.param      % All model parameters
```

---

## Parameter Description

### Physical Constants
- `F = 96487` [C/mol]: Faraday's constant
- `Rg = 8.314` [J/mol/K]: Universal gas constant

### Cell Specifications (LG INR21700)
- Capacity: 4.85 Ah
- Nominal voltage: 3.63 V
- Voltage range: 2.4 - 4.3 V
- Anode: Graphite (c_max = 29583 mol/m³)
- Cathode: NMC (c_max = 51765 mol/m³)

### Geometric Parameters (from Chen et al. 2020)
- Anode thickness (L_n): 85.2 μm [m]
- Cathode thickness (L_p): 75.6 μm [m]
- Cell area (A): 0.1037 m²
- Anode particle radius (R_sn): 3.6 μm [m]
- Cathode particle radius (R_sp): 6.25 μm [m]

### Electrochemical Parameters
These parameters are defined in `modelParameters.m`:
1. **D_sn**: Anode solid diffusion coefficient [m²/s]
2. **ε_n**: Anode active material volume fraction [-]
3. **k₀_n**: Anode reaction rate constant [m^2.5/(mol^0.5·s)]
4. **D_sp**: Cathode solid diffusion coefficient [m²/s]
5. **ε_p**: Cathode active material volume fraction [-]
6. **k₀_p**: Cathode reaction rate constant [m^2.5/(mol^0.5·s)]

### Fixed Parameters
- **R₀ = 0.02983 Ω**: Lumped resistance
- **Stoichiometry ranges**: From experimental OCP curves

---

## Usage Examples

### Example 1: Basic Simulation

```matlab
% Load identified parameters
load('./master_input_dir/output_reiden_FVM_Nr_101_params.mat');

% Set up high-resolution simulation
param_in.Nr = 101;  % Fine grid [-]
param_in.method = 'FVM';
param_in.SOC_IC = 0.8;  % Start at 80% SOC [-]

% Load drive cycle
load('./data/data_udds.mat');
I_data = output.I;  % Current [A]
t_data = output.t;  % Time [s]

% Run simulation
[sol, t_eval, param_out] = ...
    SPM_sim(param_in, t_data, I_data, 25, solver_opts);  % T = 25°C

% Extract results
V_cell = sol.V;      % Terminal voltage [V]
soc_n = sol.soc_bulk_n;  % Anode SOC [-]
soc_p = sol.soc_bulk_p;  % Cathode SOC [-]

% Plot results
figure;
subplot(2,1,1); plot(t_eval, V_cell); ylabel('Voltage [V]');
subplot(2,1,2); plot(t_eval, soc_n, t_eval, soc_p); 
ylabel('SOC [-]'); legend('Anode', 'Cathode');
```

---

## Troubleshooting

### Common Issues

**1. CasADi not found**
```
Error: Undefined function or variable 'casadi'.
```
**Solution:** Add CasADi to path:
```matlab
addpath(genpath('/path/to/casadi-3.7.2'));
```

**2. Complex voltage values**
```
Warning: V_cell contains imaginary components
```
**Solution:** 
- Check parameter bounds (may be physically unrealistic)
- Reduce time step or increase tolerance
- Check initial conditions (SOC too high/low)

---

## Advanced Configuration

### Adjusting ODE Solver Tolerances

In `run_sim_CasADi.m`:
```matlab
opts = struct('reltol', 1e-8,      % Relative tolerance
              'abstol', 1e-11,     % Absolute tolerance
              'max_step_size', 1e0);
```

### Non-uniform Initial Conditions

For concentration gradients at t=0:
```matlab
solver_opts.uniform_initial_cond = false;
```

This creates a parabolic initial profile based on SOC_IC.

---

## Citation

### References

**Cell Parameters:**
- Chen, C.-H., et al., "Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models," *Journal of The Electrochemical Society*, vol. 167, 080534, 2020.

**SPM and Battery Modeling:**
- Allam, A., and Onori, S., "An Interconnected Observer for Concurrent Estimation of Bulk and Surface Concentration in Cathode and Anode of a Lithium-Ion Battery," *IEEE Transactions on Industrial Electronics*, vol. 65, no. 9, pp. 7311–7321, 2018.

- Weaver, T., Allam, A., and Onori, S., "Novel Lithium-Ion Battery Pack Modeling Framework – Series-Connected Case Study," *American Control Conference (ACC)*, Denver, CO, USA, 2020.

- Pozzato, G., Lee, S. B., and Onori, S., "Modeling Degradation for Second-Life Battery: Preliminary Results," *IEEE Conference on Control Technology and Applications (CCTA)*, 2021, pp. 826–831.

- Fasolato, S., Allam, A., Li, X., Lee, D., Ko, J., and Onori, S., "Reduced-Order Model of Lithium-Iron Phosphate Battery Dynamics: A POD–Galerkin Approach," *IEEE Control Systems Letters*, vol. 7, pp. 1117–1122, 2023.

- Xu, L., Cooper, J., Allam, A., and Onori, S., "Comparative Analysis of Numerical Methods for Lithium-Ion Battery Electrochemical Modeling," *Journal of The Electrochemical Society*, vol. 170, no. 12, 120525, 2023.

- Lucero, J.N.E., Xu, L., and Onori, S., "Comparing Mass-Preserving Numerical Methods for the Lithium-Ion Battery Single Particle Model," *Modeling, Estimation and Control Conference (MECC)*, Chicago, IL, USA, 2024.

- Lucero, J.N.E., and Onori, S., "Comparative Nonlinear Observability Analysis of Spatial Discretization Schemes for Lithium-Ion Battery Models," *American Control Conference (ACC)*, Denver, CO, USA, 2025.

---


## License

This code is provided for academic and research purposes. Please check with the author for commercial use.

---

**Last Updated**: January 2026
