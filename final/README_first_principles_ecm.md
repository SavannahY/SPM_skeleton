# First-Principles Equivalent Circuit Model (Padé-reduced SPM)

This folder extends the ENERGY 295 midterm SPM skeleton with two additions:

1. `cell_model/fdm_matrices.m`
   - completes the missing midterm FDM discretization.
2. `cell_model/demo_first_principles_ecm.m`
   - implements a compact Padé-reduced single-particle model inspired by
     **Prasad & Rahn (2012), "Development of a First Principles Equivalent Circuit Model for a Lithium Ion Battery"**.

## What is included

- A practical 8-state MATLAB model:
  - 3 Padé states for anode surface concentration
  - 3 Padé states for cathode surface concentration
  - 2 average-concentration states for SOC bookkeeping
- The nonlinear voltage output from the midterm skeleton is preserved:
  - `U_p(theta_p)` and `U_n(theta_n)`
  - `eta_cathode(...)` and `eta_anode(...)`
  - `V = U_p - U_n + eta_p - eta_n - R0 I`
- A script `demo_first_principles_ecm.m` that runs on the supplied midterm datasets.

## Important modeling note

The original Rahn/Prasad paper is a **local linear model around one SOC operating point**. To make the model usable with the ENERGY 295 midterm dataset/skeleton, this implementation keeps the paper's **Padé diffusion reduction** but uses the skeleton's **nonlinear OCP and Butler–Volmer output equations**.

That means:
- the state reduction follows the paper,
- but the voltage reconstruction is more practical than the paper's pure small-signal impedance form.

This is deliberate, and it makes the code easier to use on HPPC / C/20 / UDDS data.

## Files added or changed

- `cell_model/fdm_matrices.m` — completed FDM midterm TODO
- `cell_model/dUdtheta_fd.m` — numerical OCP derivative helper
- `cell_model/pade_spm_coeffs.m` — Padé/linearization coefficients
- `cell_model/demo_first_principles_ecm.m` — core Padé-reduced simulator
- `demo_first_principles_ecm.m` — demo script using the midterm dataset
- `utils/n_states.m` — function name fixed to match the filename

## Quick start

In MATLAB, from this folder:

```matlab
demo_first_principles_ecm
```

By default the script uses:
- `./data/data_HPPC.mat`
- `SOC_IC = 1.0`
- `linearization_soc = 1.0`

For a different profile, edit the top of `demo_first_principles_ecm.m`:

```matlab
data_file = './data/data_co20.mat';
% or
data_file = './data/data_udds.mat';
```

## Paper-style ECM coefficients

`pade_spm_coeffs.m` also computes the paper-style linearized coefficients (`K`, `alpha1`, `alpha2`, `beta1`, `beta2`) and the Table 2 RC values.

Depending on the OCP slope at the chosen linearization SOC, some of those RC values can become non-passive (negative). That is a consequence of the local linearization and sign conventions. For simulation, use the Padé state-space model in `demo_first_principles_ecm.m`.
