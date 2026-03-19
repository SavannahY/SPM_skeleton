%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Single Particle Model Simulator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Joseph N. E. Lucero %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sai Thatipamula     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Date: 2025/11/25 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DESCRIPTION:
%   Main script for simulating battery response using previously identified parameters.
%   This script loads optimized parameters and simulates the battery under a specified
%   current profile (e.g., HPPC test, drive cycle, constant current).
%
% PURPOSE:
%   - Simulate model with identified parameters against experimental data
%   - Analyze internal battery states (concentrations, overpotentials, etc.)
%
% INPUTS:
%   - Identified parameters from modelParameters.m
%   - Experimental current/voltage profile from ./data/
%
% OUTPUTS:
%   - Complete simulation results saved to ./master_output_dir/
%   - Includes: voltage, SOC, concentrations, OCPs, overpotentials, etc.
%   - Results can be visualized using utils/plot_results.m
%
% NOTES:
%   - Can use high resolution (Nr=101) since we're not optimizing or estimating parameters
%   - Typical runtime: 5-30 seconds depending on Nr and profile length
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean workspace
restoredefaultpath;
clear all;
clc;
close all;




%% PATH SETUP
% 1. Restore the default MATLAB path to avoid path conflicts
restoredefaultpath; 
clear all; clc; close all;

% 2. Add the actual CasADi absolute path used on this machine

addpath(genpath('/Users/zhengjieyang/Documents/MATLAB/SPM_skeleton/casadi-3.7.2-osx_arm64-matlab2018b'));

% casadi_dir = 'C:\Users\86189\Desktop\Stanford courses\295\SPM_skeleton\SPM_skeleton\casadi-3.7.2-windows64-matlab2018b';
% addpath(casadi_dir);

% 3. Add the project code folders
addpath('./cell_model');
addpath('./utils');

%% LOAD EXPERIMENTAL DATA
% Define output directory for results
output_dir = './master_output_dir';

% Note: Experimental data loading handled by load_data() function below
% Data file should contain:
%   t - time vector [s]
%   I - current vector [A] (positive = discharge)
%   V - voltage vector [V]
% load('./data/data_hppc.mat')

%% SETUP TIME VECTOR AND LOAD DATA
% Sampling time for simulation [seconds]
% Data will be resampled to this time step
dt = 1.0;  % [s]

% End time for simulation [seconds]
% Set to inf to simulate entire profile
Tend = inf;  % [s]



% Load and resample experimental data
[t_data, I_data, V_data, solver_opts] = load_data('./data/data_HPPC.mat', dt, Tend);

% Specify test profile identifier for output file naming
solver_opts.dchg_type_str = 'HPPC';  % 'HPPC', 'UDDS', '1C', etc.

% Change the filename to data_Co20.mat to run the constant-current case
%[t_data, I_data, V_data, solver_opts] = load_data('./data/data_Co20.mat', dt, Tend);
%solver_opts.dchg_type_str = 'C20'; % Update the identifier as well for cleaner saved outputs






%% SET SIMULATION CONFIGURATION
% Spatial discretization method
% 'FDM' = Finite Difference Method

% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% solver_opts.method = 'FDM';

solver_opts.method = 'FVM';

% Number of discretization nodes per electrode
% Higher values = more accurate but slower (Nr=4 for identification, Nr=101 for simulation)
solver_opts.Nr = [101];  % [-] dimensionless grid points

% ODE integrator
% 'casadi' = Fast symbolic integration (recommended)
% 'matlab'  = Standard MATLAB ode15s
solver_opts.integrator = 'casadi';

% Initial condition type
% true  = Uniform concentration throughout particle
% false = Non-uniform (parabolic) profile
solver_opts.uniform_initial_cond = true;

%% LOAD IDENTIFIED PARAMETERS
% Initialize parameter structure with physical and geometric constants
% Includes optimized values for: D_s (diffusivity), epsilon (porosity), k0 (reaction rate), etc.
params = modelParameters(solver_opts.Nr);

% Nominal cell capacity [Ah]
% Used for SOC calculation and normalization
params.Q_IC = 4.8768;  % [A·h]

%% SET INITIAL CONDITIONS
% Initial State of Charge [dimensionless, 0.0 to 1.0]
% 1.0 = fully charged, 0.0 = fully discharged
params.SOC_IC = 1;  % [-]

% Ambient temperature [°C]
% Model assumes isothermal operation at this constant temperature
T_amb = 23;  % [°C]

%% CONFIGURE STATE DIMENSIONS

% Calculate number of state variables per electrode based on discretization method and Nr
solver_opts.nstates_electrode = n_states(solver_opts);  % [-] number of states per electrode

% Total state variables = 2 electrodes (anode + cathode) × states per electrode
solver_opts.nstates_tot = 2.*solver_opts.nstates_electrode;  % [-] total states
%% RUN BATTERY SIMULATION
% Start timing
tic;

% Call main Single Particle Model simulation function
% SPM_sim integrates ODEs and computes:
%   - Cell terminal voltage vs. time
%   - SOC evolution in both electrodes
%   - Solid-phase lithium concentration profiles in particles
%   - Surface concentrations (for Butler-Volmer kinetics)
%   - Open circuit potentials (OCPs) from thermodynamic curves
%   - Activation overpotentials at each electrode

[sol, ...                           % Solution structure with state trajectories
t_eval, ...                         % Time points at which solution was computed
param_out] = ...                    % Parameter structure with expanded fields
    SPM_sim(params, t_data, I_data, T_amb, solver_opts);

% Record computation time
calc_time = toc;

%% POST-PROCESS RESULTS
all_data = post_process(sol, param_out, t_eval, solver_opts);

%% SAVE RESULTS TO FILE
save_results(output_dir, all_data, solver_opts)

%% Plot results (optional)
plot_results(all_data);



%% ================== PART 2: PHYSICAL CONSTRAINTS & VALIDATION ==================
fprintf('\n--- Physical Constraints & Validation Report ---\n');

% 1. Concentration bounds check (using the sol structure)
% Check whether negative-electrode concentration stays within [0, csn_max]
csn_min = min(sol.cs_n(:)); 
csn_max_sim = max(sol.cs_n(:));
csn_limit = params.csn_max; % expected value: 29583 [mol/m^3]

% Check whether positive-electrode concentration stays within [0, csp_max]
csp_min = min(sol.cs_p(:));
csp_max_sim = max(sol.cs_p(:));
csp_limit = params.csp_max; % expected value: 51765 [mol/m^3]

fprintf('Negative electrode concentration range: [%.2f, %.2f] mol/m^3 (Limit: %d)\n', ...
    csn_min, csn_max_sim, csn_limit);
fprintf('Positive electrode concentration range: [%.2f, %.2f] mol/m^3 (Limit: %d)\n', ...
    csp_min, csp_max_sim, csp_limit);

% Validation logic
if csn_min < -1e-6 || csn_max_sim > csn_limit + 1e-6 || ...
   csp_min < -1e-6 || csp_max_sim > csp_limit + 1e-6
    fprintf('Result: [FAILED] Concentration violated physical bounds!\n');
else
    fprintf('Result: [PASSED] Concentration within physical bounds.\n');
end

% 2. Lithium conservation check
% Vol_n = params.Acell * params.Ln * params.epsn; 
% Vol_p = params.Acell * params.Lp * params.epsp;
% 
% % 【关键修正】使用 (:) 确保向量方向一致，防止生成巨型矩阵
% li_n = Vol_n * (sol.soc_bulk_n(:) * params.csn_max);
% li_p = Vol_p * (sol.soc_bulk_p(:) * params.csp_max);
% 
% n_Li_tot = li_n + li_p;
% 
% % 计算相对偏差
% conservation_error = (max(n_Li_tot) - min(n_Li_tot)) / n_Li_tot(1) * 100;

Vol_n = params.Acell * params.Ln * params.epsn;
Vol_p = params.Acell * params.Lp * params.epsp;

% compute average concentration in each electrode
cavg_n = mean(sol.cs_n,1)';
cavg_p = mean(sol.cs_p,1)';

% lithium inventory
li_n = Vol_n * cavg_n;
li_p = Vol_p * cavg_p;

n_Li_tot = li_n + li_p;

% conservation error
conservation_error = (max(n_Li_tot) - min(n_Li_tot)) / n_Li_tot(1) * 100;



fprintf('Total Lithium Conservation Error: %.4e %%\n', conservation_error);

% 3. 电压范围检查 [cite: 125]
V_sim = sol.V_cell;
fprintf('Voltage Range: [%.4f, %.4f] V\n', min(V_sim), max(V_sim));

% 计算 RMSE
rmse_val = sqrt(mean((V_sim - V_data).^2)); 
fprintf('\nSimulation RMSE: %.4f V (%.2f mV)\n', rmse_val, rmse_val*1000);




