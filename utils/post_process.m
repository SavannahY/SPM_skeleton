function all_data = post_process(sol, param_out, t_data, solver_opts)
% POST_PROCESS - Package simulation results into unified structure
%
% INPUTS:
%   sol         - Solution structure from SPM_sim
%   param_out   - Parameter structure
%   t_data      - Time vector [s]
%   solver_opts - Solver options structure
%
% OUTPUTS:
%   all_data - Structure with all results including units in comments

    % Calculate volume-averaged concentration in each electrode
    % This gives a single representative concentration for the entire particle
    [all_data.csn_ave, all_data.csp_ave] = volume_average(sol.cs_n, sol.cs_p, param_out, solver_opts.method);
    %% PACKAGE ALL RESULTS FOR SAVING
    % Time and input data
    all_data.t          = t_data;                               % Time [s]
    all_data.I          = solver_opts.I_interpolant(all_data.t);  % Current [A]
    all_data.V_ref      = solver_opts.V_interpolant(all_data.t);  % Reference voltage [V]
    all_data.SOC_ref    = max(0.0, param_out.SOC_IC - cumtrapz(t_data, all_data.I)./(param_out.Q_IC.*3600));  % Reference SOC [-]

    % Concentration data
    all_data.csp        = sol.cs_p;        % Cathode concentration profile [mol/m³]
    all_data.csn        = sol.cs_n;        % Anode concentration profile [mol/m³]
    all_data.cssurp     = sol.csp_surf;    % Cathode surface concentration [mol/m³]
    all_data.cssurn     = sol.csn_surf;    % Anode surface concentration [mol/m³]

    % Voltage and potential data
    all_data.V          = sol.V_cell;      % Simulated terminal voltage [V]
    all_data.Voc        = sol.V_oc;        % Open circuit voltage [V]
    all_data.OCP_p      = sol.ocp_p;       % Cathode OCP [V]
    all_data.OCP_n      = sol.ocp_n;       % Anode OCP [V]
    all_data.eta_p      = sol.eta_p;       % Cathode overpotential [V]
    all_data.eta_n      = sol.eta_n;       % Anode overpotential [V]

    % State of charge data
    all_data.soc_bulk_n = sol.soc_bulk_n;  % Anode bulk SOC [-]
    all_data.soc_bulk_p = sol.soc_bulk_p;  % Cathode bulk SOC [-]
    % Model parameters
    all_data.param      = param_out;   % Complete parameter structure
end