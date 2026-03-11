function [sol, t_eval, param_out] = run_sim_MATLAB(x_initial, param)
% RUN_SIM_MATLAB - Solve SPM using MATLAB's ode15s integrator
%
% INPUTS:
%   x_initial - Initial concentration state vector [mol/m³]
%   param     - Parameter structure
%
% OUTPUTS:
%   sol       - Structure with voltage [V], SOC [-], concentrations [mol/m³], etc.
%   t_eval    - Time vector [s]
%   param_out - Updated parameter structure

    tspan = param.t_data;  % [s]

    % Construct FDM matrices
    switch param.method
        case 'FDM'
            [param.A_mat, param.Bvector] = fdm_matrices(param, param.nstates_electrode);
    end

    fun_ode = @(t, x) param.A_mat * x - param.Bvector .* param.I_interpolant(t);
    initial_slope = fun_ode(tspan(1), x_initial);

    J_base = ones(param.nstates_electrode, 3);
    switch param.method
        case 'FDM'
            Jpattern_i = spdiags(J_base, [-1, 0, 1], param.nstates_electrode, param.nstates_electrode);
    end

    Jpattern = blkdiag(Jpattern_i, Jpattern_i);

    I_interpolated = param.I_interpolant(tspan);
    var_I_interpolated = var(I_interpolated);

    if var_I_interpolated < 1e-9 % constant current discharge
        reltol_target = 1e1;
    else                         % all other types of discharge
        reltol_target = 1e-8;
    end

    Jacobian = odeJacobian('SparsityPattern', Jpattern);

    F = ode('ODEFcn', fun_ode, ...
            'InitialTime', tspan(1), ...
            'InitialValue', x_initial, ...
            'Parameters', [], ...
            'EquationType', 'standard', ...
            'Jacobian', Jacobian, ...
            'EventDefinition', [], ...
            'NonNegativeVariables', [], ...
            'InitialSlope', initial_slope, ...
            'Sensitivity', [], ...
            'Solver', 'ode15s', ...
            'RelativeTolerance', reltol_target, ...
            'AbsoluteTolerance', 1e-11);
    F.SolverOptions.InitialStep = 1e-12;
    F.SolverOptions.MaxStep = 1e0;
    F.SolverOptions.MaxOrder = 5;

    S = solve(F, tspan);
    t_eval = S.Time(:);
    x_out = S.Solution;

    % ====== Separate electrochemical, thermal & aging state variables from x_out matrix ======
    % State variables are now concentrations [mol/m³]
    cs_n = x_out(1:param.nstates_electrode,:);               % Anode solid concentrations [mol/m³]
    cs_p = x_out(param.nstates_electrode+1:end,:);           % Cathode solid concentrations [mol/m³]

    if ~isreal(cs_n)                                                    
        cs_n = abs(cs_n);
    end
    if ~isreal(cs_p)                                                    
        cs_p = abs(cs_p);
    end
    
    % Convert concentrations to stoichiometry for OCP calculations
    sol.theta_s_n = cs_n ./ param.csn_max;                 % Anode stoichiometry [-]
    sol.theta_s_p = cs_p ./ param.csp_max;                 % Cathode stoichiometry [-]

    if ~isreal(sol.theta_s_n)                                                    
        sol.theta_s_n = abs(sol.theta_s_n);
    end
    if ~isreal(sol.theta_s_p)                                                    
        sol.theta_s_p = abs(sol.theta_s_p);
    end
    
    % Store concentrations in output structure
    sol.cs_n = cs_n;
    sol.cs_p = cs_p;
    
    % ==== Kinetics for each cell ====
    num_t = length(t_eval);
    % Open circuit potential and overpotential
    sol.ocp_n = zeros(num_t, 1);
    sol.ocp_p = zeros(num_t, 1);
    sol.eta_n = zeros(num_t, 1);
    sol.eta_p = zeros(num_t, 1);

    % surface concentrations
    sol.theta_s_n_surf = zeros(num_t, 1);
    sol.theta_s_p_surf = zeros(num_t, 1);
    switch param.method
        case 'FDM'

            sol.theta_s_n_surf(:) = sol.theta_s_n(end,:)';
            sol.theta_s_p_surf(:) = sol.theta_s_p(end,:)';

            % Calculate SOC
            [sol.soc_bulk_n, sol.soc_bulk_p, ~, ~] = soc_calculation(sol.theta_s_n, sol.theta_s_p, param, param.method);
    end
    
    % Store surface concentrations
    sol.csn_surf = sol.theta_s_n_surf .* param.csn_max;
    sol.csp_surf = sol.theta_s_p_surf .* param.csp_max;

    % Surface concentration -> surface stoichiometry
    theta_surf_n = sol.theta_s_n_surf(:);
    theta_surf_p = sol.theta_s_p_surf(:);   

    I_eval = param.I_interpolant(t_eval);
    for j = 1:num_t
        sol.ocp_n(j) = U_n(theta_surf_n(j), param.T_amb, param);
        sol.ocp_p(j) = U_p(theta_surf_p(j), param.T_amb, param);
        sol.eta_n(j) = eta_anode(theta_surf_n(j), param.T_amb, I_eval(j), param);  
        sol.eta_p(j) = eta_cathode(theta_surf_p(j), param.T_amb, I_eval(j), param);
    end

    sol.V_oc = sol.ocp_p - sol.ocp_n;

    % total cell potential
    sol.V_cell = sol.V_oc + sol.eta_p - sol.eta_n - param.R0.*I_eval;
        
    param_out = param;

end
