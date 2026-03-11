function [cs_initial, cs_n0, cs_p0, theta_s_n0, theta_s_p0] = conc_initial_sd(param, solver_opts)
    % Compute initial concentration values based on initial SOC
    % OUTPUTS:
    %   cs_initial   - Initial concentration vector [mol/m³] for all states
    %   cs_n0        - Anode initial concentration [mol/m³]
    %   cs_p0        - Cathode initial concentration [mol/m³]
    %   theta_s_n0   - Anode initial stoichiometry [-]
    %   theta_s_p0   - Cathode initial stoichiometry [-]

    % Calculate stoichiometry from SOC
    theta_s_n0 = (param.SOC_IC*(param.theta_n_100 - param.theta_n_0) + param.theta_n_0);  % [-]
    theta_s_p0 = (param.theta_p_0 - param.SOC_IC*(param.theta_p_0 - param.theta_p_100));  % [-] 
    
    % Convert stoichiometry to concentration [mol/m³]
    cs_n0 = theta_s_n0 * param.csn_max;  % [mol/m³]
    cs_p0 = theta_s_p0 * param.csp_max;  % [mol/m³]
    
    cs_n_initial = cs_n0 .* ones(solver_opts.nstates_electrode, 1);
    cs_p_initial = cs_p0 .* ones(solver_opts.nstates_electrode, 1);
    
    cs_initial = [cs_n_initial; cs_p_initial];

end
