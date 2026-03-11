function [sol, t_eval, param_out] = run_local_linear_pade_ecm_casadi(param, t_data, I_data, T_amb, solver_opts)
%RUN_LOCAL_LINEAR_PADE_ECM_CASADI CasADi simulation for the Local Linear Padé-ECM model.

    import casadi.*

    param.t_data = t_data(:);
    param.I_data = I_data(:);
    param.T_amb = 273.15 + T_amb;
    param.V_thermal = param.Rg .* param.T_amb ./ param.F;

    if ~isfield(solver_opts, 'linearization_soc') || isempty(solver_opts.linearization_soc)
        solver_opts.linearization_soc = param.SOC_IC;
    end

    coeff = local_linear_pade_ecm_coeffs(param, solver_opts.linearization_soc, param.T_amb);

    [~, csn0, csp0] = conc_initial_sd(param, struct('nstates_electrode',1));
    x0 = zeros(8,1);
    x0(7) = csn0;
    x0(8) = csp0;

    Bavg_n = -1 / (param.epsn * param.Acell * param.Ln * param.F);
    Bavg_p =  1 / (param.epsp * param.Acell * param.Lp * param.F);

    param.A_mat = blkdiag(coeff.n.A, coeff.p.A, 0, 0);
    param.Bvector = [coeff.n.B; coeff.p.B; Bavg_n; Bavg_p];
    % Preserve the original local-linear Padé-ECM interpretation:
    % the voltage model is reconstructed as a perturbation around a chosen
    % operating point surface concentration, not from the evolving average state.
    param.Csurf_n = [coeff.n.C, 0, 0, 0, 0, 0];
    param.Csurf_p = [0, 0, 0, coeff.p.C, 0, 0];

    tspan = param.t_data;
    x = MX.sym('x', 8);
    I_sym = MX.sym('I');

    dx = MX.zeros(8,1);
    dx(1:3) = coeff.n.A * x(1:3) + coeff.n.B * I_sym;
    dx(4:6) = coeff.p.A * x(4:6) + coeff.p.B * I_sym;
    dx(7)   = Bavg_n * I_sym;
    dx(8)   = Bavg_p * I_sym;

    ode = struct('x', x, 'p', I_sym, 'ode', dx);

    dt = tspan(2) - tspan(1);
    F = integrator('F_ll_pade_ecm', 'cvodes', ode, 0, dt, struct());
    F_map = F.mapaccum('F_ll_pade_map', numel(tspan)-1);

    res = F_map('x0', x0, 'p', I_data(1:end-1)');
    x_res = [x0, full(res.xf)].';

    t_eval = tspan(:);
    I_eval = I_data(:);

    delta_csn = (coeff.n.C * x_res(:,1:3).').';
    delta_csp = (coeff.p.C * x_res(:,4:6).').';

    csn_avg = x_res(:,7);
    csp_avg = x_res(:,8);
    csn_surf = coeff.csn0 + delta_csn;
    csp_surf = coeff.csp0 + delta_csp;

    theta_avg_n = csn_avg ./ param.csn_max;
    theta_avg_p = csp_avg ./ param.csp_max;
    theta_surf_n = min(max(csn_surf ./ param.csn_max, 1e-8), 1 - 1e-8);
    theta_surf_p = min(max(csp_surf ./ param.csp_max, 1e-8), 1 - 1e-8);

    sol.ocp_n = U_n(theta_surf_n, param.T_amb, param);
    sol.ocp_p = U_p(theta_surf_p, param.T_amb, param);
    sol.eta_n = eta_anode(theta_surf_n, param.T_amb, I_eval, param);
    sol.eta_p = eta_cathode(theta_surf_p, param.T_amb, I_eval, param);

    sol.V_oc = sol.ocp_p - sol.ocp_n;
    sol.V_cell = sol.V_oc + sol.eta_p - sol.eta_n - param.R0 .* I_eval;

    sol.soc_bulk_n = (theta_avg_n - param.theta_n_0) ./ (param.theta_n_100 - param.theta_n_0);
    sol.soc_bulk_p = (theta_avg_p - param.theta_p_0) ./ (param.theta_p_100 - param.theta_p_0);

    sol.theta_s_n = theta_avg_n.';
    sol.theta_s_p = theta_avg_p.';
    sol.cs_n = csn_avg.';
    sol.cs_p = csp_avg.';
    sol.csn_surf = csn_surf;
    sol.csp_surf = csp_surf;

    param_out = param;
end
