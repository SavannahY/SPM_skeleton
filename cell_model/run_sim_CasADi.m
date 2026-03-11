function [output, t_eval, param_out] = run_sim_CasADi(x_initial, param, solver_opts)
% RUN_SIM_CASADI - Solve SPM using CasADi integrator
%
% Supports:
%   FDM
%   FVM (S0/S1/S2)
%   PADE2
%   PADE3

import casadi.*

tspan = param.t_data;

%% Unpack options
method = upper(solver_opts.method);
nstates_electrode = solver_opts.nstates_electrode;

if ~isfield(solver_opts,'fvm_scheme')
    solver_opts.fvm_scheme = 'S0';
end

%% Build system matrices

switch method

    case 'FDM'

        [param.A_mat, param.Bvector] = ...
            fdm_matrices(param, nstates_electrode);

    case 'FVM'

        [param.A_mat, param.Bvector] = ...
            fvm_matrices(param, nstates_electrode);

    case {'PADE2','PADE3'}

        pade = pade_matrices(param, method);

        param.A_mat   = pade.A_mat;
        param.Bvector = pade.Bvector;

        param.Csurf_n = pade.Csurf_n;
        param.Dsurf_n = pade.Dsurf_n;

        param.Csurf_p = pade.Csurf_p;
        param.Dsurf_p = pade.Dsurf_p;

        param.pade_order = pade.order;

    otherwise

        error('Unknown discretization method')

end

%% Solve ODE system

sundials_t = MX.sym('t');
sundials_x = MX.sym('sundials_x', size(x_initial));

I_interpolated = solver_opts.I_interpolant(tspan);
interp_current = interpolant('interp_current','linear',{tspan},I_interpolated);

I_cell = interp_current(sundials_t);

odefun = param.A_mat * sundials_x - param.Bvector .* I_cell;

odesys = struct('x',sundials_x,'ode',odefun,'alg',[],'quad',[],'t',sundials_t);

opts = struct( ...
    'reltol',1e-8,...
    'abstol',1e-11,...
    'linear_solver','csparse',...
    'newton_scheme','gmres',...
    'linear_multistep_method','bdf',...
    'max_multistep_order',5,...
    'max_step_size',1e0,...
    'always_recalculate_jacobian',false,...
    'show_eval_warnings',false);

F = integrator('F','cvodes',odesys,tspan(1),tspan,opts);

x0 = x_initial(:);

sol = F('x0',x0);

x_out = full(sol.xf);

t_eval = tspan;

%% Separate states

switch method

    case {'FDM','FVM'}

        cs_n = x_out(1:nstates_electrode,:);
        cs_p = x_out(nstates_electrode+1:end,:);

        theta_s_n = cs_n ./ param.csn_max;
        theta_s_p = cs_p ./ param.csp_max;

    case {'PADE2','PADE3'}

        ord = param.pade_order;

        ne = 1 + ord;

        xn = x_out(1:ne,:);
        xp = x_out(ne+1:end,:);

        csn_avg = xn(1,:);
        csp_avg = xp(1,:);

        xdev_n = xn(2:end,:);
        xdev_p = xp(2:end,:);

        cs_n = csn_avg;
        cs_p = csp_avg;

        theta_s_n = csn_avg ./ param.csn_max;
        theta_s_p = csp_avg ./ param.csp_max;

    otherwise

        error('Unknown method')

end

%% Initialize arrays

num_t = length(t_eval);

ocp_n = zeros(num_t,1);
ocp_p = zeros(num_t,1);

eta_n = zeros(num_t,1);
eta_p = zeros(num_t,1);

theta_s_n_surf = zeros(num_t,1);
theta_s_p_surf = zeros(num_t,1);

I_eval = solver_opts.I_interpolant(t_eval);

%% Surface concentration extraction

switch method

    case 'FDM'

        theta_s_n_surf(:) = theta_s_n(end,:)';
        theta_s_p_surf(:) = theta_s_p(end,:)';

        [soc_bulk_n,soc_bulk_p,~,~] = ...
            soc_calculation(theta_s_n,theta_s_p,param,method);

    case 'FVM'

        switch upper(solver_opts.fvm_scheme)

            case 'S0'

                cs_n_surf = cs_n(end,:);
                cs_p_surf = cs_p(end,:);

            case 'S1'

                cs_n_surf = 1.5*cs_n(end,:) - 0.5*cs_n(end-1,:);
                cs_p_surf = 1.5*cs_p(end,:) - 0.5*cs_p(end-1,:);

            case 'S2'

                cs_n_surf = cs_n(end,:) ...
                    + 0.5*(cs_n(end,:) - cs_n(end-1,:)) ...
                    + (1/6)*(cs_n(end,:) - 2*cs_n(end-1,:) + cs_n(end-2,:));

                cs_p_surf = cs_p(end,:) ...
                    + 0.5*(cs_p(end,:) - cs_p(end-1,:)) ...
                    + (1/6)*(cs_p(end,:) - 2*cs_p(end-1,:) + cs_p(end-2,:));

            otherwise

                error('Unknown FVM surface scheme')

        end

        theta_s_n_surf = (cs_n_surf ./ param.csn_max)';
        theta_s_p_surf = (cs_p_surf ./ param.csp_max)';

        [soc_bulk_n,soc_bulk_p,~,~] = ...
            soc_calculation(theta_s_n,theta_s_p,param,method);

    case {'PADE2','PADE3'}

        cs_n_surf = csn_avg + param.Csurf_n(2:end)*xdev_n + param.Dsurf_n .* I_eval';
        cs_p_surf = csp_avg + param.Csurf_p(2:end)*xdev_p + param.Dsurf_p .* I_eval';

        theta_s_n_surf = (cs_n_surf ./ param.csn_max)';
        theta_s_p_surf = (cs_p_surf ./ param.csp_max)';

        soc_bulk_n = theta_s_n(:);
        soc_bulk_p = theta_s_p(:);

end

%% Clamp stoichiometry

eps_theta = 1e-6;

theta_surf_n = min(max(theta_s_n_surf,eps_theta),1-eps_theta);
theta_surf_p = min(max(theta_s_p_surf,eps_theta),1-eps_theta);

%% Voltage calculation

for j = 1:num_t

    ocp_n(j) = U_n(theta_surf_n(j),param.T_amb,param);
    ocp_p(j) = U_p(theta_surf_p(j),param.T_amb,param);

    eta_n(j) = eta_anode(theta_surf_n(j),param.T_amb,I_eval(j),param);
    eta_p(j) = eta_cathode(theta_surf_p(j),param.T_amb,I_eval(j),param);

end

V_oc = ocp_p - ocp_n;

V_cell = V_oc + eta_p - eta_n - param.R0 .* I_eval;

%% Save outputs

output.V_cell = V_cell;

output.ocp_n = ocp_n;
output.ocp_p = ocp_p;

output.V_oc = V_oc;

output.eta_n = eta_n;
output.eta_p = eta_p;

output.soc_bulk_n = soc_bulk_n;
output.soc_bulk_p = soc_bulk_p;

output.theta_s_n = theta_s_n;
output.theta_s_p = theta_s_p;

output.cs_p = cs_p;
output.cs_n = cs_n;

output.csp_surf = theta_s_p_surf .* param.csp_max;
output.csn_surf = theta_s_n_surf .* param.csn_max;

param_out = param;

end