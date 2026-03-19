%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORT_OBSERVABILITY_DEEP_DIVE_METRICS
% Recompute representative observability data with full singular-value spectra
% for higher-quality visualization.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clearvars;
clc;
close all;

addpath('./cell_model');
addpath('./utils');
addpath(genpath('./casadi-3.7.2-osx_arm64-matlab2018b'));

output_dir = './master_output_dir/option1';
if exist(output_dir,'dir') ~= 7
    mkdir(output_dir);
end

profile_file = './data/data_HPPC.mat';
profile_name = 'HPPC';
dt = 1.0;
Tend = inf;
T_amb = 23;
SOC_IC = 1.0;
Q_IC = 4.8768;
compare_nr = 12;

solver_opts.integrator = 'casadi';
solver_opts.uniform_initial_cond = true;

[t_data,I_data,V_data,data_opts] = load_data(profile_file,dt,Tend);

method_specs = [ ...
    struct('label','FDM','code','FDM','scheme',''), ...
    struct('label','FVM-S0','code','FVM','scheme','S0'), ...
    struct('label','FVM-S1','code','FVM','scheme','S1'), ...
    struct('label','FVM-S2','code','FVM','scheme','S2'), ...
    struct('label','Nonlinear SPM-Padé 2','code','PADE2','scheme',''), ...
    struct('label','Nonlinear SPM-Padé 3','code','PADE3','scheme',''), ...
    struct('label','Local Linear Padé-ECM','code','LLPADEECM','scheme','') ...
];

deep_metrics = struct([]);

for i_method = 1:numel(method_specs)
    spec = method_specs(i_method);
    fprintf('Recomputing observability spectra for %s ...\n', spec.label);

    case_result = run_single_case(spec.code,spec.scheme,compare_nr,t_data,I_data,V_data,data_opts,...
        T_amb,SOC_IC,Q_IC,solver_opts,profile_name);

    deep_metrics = [deep_metrics; compute_observability_spectra(case_result, spec.label)]; %#ok<AGROW>
end

save(fullfile(output_dir,'observability_deep_dive_metrics.mat'),'deep_metrics','-v7.3');
write_deep_metrics_csv(deep_metrics, fullfile(output_dir,'observability_deep_dive_timeseries.csv'));
plot_pade_normalized_vs_raw(deep_metrics, output_dir);
fprintf('Saved deep-dive observability metrics to %s\n', fullfile(output_dir,'observability_deep_dive_metrics.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function case_result = run_single_case(method,scheme,Nr,t_data,I_data,V_data,data_opts,...
    T_amb,SOC_IC,Q_IC,solver_opts_template,profile_name)

solver_opts = solver_opts_template;
solver_opts.method = method;
solver_opts.Nr = Nr;
solver_opts.dchg_type_str = profile_name;
solver_opts.I_interpolant = data_opts.I_interpolant;
solver_opts.V_interpolant = data_opts.V_interpolant;
solver_opts.fvm_scheme = scheme;
solver_opts.linearization_soc = SOC_IC;

params = modelParameters(Nr);
params.SOC_IC = SOC_IC;
params.Q_IC = Q_IC;
params.method = method;

if strcmpi(method,'LLPADEECM')
    solver_opts.nstates_electrode = 4;
    solver_opts.nstates_tot = 8;
    [sol,t_eval,param_out] = run_local_linear_pade_ecm_casadi(params,t_data,I_data,T_amb,solver_opts);
else
    solver_opts.nstates_electrode = n_states(solver_opts);
    solver_opts.nstates_tot = 2*solver_opts.nstates_electrode;
    [sol,t_eval,param_out] = SPM_sim(params,t_data,I_data,T_amb,solver_opts);
end

all_data = post_process(sol,param_out,t_eval,solver_opts);

case_result.method = method;
case_result.Nr = Nr;
case_result.sol = sol;
case_result.t_eval = t_eval;
case_result.param_out = param_out;
case_result.all_data = all_data;
case_result.solver_opts = solver_opts;
end

function deep = compute_observability_spectra(case_result, method_label)

param = case_result.param_out;
all_data = case_result.all_data;
A = param.A_mat;
n = size(A,1);

[surf_map_n, surf_map_p] = build_surface_maps(case_result);

num_t = numel(all_data.t);
rank_hist = zeros(num_t,1);
sigma_min_hist = zeros(num_t,1);
sigma_min_nz_hist = zeros(num_t,1);
effective_rank_hist = zeros(num_t,1);
rank_raw_hist = zeros(num_t,1);
sigma_min_raw_hist = zeros(num_t,1);
sigma_min_nz_raw_hist = zeros(num_t,1);
effective_rank_raw_hist = zeros(num_t,1);
ocp_slope_metric = zeros(num_t,1);
singular_values = zeros(num_t,n);
singular_values_raw = zeros(num_t,n);

for k = 1:num_t
    I_k = all_data.I(k);
    theta_n = clamp_theta(all_data.cssurn(k) / param.csn_max);
    theta_p = clamp_theta(all_data.cssurp(k) / param.csp_max);

    dUn_dtheta = finite_difference_scalar(@(theta) U_n(theta,param.T_amb,param), theta_n);
    dUp_dtheta = finite_difference_scalar(@(theta) U_p(theta,param.T_amb,param), theta_p);
    dEtaN_dtheta = finite_difference_scalar(@(theta) eta_anode(theta,param.T_amb,I_k,param), theta_n);
    dEtaP_dtheta = finite_difference_scalar(@(theta) eta_cathode(theta,param.T_amb,I_k,param), theta_p);

    dV_dcsn = -(dUn_dtheta + dEtaN_dtheta) / param.csn_max;
    dV_dcsp =  (dUp_dtheta + dEtaP_dtheta) / param.csp_max;
    Ck = dV_dcsn .* surf_map_n + dV_dcsp .* surf_map_p;

    O = zeros(n,n);
    row_vec = Ck;
    O(1,:) = normalize_observability_row(row_vec);
    for i_row = 2:n
        row_vec = O(i_row-1,:) * A;
        O(i_row,:) = normalize_observability_row(row_vec);
    end

    s = svd(O);
    singular_values(k,:) = s(:).';

    tol = max(size(O)) * eps(max(s));
    rank_hist(k) = sum(s > tol);
    sigma_min_hist(k) = s(end);

    nonzero_s = s(s > tol);
    if isempty(nonzero_s)
        sigma_min_nz_hist(k) = 0.0;
    else
        sigma_min_nz_hist(k) = nonzero_s(end);
    end

    if sum(s) <= 0
        effective_rank_hist(k) = 0.0;
    else
        p = s / sum(s);
        p = p(p > 0);
        effective_rank_hist(k) = exp(-sum(p .* log(p)));
    end

    O_raw = zeros(n,n);
    row_vec_raw = full(Ck);
    row_vec_raw(~isfinite(row_vec_raw)) = 0.0;
    O_raw(1,:) = row_vec_raw;
    for i_row = 2:n
        row_vec_raw = O_raw(i_row-1,:) * A;
        row_vec_raw = full(row_vec_raw);
        row_vec_raw(~isfinite(row_vec_raw)) = 0.0;
        O_raw(i_row,:) = row_vec_raw;
    end

    s_raw = svd(O_raw);
    singular_values_raw(k,:) = s_raw(:).';

    tol_raw = max(size(O_raw)) * eps(max(s_raw));
    rank_raw_hist(k) = sum(s_raw > tol_raw);
    sigma_min_raw_hist(k) = s_raw(end);

    nonzero_s_raw = s_raw(s_raw > tol_raw);
    if isempty(nonzero_s_raw)
        sigma_min_nz_raw_hist(k) = 0.0;
    else
        sigma_min_nz_raw_hist(k) = nonzero_s_raw(end);
    end

    if sum(s_raw) <= 0
        effective_rank_raw_hist(k) = 0.0;
    else
        p_raw = s_raw / sum(s_raw);
        p_raw = p_raw(p_raw > 0);
        effective_rank_raw_hist(k) = exp(-sum(p_raw .* log(p_raw)));
    end

    ocp_slope_metric(k) = abs(dUp_dtheta - dUn_dtheta);
end

deep.method = method_label;
deep.Nr = case_result.Nr;
deep.total_states = n;
deep.t = all_data.t(:);
deep.soc_ref = all_data.SOC_ref(:);
deep.rank = rank_hist;
deep.sigma_min = sigma_min_hist;
deep.sigma_min_nonzero = sigma_min_nz_hist;
deep.effective_rank = effective_rank_hist;
deep.rank_raw = rank_raw_hist;
deep.sigma_min_raw = sigma_min_raw_hist;
deep.sigma_min_nonzero_raw = sigma_min_nz_raw_hist;
deep.effective_rank_raw = effective_rank_raw_hist;
deep.ocp_slope_metric = ocp_slope_metric;
deep.singular_values = singular_values;
deep.singular_values_raw = singular_values_raw;
end

function [surf_map_n, surf_map_p] = build_surface_maps(case_result)
n_e = case_result.solver_opts.nstates_electrode;
n_tot = case_result.solver_opts.nstates_tot;

surf_map_n = zeros(1,n_tot);
surf_map_p = zeros(1,n_tot);

switch upper(case_result.method)
    case 'FDM'
        surf_map_n(n_e) = 1.0;
        surf_map_p(end) = 1.0;

    case 'FVM'
        switch upper(case_result.solver_opts.fvm_scheme)
            case 'S0'
                coeffs = 1.0;
                idx_local = n_e;
            case 'S1'
                coeffs = [-0.5, 1.5];
                idx_local = [n_e-1, n_e];
            case 'S2'
                coeffs = [1/6, -5/6, 5/3];
                idx_local = [n_e-2, n_e-1, n_e];
            otherwise
                error('Unknown FVM surface scheme')
        end
        surf_map_n(idx_local) = coeffs;
        surf_map_p(n_e + idx_local) = coeffs;

    case {'PADE2','PADE3'}
        surf_map_n(1:n_e) = full(case_result.param_out.Csurf_n);
        surf_map_p(n_e+1:end) = full(case_result.param_out.Csurf_p);

    case 'LLPADEECM'
        surf_map_n = full(case_result.param_out.Csurf_n);
        surf_map_p = full(case_result.param_out.Csurf_p);

    otherwise
        error('Unknown method for observability analysis')
end
end

function row_out = normalize_observability_row(row_in)
row_out = full(row_in);
row_out(~isfinite(row_out)) = 0.0;
row_scale = norm(row_out,2);
if row_scale > 0
    row_out = row_out ./ row_scale;
end
end

function deriv = finite_difference_scalar(fun_handle, x0)
eps_theta = 1e-6;
x_lo = max(x0 - eps_theta, eps_theta);
x_hi = min(x0 + eps_theta, 1 - eps_theta);
if x_hi <= x_lo
    deriv = 0.0;
    return
end
deriv = (fun_handle(x_hi) - fun_handle(x_lo)) / (x_hi - x_lo);
end

function theta = clamp_theta(theta)
theta = min(max(theta,1e-6),1-1e-6);
end

function write_deep_metrics_csv(deep_metrics, csv_path)

method_col = {};
t_col = [];
soc_col = [];
rank_col = [];
sigma_min_col = [];
sigma_min_nz_col = [];
effective_rank_col = [];
rank_raw_col = [];
sigma_min_raw_col = [];
sigma_min_nz_raw_col = [];
effective_rank_raw_col = [];
ocp_slope_col = [];

for i = 1:numel(deep_metrics)
    dm = deep_metrics(i);
    n = numel(dm.t);

    method_col = [method_col; repmat({dm.method}, n, 1)]; %#ok<AGROW>
    t_col = [t_col; dm.t(:)]; %#ok<AGROW>
    soc_col = [soc_col; dm.soc_ref(:)]; %#ok<AGROW>
    rank_col = [rank_col; dm.rank(:)]; %#ok<AGROW>
    sigma_min_col = [sigma_min_col; dm.sigma_min(:)]; %#ok<AGROW>
    sigma_min_nz_col = [sigma_min_nz_col; dm.sigma_min_nonzero(:)]; %#ok<AGROW>
    effective_rank_col = [effective_rank_col; dm.effective_rank(:)]; %#ok<AGROW>
    rank_raw_col = [rank_raw_col; dm.rank_raw(:)]; %#ok<AGROW>
    sigma_min_raw_col = [sigma_min_raw_col; dm.sigma_min_raw(:)]; %#ok<AGROW>
    sigma_min_nz_raw_col = [sigma_min_nz_raw_col; dm.sigma_min_nonzero_raw(:)]; %#ok<AGROW>
    effective_rank_raw_col = [effective_rank_raw_col; dm.effective_rank_raw(:)]; %#ok<AGROW>
    ocp_slope_col = [ocp_slope_col; dm.ocp_slope_metric(:)]; %#ok<AGROW>
end

T = table(method_col, t_col, soc_col, rank_col, sigma_min_col, sigma_min_nz_col, effective_rank_col, ...
    rank_raw_col, sigma_min_raw_col, sigma_min_nz_raw_col, effective_rank_raw_col, ocp_slope_col, ...
    'VariableNames', {'method','t','soc_ref','rank','sigma_min','sigma_min_nonzero','effective_rank', ...
    'rank_raw','sigma_min_raw','sigma_min_nonzero_raw','effective_rank_raw','ocp_slope_metric'});
writetable(T, csv_path);
end

function plot_pade_normalized_vs_raw(deep_metrics, output_dir)
fig_dir = fullfile(output_dir, 'observability_deep_dive');
if exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

pade_methods = { ...
    'Local Linear Padé-ECM', ...
    'Nonlinear SPM-Padé 2', ...
    'Nonlinear SPM-Padé 3'};
colors = [ ...
    0.00, 0.60, 0.45; ...
    0.30, 0.75, 0.90; ...
    0.60, 0.15, 0.30];

fig = figure('Color', 'w', 'Position', [100, 100, 1200, 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
hold on;
for i = 1:numel(pade_methods)
    dm = select_method(deep_metrics, pade_methods{i});
    plot(dm.soc_ref, dm.rank, '-', 'LineWidth', 2.0, 'Color', colors(i,:));
    plot(dm.soc_ref, dm.rank_raw, '--', 'LineWidth', 2.0, 'Color', colors(i,:));
end
set(gca, 'XDir', 'reverse', 'FontSize', 12, 'LineWidth', 1.0);
xlabel('SOC [-]', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Observability rank', 'FontSize', 16, 'FontWeight', 'bold');
title('Padé Family: Normalized vs Raw Rank', 'FontSize', 20, 'FontWeight', 'bold');
grid on;
legend({ ...
    'Local Linear Padé-ECM (normalized)', ...
    'Local Linear Padé-ECM (raw)', ...
    'Nonlinear SPM-Padé 2 (normalized)', ...
    'Nonlinear SPM-Padé 2 (raw)', ...
    'Nonlinear SPM-Padé 3 (normalized)', ...
    'Nonlinear SPM-Padé 3 (raw)'}, ...
    'Location', 'eastoutside', 'FontSize', 11);

nexttile;
hold on;
for i = 1:numel(pade_methods)
    dm = select_method(deep_metrics, pade_methods{i});
    plot(dm.soc_ref, dm.effective_rank, '-', 'LineWidth', 2.0, 'Color', colors(i,:));
    plot(dm.soc_ref, dm.effective_rank_raw, '--', 'LineWidth', 2.0, 'Color', colors(i,:));
end
set(gca, 'XDir', 'reverse', 'FontSize', 12, 'LineWidth', 1.0);
xlabel('SOC [-]', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Effective rank', 'FontSize', 16, 'FontWeight', 'bold');
title('Padé Family: Normalized vs Raw Effective Rank', 'FontSize', 20, 'FontWeight', 'bold');
grid on;

exportgraphics(fig, fullfile(fig_dir, 'obs_appendix_pade_normalized_vs_raw.png'), 'Resolution', 200);
close(fig);
end

function dm = select_method(deep_metrics, method_label)
idx = find(strcmp({deep_metrics.method}, method_label), 1, 'first');
if isempty(idx)
    error('Method %s not found in deep metrics.', method_label);
end
dm = deep_metrics(idx);
end
