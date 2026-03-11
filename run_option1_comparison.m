%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTION 1 DRIVER: FDM vs FVM vs PADE PERFORMANCE COMPARISON
% Science/Nature style figure generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clearvars;
clc;
close all;
set(groot,'defaultFigureToolBar','none');

addpath('./cell_model');
addpath('./utils');

addpath(genpath('/Users/zhengjieyang/Documents/MATLAB/SPM_skeleton/casadi-3.7.2-osx_arm64-matlab2018b'));

%% USER SETTINGS
profile_file = './data/data_HPPC.mat';
profile_name = 'HPPC';

solver_opts.integrator = 'casadi';
solver_opts.uniform_initial_cond = true;

Nr_list = [4 6 8 12 16 24 32 48 64 101];

method_list = {'FDM','FVM-S0','FVM-S1','FVM-S2','PADE2','PADE3'};

Nr_ref = 101;
compare_nr = 12;

dt = 1.0;
Tend = 4000;

T_amb = 23;
SOC_IC = 1.0;
Q_IC = 4.8768;

output_dir = './master_output_dir/option1';
if exist(output_dir,'dir') ~= 7
    mkdir(output_dir);
end

%% LOAD DATA
[t_data,I_data,V_data,data_opts] = load_data(profile_file,dt,Tend);

%% DENSE FDM REFERENCE
fprintf('Running dense FDM reference (Nr = %d)...\n',Nr_ref);

ref_case = run_single_case('FDM','',Nr_ref,t_data,I_data,V_data,data_opts,...
    T_amb,SOC_IC,Q_IC,solver_opts,profile_name);

V_ref_dense = ref_case.all_data.V;

%% RUN ALL CASES
num_cases = numel(method_list)*numel(Nr_list);

results(num_cases,1)=struct();
metrics(num_cases,1)=struct();

row=0;

for i_method = 1:numel(method_list)

    method_full = method_list{i_method};

    if startsWith(method_full,'FVM')

        method = 'FVM';
        scheme = extractAfter(method_full,'FVM-');

    elseif startsWith(method_full,'PADE')

        method = method_full;
        scheme = '';

    else

        method = 'FDM';
        scheme = '';

    end

    for i_nr = 1:numel(Nr_list)

        Nr = Nr_list(i_nr);
        row=row+1;

        fprintf('Running %s, Nr = %d ...\n',method_full,Nr);

        solver_opts.fvm_scheme = scheme;

        case_result = run_single_case(method,scheme,Nr,t_data,I_data,V_data,data_opts,...
            T_amb,SOC_IC,Q_IC,solver_opts,profile_name);

        inv = compute_lithium_inventory(case_result.all_data,method);
        rel_inv_drift = (inv-inv(1))./inv(1);

        metrics(row).method = method_full;
        metrics(row).Nr = Nr;
        metrics(row).states_per_electrode = case_result.solver_opts.nstates_electrode;
        metrics(row).total_states = case_result.solver_opts.nstates_tot;
        metrics(row).runtime_s = case_result.runtime_s;

        metrics(row).rmse_vs_exp_V = sqrt(mean((case_result.all_data.V - case_result.all_data.V_ref).^2));
        metrics(row).rmse_vs_dense_fdm_V = sqrt(mean((case_result.all_data.V - V_ref_dense).^2));
        metrics(row).max_rel_inventory_drift = max(abs(rel_inv_drift));

        results(row).case_result = case_result;

    end
end

%% OBSERVABILITY SUMMARY ACROSS ALL CASES
obs_sample_count = 7;
observability_sweep = analyze_observability_sweep(results, metrics, obs_sample_count);

for k = 1:numel(metrics)

    metrics(k).obs_sample_count = observability_sweep(k).sample_count;
    metrics(k).obs_min_rank = observability_sweep(k).min_rank;
    metrics(k).obs_max_rank = observability_sweep(k).max_rank;
    metrics(k).obs_min_rank_fraction = observability_sweep(k).min_rank_fraction;
    metrics(k).obs_finite_cond_fraction = observability_sweep(k).finite_condition_fraction;
    metrics(k).obs_median_finite_condition_number = ...
        observability_sweep(k).median_finite_condition_number;
    metrics(k).obs_min_sigma_min = observability_sweep(k).min_sigma_min;
    metrics(k).obs_min_ocp_slope_metric = observability_sweep(k).min_ocp_slope_metric;

end

metrics_table = struct2table(metrics);
observability_sweep_table = struct2table(observability_sweep);

%% OBSERVABILITY ANALYSIS FOR REPRESENTATIVE CASES
[observability_results, observability_summary] = ...
    analyze_selected_observability(results, metrics, compare_nr);

observability_summary_table = struct2table(observability_summary);

%% SAVE METRICS
save(fullfile(output_dir,'option1_metrics_and_results.mat'),...
    'metrics_table','metrics','results','ref_case', ...
    'observability_results','observability_summary_table', ...
    'observability_sweep_table');

writetable(metrics_table,fullfile(output_dir,'option1_metrics.csv'));
writetable(observability_summary_table, ...
    fullfile(output_dir,'option1_observability_summary.csv'));
writetable(observability_sweep_table, ...
    fullfile(output_dir,'option1_observability_sweep_summary.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% COLOR PALETTE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

color_map = containers.Map( ...
    {'FDM','FVM-S0','FVM-S1','FVM-S2','PADE2','PADE3'}, ...
    { [0.00 0.45 0.74], ...
      [0.85 0.33 0.10], ...
      [0.93 0.69 0.13], ...
      [0.49 0.18 0.56], ...
      [0.30 0.75 0.93], ...
      [0.64 0.08 0.18] });

set_plot_style = @(ax) set(ax,...
    'FontSize',16,...
    'LineWidth',1.3,...
    'Box','on',...
    'GridColor',[0.85 0.85 0.85]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VOLTAGE VS TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',[200 200 900 450])
hold on
grid on

set_plot_style(gca)

plot(t_data,V_data,'k','LineWidth',3,'DisplayName','Experimental')

for k=1:numel(results)

    this_case = results(k).case_result;

    if this_case.Nr == compare_nr

        method_label = metrics(k).method;
        color = color_map(method_label);

        plot(this_case.all_data.t,...
             this_case.all_data.V,...
             'LineWidth',2.5,...
             'Color',color,...
             'DisplayName',method_label)

    end
end

xlabel('Time (s)','FontSize',18,'FontWeight','bold')
ylabel('Voltage (V)','FontSize',18,'FontWeight','bold')

title('Voltage vs Time','FontSize',18,'FontWeight','bold','Color','k')

legend('Location','best','FontSize',14)

exportgraphics(gcf,fullfile(output_dir,'fig_voltage_vs_time.png'),'Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RMSE VS EXPERIMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',[200 200 600 450])
hold on
grid on

set_plot_style(gca)

plot_metric_vs_states(metrics_table,'FDM','rmse_vs_exp_V','o-')
plot_metric_vs_states(metrics_table,'FVM-S0','rmse_vs_exp_V','s-')
plot_metric_vs_states(metrics_table,'FVM-S1','rmse_vs_exp_V','d-')
plot_metric_vs_states(metrics_table,'FVM-S2','rmse_vs_exp_V','^-')
plot_metric_vs_states(metrics_table,'PADE2','rmse_vs_exp_V','v-')
plot_metric_vs_states(metrics_table,'PADE3','rmse_vs_exp_V','>-')

set(gca,'XScale','log')
xtickformat('%.0f')

xlabel('Total States','FontSize',18,'FontWeight','bold')
ylabel('RMSE vs Experiment (V)','FontSize',18,'FontWeight','bold')

title('Model Accuracy','FontSize',18,'FontWeight','bold','Color','k')

legend('Location','best')

exportgraphics(gcf,fullfile(output_dir,'fig_rmse_vs_experiment.png'),'Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RMSE VS DENSE FDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',[200 200 600 450])
hold on
grid on

set_plot_style(gca)

plot_metric_vs_states(metrics_table,'FDM','rmse_vs_dense_fdm_V','o-')
plot_metric_vs_states(metrics_table,'FVM-S0','rmse_vs_dense_fdm_V','s-')
plot_metric_vs_states(metrics_table,'FVM-S1','rmse_vs_dense_fdm_V','d-')
plot_metric_vs_states(metrics_table,'FVM-S2','rmse_vs_dense_fdm_V','^-')
plot_metric_vs_states(metrics_table,'PADE2','rmse_vs_dense_fdm_V','v-')
plot_metric_vs_states(metrics_table,'PADE3','rmse_vs_dense_fdm_V','>-')

set(gca,'XScale','log')
xtickformat('%.0f')

xlabel('Total States','FontSize',18,'FontWeight','bold')
ylabel('RMSE vs Dense FDM (V)','FontSize',18,'FontWeight','bold')

title('Numerical Convergence','FontSize',18,'FontWeight','bold','Color','k')

legend('Location','best')

exportgraphics(gcf,fullfile(output_dir,'fig_rmse_vs_dense_fdm.png'),'Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUNTIME VS STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',[200 200 600 450])
hold on
grid on

set_plot_style(gca)

plot_metric_vs_states(metrics_table,'FDM','runtime_s','o-')
plot_metric_vs_states(metrics_table,'FVM-S0','runtime_s','s-')
plot_metric_vs_states(metrics_table,'FVM-S1','runtime_s','d-')
plot_metric_vs_states(metrics_table,'FVM-S2','runtime_s','^-')
plot_metric_vs_states(metrics_table,'PADE2','runtime_s','v-')
plot_metric_vs_states(metrics_table,'PADE3','runtime_s','>-')

set(gca,'XScale','log')
xtickformat('%.0f')

xlabel('Total States','FontSize',18,'FontWeight','bold')
ylabel('Runtime (s)','FontSize',18,'FontWeight','bold')

title('Runtime Scaling','FontSize',18,'FontWeight','bold','Color','k')

legend('Location','best')

exportgraphics(gcf,fullfile(output_dir,'fig_runtime_scaling.png'),'Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LITHIUM CONSERVATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',[200 200 600 450])
hold on
grid on

set_plot_style(gca)

plot_metric_vs_states(metrics_table,'FDM','max_rel_inventory_drift','o-')
plot_metric_vs_states(metrics_table,'FVM-S0','max_rel_inventory_drift','s-')
plot_metric_vs_states(metrics_table,'FVM-S1','max_rel_inventory_drift','d-')
plot_metric_vs_states(metrics_table,'FVM-S2','max_rel_inventory_drift','^-')
plot_metric_vs_states(metrics_table,'PADE2','max_rel_inventory_drift','v-')
plot_metric_vs_states(metrics_table,'PADE3','max_rel_inventory_drift','>-')

set(gca,'XScale','log')
xtickformat('%.0f')

xlabel('Total States','FontSize',18,'FontWeight','bold')
ylabel('Inventory Drift','FontSize',18,'FontWeight','bold')

title('Lithium Conservation','FontSize',18,'FontWeight','bold','Color','k')

legend('Location','best')

exportgraphics(gcf,fullfile(output_dir,'fig_lithium_conservation.png'),'Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ACCURACY VS RUNTIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',[200 200 600 450])
hold on
grid on

set_plot_style(gca)

scatter(metrics_table.runtime_s,...
        metrics_table.rmse_vs_dense_fdm_V,...
        90,'filled')

xlabel('Runtime (s)','FontSize',18,'FontWeight','bold')
ylabel('RMSE vs Dense FDM (V)','FontSize',18,'FontWeight','bold')

title('Accuracy vs Runtime Tradeoff','FontSize',18,'FontWeight','bold','Color','k')

exportgraphics(gcf,fullfile(output_dir,'fig_accuracy_vs_runtime.png'),'Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OBSERVABILITY SUMMARY VS STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Color','w','Position',[150 120 700 980])
tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile
hold on
grid on
set_plot_style(gca)

plot_metric_vs_states(observability_sweep_table,'FDM','min_rank_fraction','o-')
plot_metric_vs_states(observability_sweep_table,'FVM-S0','min_rank_fraction','s-')
plot_metric_vs_states(observability_sweep_table,'FVM-S1','min_rank_fraction','d-')
plot_metric_vs_states(observability_sweep_table,'FVM-S2','min_rank_fraction','^-')
plot_metric_vs_states(observability_sweep_table,'PADE2','min_rank_fraction','v-')
plot_metric_vs_states(observability_sweep_table,'PADE3','min_rank_fraction','>-')

set(gca,'XScale','log')
xtickformat('%.0f')
ylabel('Min rank / n','FontSize',16,'FontWeight','bold')
title(sprintf('Sampled Observability Summary (%d SOC samples per case)',obs_sample_count),...
    'FontSize',18,'FontWeight','bold','Color','k')
legend('Location','bestoutside')

nexttile
hold on
grid on
set_plot_style(gca)

plot_metric_vs_states(observability_sweep_table,'FDM','finite_condition_fraction','o-')
plot_metric_vs_states(observability_sweep_table,'FVM-S0','finite_condition_fraction','s-')
plot_metric_vs_states(observability_sweep_table,'FVM-S1','finite_condition_fraction','d-')
plot_metric_vs_states(observability_sweep_table,'FVM-S2','finite_condition_fraction','^-')
plot_metric_vs_states(observability_sweep_table,'PADE2','finite_condition_fraction','v-')
plot_metric_vs_states(observability_sweep_table,'PADE3','finite_condition_fraction','>-')

set(gca,'XScale','log')
xtickformat('%.0f')
ylabel('Finite cond fraction','FontSize',16,'FontWeight','bold')

nexttile
hold on
grid on
set_plot_style(gca)

plot_metric_vs_states(observability_sweep_table,'FDM','median_finite_condition_number','o-')
plot_metric_vs_states(observability_sweep_table,'FVM-S0','median_finite_condition_number','s-')
plot_metric_vs_states(observability_sweep_table,'FVM-S1','median_finite_condition_number','d-')
plot_metric_vs_states(observability_sweep_table,'FVM-S2','median_finite_condition_number','^-')
plot_metric_vs_states(observability_sweep_table,'PADE2','median_finite_condition_number','v-')
plot_metric_vs_states(observability_sweep_table,'PADE3','median_finite_condition_number','>-')

set(gca,'XScale','log','YScale','log')
xtickformat('%.0f')
xlabel('Total States','FontSize',16,'FontWeight','bold')
ylabel('$\mathrm{Median\ finite}\ \mathrm{cond}(\mathcal{O})$',...
    'Interpreter','latex','FontSize',16,'FontWeight','bold')

title(tl,'Observability Across Model Orders','FontSize',18,'FontWeight','bold')

exportgraphics(gcf,fullfile(output_dir,'fig_observability_summary_vs_states.png'),'Resolution',300)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OBSERVABILITY VS SOC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(observability_results)

    figure('Color','w','Position',[150 120 820 860])
    tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

    nexttile
    hold on
    grid on
    set_plot_style(gca)

    for k = 1:numel(observability_results)

        obs_case = observability_results(k);
        color = color_map(obs_case.method);

        plot(obs_case.soc_ref,...
             obs_case.ocp_slope_metric,...
             'LineWidth',2.2,...
             'Color',color,...
             'DisplayName',obs_case.method)

    end

    set(gca,'XDir','reverse')
    ylabel('|dU_p/d\theta_p - dU_n/d\theta_n|','FontSize',16,'FontWeight','bold')
    title(sprintf('OCP Slope and Observability at Nr = %d',compare_nr),...
        'FontSize',18,'FontWeight','bold','Color','k')
    legend('Location','bestoutside')

    nexttile
    hold on
    grid on
    set_plot_style(gca)

    for k = 1:numel(observability_results)

        obs_case = observability_results(k);
        color = color_map(obs_case.method);

        plot(obs_case.soc_ref,...
             obs_case.rank,...
             'LineWidth',2.2,...
             'Color',color,...
             'DisplayName',obs_case.method)

    end

    set(gca,'XDir','reverse')
    ylabel('$\mathrm{rank}(\mathcal{O})$','Interpreter','latex','FontSize',16,'FontWeight','bold')

    nexttile
    hold on
    grid on
    set_plot_style(gca)

    for k = 1:numel(observability_results)

        obs_case = observability_results(k);
        color = color_map(obs_case.method);
        cond_to_plot = obs_case.cond_num;
        cond_to_plot(~isfinite(cond_to_plot)) = NaN;

        semilogy(obs_case.soc_ref,...
                 cond_to_plot,...
                 'LineWidth',2.2,...
                 'Color',color,...
                 'DisplayName',obs_case.method)

    end

    set(gca,'XDir','reverse')
    xlabel('SOC [-]','FontSize',16,'FontWeight','bold')
    ylabel('$\mathrm{cond}(\mathcal{O})$','Interpreter','latex','FontSize',16,'FontWeight','bold')

    title(tl,'Representative Observability Comparison','FontSize',18,'FontWeight','bold')

    exportgraphics(gcf,fullfile(output_dir,'fig_observability_vs_soc.png'),'Resolution',300)

end

fprintf('\nAll scientific-style figures exported.\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS
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

solver_opts.nstates_electrode = n_states(solver_opts);
solver_opts.nstates_tot = 2*solver_opts.nstates_electrode;

params = modelParameters(Nr);

params.SOC_IC = SOC_IC;
params.Q_IC = Q_IC;
params.method = method;

tic
[sol,t_eval,param_out] = SPM_sim(params,t_data,I_data,T_amb,solver_opts);
runtime_s = toc;

all_data = post_process(sol,param_out,t_eval,solver_opts);

case_result.method = method;
case_result.Nr = Nr;
case_result.runtime_s = runtime_s;
case_result.sol = sol;
case_result.t_eval = t_eval;
case_result.param_out = param_out;
case_result.all_data = all_data;
case_result.solver_opts = solver_opts;

end

function plot_metric_vs_states(T,method,field_name,line_spec)

rows = strcmpi(T.method,method);
Tm = sortrows(T(rows,:), 'total_states');

plot(Tm.total_states,Tm.(field_name),line_spec,...
    'LineWidth',2,'MarkerSize',7,'DisplayName',method)

end

function [obs_results, obs_summary] = analyze_selected_observability(results, metrics, compare_nr)

method_names = unique({metrics.method},'stable');

obs_results = struct([]);
obs_summary = struct([]);

for i_method = 1:numel(method_names)

    method_label = method_names{i_method};

    idx = find(strcmpi({metrics.method},method_label) & [metrics.Nr] == compare_nr, 1, 'first');

    if isempty(idx)
        continue
    end

    obs_case = compute_case_observability(results(idx).case_result, method_label, []);

    obs_results = [obs_results; obs_case]; %#ok<AGROW>

    obs_summary(end+1,1).method = method_label; %#ok<AGROW>
    obs_summary(end).Nr = obs_case.Nr;
    obs_summary(end).total_states = obs_case.total_states;
    obs_summary(end).min_rank = min(obs_case.rank);
    obs_summary(end).max_rank = max(obs_case.rank);
    obs_summary(end).full_rank_fraction = mean(obs_case.rank == obs_case.total_states);
    obs_summary(end).min_sigma_min = min(obs_case.sigma_min);
    finite_cond = obs_case.cond_num(isfinite(obs_case.cond_num));
    if isempty(finite_cond)
        obs_summary(end).median_condition_number = Inf;
    else
        obs_summary(end).median_condition_number = median(finite_cond);
    end
    obs_summary(end).max_condition_number = max(obs_case.cond_num);

end

end

function obs_sweep = analyze_observability_sweep(results, metrics, sample_count)

num_cases = numel(results);
obs_sweep(num_cases,1) = struct();

for i_case = 1:num_cases

    method_label = metrics(i_case).method;
    sample_idx = pick_sample_indices(results(i_case).case_result.all_data.t, sample_count);
    obs_case = compute_case_observability(results(i_case).case_result, method_label, sample_idx);

    finite_cond = obs_case.cond_num(isfinite(obs_case.cond_num));

    obs_sweep(i_case).method = method_label;
    obs_sweep(i_case).Nr = metrics(i_case).Nr;
    obs_sweep(i_case).states_per_electrode = metrics(i_case).states_per_electrode;
    obs_sweep(i_case).total_states = metrics(i_case).total_states;
    obs_sweep(i_case).sample_count = numel(sample_idx);
    obs_sweep(i_case).min_rank = min(obs_case.rank);
    obs_sweep(i_case).max_rank = max(obs_case.rank);
    obs_sweep(i_case).min_rank_fraction = min(obs_case.rank) / obs_case.total_states;
    obs_sweep(i_case).finite_condition_fraction = mean(isfinite(obs_case.cond_num));
    obs_sweep(i_case).min_sigma_min = min(obs_case.sigma_min);
    obs_sweep(i_case).min_ocp_slope_metric = min(obs_case.ocp_slope_metric);

    if isempty(finite_cond)
        obs_sweep(i_case).median_finite_condition_number = Inf;
        obs_sweep(i_case).max_finite_condition_number = Inf;
    else
        obs_sweep(i_case).median_finite_condition_number = median(finite_cond);
        obs_sweep(i_case).max_finite_condition_number = max(finite_cond);
    end

end

end

function obs_case = compute_case_observability(case_result, method_label, sample_idx)

param = case_result.param_out;
all_data = case_result.all_data;

A = param.A_mat;
n = size(A,1);

[surf_map_n, surf_map_p] = build_surface_maps(case_result);

if isempty(sample_idx)
    sample_idx = (1:numel(all_data.t)).';
else
    sample_idx = sample_idx(:);
end

num_t = numel(sample_idx);

rank_hist = zeros(num_t,1);
cond_hist = zeros(num_t,1);
sigma_min_hist = zeros(num_t,1);
ocp_slope_metric = zeros(num_t,1);

for kk = 1:num_t

    k = sample_idx(kk);

    I_k = all_data.I(k);

    theta_n = clamp_theta(all_data.cssurn(k) / param.csn_max);
    theta_p = clamp_theta(all_data.cssurp(k) / param.csp_max);

    dUn_dtheta = finite_difference_scalar(@(theta) U_n(theta,param.T_amb,param), theta_n);
    dUp_dtheta = finite_difference_scalar(@(theta) U_p(theta,param.T_amb,param), theta_p);

    dEtaN_dtheta = finite_difference_scalar( ...
        @(theta) eta_anode(theta,param.T_amb,I_k,param), theta_n);

    dEtaP_dtheta = finite_difference_scalar( ...
        @(theta) eta_cathode(theta,param.T_amb,I_k,param), theta_p);

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

    singular_vals = svd(O);
    tol = max(size(O)) * eps(max(singular_vals));

    rank_hist(kk) = sum(singular_vals > tol);
    sigma_min_hist(kk) = singular_vals(end);

    if singular_vals(end) <= tol
        cond_hist(kk) = Inf;
    else
        cond_hist(kk) = singular_vals(1) / singular_vals(end);
    end

    ocp_slope_metric(kk) = abs(dUp_dtheta - dUn_dtheta);

end

obs_case.method = method_label;
obs_case.Nr = case_result.Nr;
obs_case.total_states = n;
obs_case.t = all_data.t(sample_idx);
obs_case.soc_ref = all_data.SOC_ref(sample_idx);
obs_case.rank = rank_hist;
obs_case.cond_num = cond_hist;
obs_case.sigma_min = sigma_min_hist;
obs_case.ocp_slope_metric = ocp_slope_metric;

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

    otherwise

        error('Unknown method for observability analysis')

end

end

function sample_idx = pick_sample_indices(t_vec, sample_count)

if numel(t_vec) <= sample_count
    sample_idx = (1:numel(t_vec)).';
else
    sample_idx = unique(round(linspace(1,numel(t_vec),sample_count))).';
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