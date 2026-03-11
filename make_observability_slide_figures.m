%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE_OBSERVABILITY_SLIDE_FIGURES
% Generate slide-ready two-panel figures for the observability section.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

restoredefaultpath;
clearvars;
clc;
close all;
set(groot,'defaultFigureToolBar','none');

addpath('./cell_model');
addpath('./utils');

output_dir = './master_output_dir/option1';
slide_dir = fullfile(output_dir,'slide_deck_figures');

if exist(slide_dir,'dir') ~= 7
    mkdir(slide_dir);
end

data = load(fullfile(output_dir,'option1_metrics_and_results.mat'));

metrics_table = data.metrics_table;
observability_results = data.observability_results;
observability_sweep_table = data.observability_sweep_table;

param = modelParameters(101);
param.T_amb = 273.15 + 23;

color_map = containers.Map( ...
    {'FDM','FVM-S0','FVM-S1','FVM-S2','PADE2','PADE3'}, ...
    { [0.00 0.45 0.74], ...
      [0.85 0.33 0.10], ...
      [0.93 0.69 0.13], ...
      [0.49 0.18 0.56], ...
      [0.30 0.75 0.93], ...
      [0.64 0.08 0.18] });

marker_map = containers.Map( ...
    {'FDM','FVM-S0','FVM-S1','FVM-S2','PADE2','PADE3'}, ...
    {'o','s','d','^','v','>'});

focus_methods = {'FDM','FVM-S2','PADE2','PADE3'};

set_plot_style = @(ax) set(ax,...
    'FontName','Arial',...
    'FontSize',16,...
    'LineWidth',1.4,...
    'Box','on',...
    'TickDir','out',...
    'XColor',[0.10 0.10 0.10],...
    'YColor',[0.10 0.10 0.10],...
    'GridColor',[0.82 0.82 0.82],...
    'GridAlpha',0.35);

%% Slide 1: OCP and slope sensitivity
soc_grid = linspace(0,1,500)';
theta_n_grid = param.theta_n_0 + soc_grid .* (param.theta_n_100 - param.theta_n_0);
theta_p_grid = param.theta_p_0 + soc_grid .* (param.theta_p_100 - param.theta_p_0);

U_n_grid = arrayfun(@(x) U_n(x,param.T_amb,param), theta_n_grid);
U_p_grid = arrayfun(@(x) U_p(x,param.T_amb,param), theta_p_grid);
Voc_grid = U_p_grid - U_n_grid;

dU_n = numerical_derivative(@(x) U_n(x,param.T_amb,param), theta_n_grid);
dU_p = numerical_derivative(@(x) U_p(x,param.T_amb,param), theta_p_grid);
slope_metric = abs(dU_p - dU_n);

fig = figure('Color','w','Position',[100 120 1600 700],'ToolBar','none','MenuBar','none');
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
hold(ax1,'on');
grid(ax1,'on');
set_plot_style(ax1);
plot(ax1,soc_grid,Voc_grid,'k-','LineWidth',3.0,'DisplayName','Cell OCV = U_p - U_n');
plot(ax1,soc_grid,U_p_grid,'-','LineWidth',2.6,'Color',[0.75 0.15 0.12],'DisplayName','Cathode OCP, U_p');
plot(ax1,soc_grid,U_n_grid,'-','LineWidth',2.6,'Color',[0.00 0.35 0.70],'DisplayName','Anode OCP, U_n');
xlabel(ax1,'SOC [-]','FontWeight','bold');
ylabel(ax1,'Voltage [V]','FontWeight','bold');
title(ax1,'OCP Curves Across SOC','FontSize',20,'FontWeight','bold');
legend(ax1,'Location','best','FontSize',14);

ax2 = nexttile;
hold(ax2,'on');
grid(ax2,'on');
set_plot_style(ax2);
plot(ax2,soc_grid,abs(dU_p),'-','LineWidth',2.6,'Color',[0.75 0.15 0.12],'DisplayName','|dU_p / d\theta_p|');
plot(ax2,soc_grid,abs(dU_n),'-','LineWidth',2.6,'Color',[0.00 0.35 0.70],'DisplayName','|dU_n / d\theta_n|');
plot(ax2,soc_grid,slope_metric,'k-','LineWidth',3.0,'DisplayName','|dU_p/d\theta_p - dU_n/d\theta_n|');
xlabel(ax2,'SOC [-]','FontWeight','bold');
ylabel(ax2,'Slope Magnitude [V]','FontWeight','bold');
title(ax2,'Voltage Sensitivity to Stoichiometry','FontSize',20,'FontWeight','bold');
legend(ax2,'Location','best','FontSize',14);

title(tl,'Slide 1: Why Observability Changes With SOC','FontSize',22,'FontWeight','bold');
exportgraphics(fig,fullfile(slide_dir,'slide1_ocp_and_slope.png'),'Resolution',450);

%% Slide 2: Representative observability vs SOC
fig = figure('Color','w','Position',[100 120 1600 700],'ToolBar','none','MenuBar','none');
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
hold(ax1,'on');
grid(ax1,'on');
set_plot_style(ax1);

ax2 = nexttile;
hold(ax2,'on');
grid(ax2,'on');
set_plot_style(ax2);

for i = 1:numel(focus_methods)
    method = focus_methods{i};
    obs_case = get_obs_case(observability_results, method);
    color = color_map(method);

    plot(ax1,obs_case.soc_ref,obs_case.rank,...
        'LineWidth',2.8,'Color',color,'DisplayName',method);

    cond_plot = obs_case.cond_num;
    cond_plot(~isfinite(cond_plot)) = NaN;
    semilogy(ax2,obs_case.soc_ref,cond_plot,...
        'LineWidth',2.8,'Color',color,'DisplayName',method);
end

set(ax1,'XDir','reverse');
set(ax2,'XDir','reverse');
xlabel(ax1,'SOC [-]','FontWeight','bold');
ylabel(ax1,'$\mathrm{rank}(\mathcal{O})$','Interpreter','latex','FontWeight','bold');
title(ax1,'Observability Rank vs SOC','FontSize',20,'FontWeight','bold');
legend(ax1,'Location','best','FontSize',14);

xlabel(ax2,'SOC [-]','FontWeight','bold');
ylabel(ax2,'$\mathrm{cond}(\mathcal{O})$','Interpreter','latex','FontWeight','bold');
title(ax2,'Condition Number vs SOC','FontSize',20,'FontWeight','bold');
legend(ax2,'Location','best','FontSize',14);

title(tl,'Slide 2: Representative Local Observability (Nr = 12)','FontSize',22,'FontWeight','bold');
exportgraphics(fig,fullfile(slide_dir,'slide2_local_observability.png'),'Resolution',450);

%% Slide 3: Observability vs model order
fig = figure('Color','w','Position',[100 120 1600 700],'ToolBar','none','MenuBar','none');
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
hold(ax1,'on');
grid(ax1,'on');
set_plot_style(ax1);

ax2 = nexttile;
hold(ax2,'on');
grid(ax2,'on');
set_plot_style(ax2);

for i = 1:numel(focus_methods)
    method = focus_methods{i};
    rows = strcmpi(observability_sweep_table.method,method);
    Tm = sortrows(observability_sweep_table(rows,:), 'total_states');
    color = color_map(method);
    marker = marker_map(method);

    plot(ax1,Tm.total_states,Tm.min_rank_fraction,...
        '-','LineWidth',2.6,'Color',color,...
        'Marker',marker,'MarkerSize',8,'MarkerFaceColor',color,...
        'DisplayName',method);

    semilogy(ax2,Tm.total_states,max(Tm.min_sigma_min,realmin),...
        '-','LineWidth',2.6,'Color',color,...
        'Marker',marker,'MarkerSize',8,'MarkerFaceColor',color,...
        'DisplayName',method);
end

set(ax1,'XScale','log');
set(ax2,'XScale','log');
xtickformat(ax1,'%.0f');
xtickformat(ax2,'%.0f');
xlabel(ax1,'Total States','FontWeight','bold');
ylabel(ax1,'Minimum rank / n','FontWeight','bold');
title(ax1,'Normalized Rank vs Model Order','FontSize',20,'FontWeight','bold');
legend(ax1,'Location','best','FontSize',14);

xlabel(ax2,'Total States','FontWeight','bold');
ylabel(ax2,'$\mathrm{Minimum}\ \sigma_{\min}(\mathcal{O})$','Interpreter','latex','FontWeight','bold');
title(ax2,'Weakest Observable Direction vs Model Order','FontSize',20,'FontWeight','bold');
legend(ax2,'Location','best','FontSize',14);

title(tl,'Slide 3: Global Observability Trends Across Model Order','FontSize',22,'FontWeight','bold');
exportgraphics(fig,fullfile(slide_dir,'slide3_model_order_observability.png'),'Resolution',450);

%% Slide 4: Accuracy and observability tradeoff
fig = figure('Color','w','Position',[100 120 1600 700],'ToolBar','none','MenuBar','none');
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile;
hold(ax1,'on');
grid(ax1,'on');
set_plot_style(ax1);

ax2 = nexttile;
hold(ax2,'on');
grid(ax2,'on');
set_plot_style(ax2);

for i = 1:numel(focus_methods)
    method = focus_methods{i};
    rows = strcmpi(metrics_table.method,method);
    Tm = sortrows(metrics_table(rows,:), 'total_states');
    color = color_map(method);
    marker = marker_map(method);

    plot(ax1,Tm.total_states,1e3*Tm.rmse_vs_exp_V,...
        '-','LineWidth',2.6,'Color',color,...
        'Marker',marker,'MarkerSize',8,'MarkerFaceColor',color,...
        'DisplayName',method);

    plot(ax2,Tm.obs_min_rank_fraction,1e3*Tm.rmse_vs_exp_V,...
        '-','LineWidth',2.6,'Color',color,...
        'Marker',marker,'MarkerSize',8,'MarkerFaceColor',color,...
        'DisplayName',method);
end

set(ax1,'XScale','log');
xtickformat(ax1,'%.0f');
xlabel(ax1,'Total States','FontWeight','bold');
ylabel(ax1,'Voltage RMSE vs Experiment [mV]','FontWeight','bold');
title(ax1,'Accuracy vs Model Order','FontSize',20,'FontWeight','bold');
legend(ax1,'Location','best','FontSize',14);

xlabel(ax2,'Minimum rank / n','FontWeight','bold');
ylabel(ax2,'Voltage RMSE vs Experiment [mV]','FontWeight','bold');
title(ax2,'Accuracy vs Observability Strength','FontSize',20,'FontWeight','bold');
legend(ax2,'Location','best','FontSize',14);

title(tl,'Slide 4: Practical Tradeoff Between Observability and Accuracy','FontSize',22,'FontWeight','bold');
exportgraphics(fig,fullfile(slide_dir,'slide4_accuracy_observability_tradeoff.png'),'Resolution',450);

fprintf('\nSlide-ready observability figures exported to:\n%s\n\n', slide_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydx = numerical_derivative(fun_handle, x_grid)

dydx = zeros(size(x_grid));
h = 1e-6;

for k = 1:numel(x_grid)
    x0 = min(max(x_grid(k),h),1-h);
    x_lo = max(x0-h,1e-8);
    x_hi = min(x0+h,1-1e-8);
    dydx(k) = (fun_handle(x_hi) - fun_handle(x_lo)) / (x_hi - x_lo);
end

end

function obs_case = get_obs_case(observability_results, method_name)

idx = find(strcmpi({observability_results.method},method_name),1,'first');
if isempty(idx)
    error('Method %s not found in observability results.', method_name)
end

obs_case = observability_results(idx);

end
