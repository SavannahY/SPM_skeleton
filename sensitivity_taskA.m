%% Task A – Local Sensitivity Analysis (HPPC)
% ENERGY 295 – Electrochemical Energy Storage Systems
% Local (±10%) sensitivity analysis using HPPC profile
%
% Uses correct post_process outputs for this codebase:
%   voltage field = out.V
% NOTE: CasADi may warn about internal evaluation times; verified full output grid is returned.


restoredefaultpath; clear; clc; close all;

%% --- Paths ---
thisdir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisdir,'cell_model'));
addpath(fullfile(thisdir,'utils'));
addpath(genpath('/Users/zhengjieyang/Documents/MATLAB/SPM_skeleton/casadi-3.7.2-osx_arm64-matlab2018b'));

%% --- Load HPPC data ---
dt   = 1.0;     % [s]
Tend = inf;     % simulate full profile
[data_t, I_data, V_data, solver_opts] = ...
    load_data(fullfile(thisdir,'data','data_HPPC.mat'), dt, Tend);

%% --- Solver options ---
solver_opts.method = 'FDM';
solver_opts.integrator = 'casadi';
solver_opts.uniform_initial_cond = true;

Nr = 101;
solver_opts.Nr = Nr;
solver_opts.nstates_electrode = n_states(solver_opts);
solver_opts.nstates_tot = 2*solver_opts.nstates_electrode;

%% --- Nominal parameters ---
param0 = modelParameters(Nr);
param0.SOC_IC = 1.0;
param0.Q_IC   = 4.8768;    % [Ah]
T_amb = 23;                % [C]

%% --- Baseline simulation ---
[sol0, t_eval, p0] = SPM_sim(param0, data_t, I_data, T_amb, solver_opts);
out0 = post_process(sol0, p0, t_eval, solver_opts);

Vnom = out0.V(:);                  % <-- correct field
Vnom(abs(Vnom) < 1e-6) = 1e-6;      % avoid divide-by-zero
t_base = out0.t(:);
SOC_base = out0.soc_bulk_n(:);      % use anode SOC for SOC-dependent plots

%% --- Sensitivity setup ---
p = 0.10;   % ±10%
param_list = {'Dsn','Dsp','k0n','k0p','epsn','epsp','R0'};

S_time = struct();
S_max  = zeros(numel(param_list),1);
S_rms  = zeros(numel(param_list),1);

%% --- Sensitivity loop (Task A) ---
for k = 1:numel(param_list)
    pname = param_list{k};
    fprintf('\n=== Sensitivity for %s ===\n', pname);

    % --- ±10% perturbations (structure-aware) ---
    if ismember(pname, {'epsn','epsp'})
        % Rebuild model to refresh FDM matrices
        pm = modelParameters(Nr);
        pp = modelParameters(Nr);
    
        % Copy ICs and any non-geometry settings
        pm.SOC_IC = param0.SOC_IC;
        pm.Q_IC   = param0.Q_IC;
        pp.SOC_IC = param0.SOC_IC;
        pp.Q_IC   = param0.Q_IC;
    
        % Apply perturbation
        pm.(pname) = (1-p)*param0.(pname);
        pp.(pname) = (1+p)*param0.(pname);
    else
        % Safe to perturb algebraic parameters directly
        pm = param0;
        pp = param0;
    
        pm.(pname) = (1-p)*param0.(pname);
        pp.(pname) = (1+p)*param0.(pname);
    end


    % Run simulations
    [solm, tm, pm_out] = SPM_sim(pm, data_t, I_data, T_amb, solver_opts);
    [solp, tp, pp_out] = SPM_sim(pp, data_t, I_data, T_amb, solver_opts);

    outm = post_process(solm, pm_out, tm, solver_opts);
    outp = post_process(solp, pp_out, tp, solver_opts);

    Vm = outm.V(:);
    Vp = outp.V(:);

    % Align lengths (in case solver truncates)
    N = min([numel(Vnom), numel(Vm), numel(Vp)]);
    Vn = Vnom(1:N);
    Vm = Vm(1:N);
    Vp = Vp(1:N);
    tt = t_base(1:N);
    SOC = SOC_base(1:N);

    % --- Normalized sensitivity (central difference) ---
    % S_theta(t) = (V+ - V-) / (0.2 * V_nom)
    S = (Vp - Vm) ./ (2*p*Vn);

    % Store
    S_time.(pname).t   = tt;
    S_time.(pname).SOC = SOC;
    S_time.(pname).S   = S;

    S_max(k) = max(abs(S));
    S_rms(k) = sqrt(mean(S.^2));

    %% --- Plots ---
    % Voltage overlay
    figure;
    plot(tt, Vn, 'k', tt, Vm, '--b', tt, Vp, '--r', 'LineWidth', 1.2);
    grid on;
    xlabel('Time [s]');
    ylabel('Voltage [V]');
    legend('Nominal','-10%','+10%','Location','best');
    title(['HPPC Voltage Overlay – ' pname]);

    % Sensitivity vs time
    figure;
    plot(tt, S, 'LineWidth', 1.5);
    grid on;
    xlabel('Time [s]');
    ylabel('S_\theta(t)');
    title(['Normalized Sensitivity vs Time – ' pname]);

    % Sensitivity vs SOC
    figure;
    plot(SOC, abs(S), '.', 'MarkerSize', 6);
    set(gca,'XDir','reverse');
    grid on;
    xlabel('SOC [-]');
    ylabel('|S_\theta|');
    title(['SOC-Dependent Sensitivity – ' pname]);
end

%% --- Ranking table & bar chart ---
T = table(string(param_list(:)), S_max, S_rms, ...
    'VariableNames',{'Parameter','Smax','Srms'});
T = sortrows(T,'Srms','descend');
disp(T);

figure;
bar(T.Srms);
set(gca,'XTickLabel',T.Parameter);
xtickangle(45);
ylabel('Mean |S_\theta|');
title('Task A – Parameter Sensitivity Ranking (HPPC)');
grid on;

%% --- Analytic sanity check for R0 ---
if ismember('R0',param_list)
    idx = find(strcmp(param_list,'R0'));
    S_R0_num = S_time.R0.S;
    I_use = I_data(1:numel(S_R0_num));
    S_R0_analytic = -(I_use*param0.R0) ./ Vnom(1:numel(S_R0_num));

    fprintf('Max error in R0 analytic sensitivity = %.3e\n', ...
        max(abs(S_R0_num - S_R0_analytic)));
end


%% Export data for plotting (non-MATLAB figures)
baseline.t   = t_base(:);
baseline.V   = Vnom(:);

% Make I same length as baseline time/voltage
N0 = numel(baseline.t);
baseline.I   = I_data(1:N0);
baseline.SOC = SOC_base(1:N0);

save('sensitivity_export.mat', 'baseline', 'S_time', 'T');

