clear; clc;

%% ensure MATLAB sees all folders
addpath(genpath(pwd))

%% load dataset
data = load('data/data_HPPC.mat');

% automatically detect variables
vars = fieldnames(data);

t_data = data.(vars{1});
I_data = data.(vars{2});
V_exp  = data.(vars{3});

%% simulation parameters
dt = t_data(2) - t_data(1);
N  = length(t_data);

solver_opts = modelParameters;
solver_opts.n_states = 4;   % reduced-order SPM

%% run SPM simulation
[Vsim, soc_n, soc_p] = SPM_sim([], t_data, I_data, 25, solver_opts);

SOC_true = soc_calculation(soc_n, soc_p, solver_opts);

%% EKF initialization

x = SOC_true(1);      % initial SOC
P = 1e-4;
Q = 1e-6;
R = 1e-4;

SOC_est = zeros(N,1);

%% EKF loop

for k = 2:N
    
    I = I_data(k);
    
    %% prediction step
    x_pred = x - dt * I / solver_opts.Q_nom;
    P_pred = P + Q;
    
    %% voltage model
    theta_n = solver_opts.theta_n_0 + ...
        x_pred*(solver_opts.theta_n_100 - solver_opts.theta_n_0);
    
    theta_p = solver_opts.theta_p_0 + ...
        x_pred*(solver_opts.theta_p_100 - solver_opts.theta_p_0);
    
    V_pred = U_p(theta_p) - U_n(theta_n) - solver_opts.R0 * I;
    
    %% measurement model linearization
    dVdSOC = (U_p(theta_p) - U_n(theta_n)) / max(x_pred,1e-6);
    
    %% Kalman gain
    K = P_pred * dVdSOC / (dVdSOC * P_pred * dVdSOC + R);
    
    %% update step
    x = x_pred + K * (V_exp(k) - V_pred);
    P = (1 - K*dVdSOC) * P_pred;
    
    SOC_est(k) = x;
    
end

%% plot results

figure
plot(t_data, SOC_true, 'k','linewidth',2)
hold on
plot(t_data, SOC_est,'r--','linewidth',2)

xlabel('Time (s)')
ylabel('SOC')
title('Option 2: EKF SOC Estimation')
legend('True SOC','EKF Estimate')
grid on

%% compute error

RMSE = sqrt(mean((SOC_est - SOC_true).^2));

fprintf('EKF SOC RMSE = %.4f\n', RMSE)