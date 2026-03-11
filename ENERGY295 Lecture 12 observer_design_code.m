%% ENERGY 295 - Closed-Loop Observer Design
% Design a closed-loop observer for the system:
%   d^2x/dt^2 + dx/dt + x - u = 0
% with x(0) = 0, dx/dt(0) = 0, u(t) = step of magnitude 2
% y = x (measured via sensor)
% Incorrect initial observer state: xhat(0) = [-2; -1]

clear; clc; close all;
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultAxesFontSize', 18);

%% System Matrices
% State variables: x1 = x, x2 = dx/dt
% x1_dot = x2
% x2_dot = -x1 - x2 + u
A = [0  1;
    -1 -1];
B = [0; 1];
C = [1 0];
D = 0;

%% Check Observability
V = [C; C*A];  % Observability matrix
fprintf('Observability matrix V:\n');
disp(V);
fprintf('rank(V) = %d\n\n', rank(V));

%% Observer Gain Design
% Place observer poles at s1 = s2 = -3
% Assign L from analytical solution

L = [5; 3];
fprintf('Observer gain L:\n');
disp(L);

%% Augmented State-Space System
% x_aug = [x; xhat] = [x1; x2; xhat1; xhat2]
% x_aug_dot = A_aug * x_aug + B_aug * u
% y_aug = C_aug * x_aug + D_aug * u
%
% A_aug = [A,      0;
%          LC,   A-LC]
%
% B_aug = [B; B]
%
% C_aug = [C,    0;
%          0,    C]
%
% D_aug = [0; 0]

A_aug = [A,        zeros(2);
         L*C,      A - L*C];
B_aug = [B; B];
C_aug = [C,        zeros(1,2);
         zeros(1,2), C];
D_aug = [D; D];

%% Simulation - Three cases of incorrect initial observer states
% True initial conditions (same for all cases)
x_0 = [0; 0];

% Case 1: xhat(0) = [-2; -1]
% Case 2: xhat(0) = [-6; -10]
% Case 3: xhat(0) = [-15; -20]
xhat_cases = {[-2; -1], [-6; -10], [-15; -20]};

% Time vector
t = linspace(0, 10, 5000)';

% Input: step of magnitude 2
u = 2 * ones(size(t));

% Create augmented state-space system
aug_syst = ss(A_aug, B_aug, C_aug, D_aug);

% Simulate each case
x1_true  = [];  x2_true  = [];
x1_hat   = {};  x2_hat   = {};

for k = 1:3
    xAug_0 = [x_0; xhat_cases{k}];
    [~, t_out, x_aug] = lsim(aug_syst, u, t, xAug_0);
    
    x1_true  = x_aug(:, 1);       % True x1 (same for all cases)
    x2_true  = x_aug(:, 2);       % True x2 (same for all cases)
    x1_hat{k} = x_aug(:, 3);      % Estimated x1
    x2_hat{k} = x_aug(:, 4);      % Estimated x2
end

%% Plot: True states and input
figure('Position', [368,427,498,370]);
yyaxis left
plot(t_out, x1_true, 'b-', 'LineWidth', 2); hold on;
plot(t_out, x2_true, 'b--', 'LineWidth', 1.5);
hold off;
ylabel('States');
yyaxis right
plot(t_out, u, 'r-', 'LineWidth', 1.5);
ylabel('Input $u(t)$');
xlabel('Time [s]');
title('System States and Input');
legend('$x_1$ (true)', '$x_2$ (true)', '$u(t)$', 'Location', 'best');
grid on;
xlim([0 10]);

%% Plot: Estimation of x1 and Error dynamics e1
figure('Position', [355,535,511,412]);

subplot(2,1,1);
plot(t_out, x1_true, 'b-', 'LineWidth', 2); hold on;
plot(t_out, x1_hat{1}, 'r--', 'LineWidth', 1.5);
plot(t_out, x1_hat{2}, 'g-.', 'LineWidth', 1.5);
plot(t_out, x1_hat{3}, 'k:', 'LineWidth', 1.5);
hold off;
ylabel('$x_1$');
title('Estimation of $x_1$');
legend('True value', 'Estimated (Case 1)', 'Estimated (Case 2)', ...
    'Estimated (Case 3)', 'Location', 'best');
grid on;
xlim([0 10]);

subplot(2,1,2);
for k = 1:3
    e1{k} = x1_true - x1_hat{k};
end
plot(t_out, e1{1}, 'r--', 'LineWidth', 1.5); hold on;
plot(t_out, e1{2}, 'g-.', 'LineWidth', 1.5);
plot(t_out, e1{3}, 'k:', 'LineWidth', 1.5);
yline(0, 'b-', 'LineWidth', 1);
hold off;
xlabel('Time [s]');
ylabel('$e_1 = x_1 - \hat{x}_1$');
title('Error Dynamics: $e_1$');
legend('Case 1', 'Case 2', 'Case 3', 'Location', 'best');
grid on;
xlim([0 10]);

%% Plot: Estimation of x2 and Error dynamics e2
figure('Position', [355,535,511,412]);

subplot(2,1,1);
plot(t_out, x2_true, 'b-', 'LineWidth', 2); hold on;
plot(t_out, x2_hat{1}, 'r--', 'LineWidth', 1.5);
plot(t_out, x2_hat{2}, 'g-.', 'LineWidth', 1.5);
plot(t_out, x2_hat{3}, 'k:', 'LineWidth', 1.5);
hold off;
ylabel('$x_2$');
title('Estimation of $x_2$');
legend('True value', 'Estimated (Case 1)', 'Estimated (Case 2)', ...
    'Estimated (Case 3)', 'Location', 'best');
grid on;
xlim([0 10]);

subplot(2,1,2);
for k = 1:3
    e2{k} = x2_true - x2_hat{k};
end
plot(t_out, e2{1}, 'r--', 'LineWidth', 1.5); hold on;
plot(t_out, e2{2}, 'g-.', 'LineWidth', 1.5);
plot(t_out, e2{3}, 'k:', 'LineWidth', 1.5);
yline(0, 'b-', 'LineWidth', 1);
hold off;
xlabel('Time [s]');
ylabel('$e_2 = x_2 - \hat{x}_2$');
title('Error Dynamics: $e_2$');
legend('Case 1', 'Case 2', 'Case 3', 'Location', 'best');
grid on;
xlim([0 10]);
