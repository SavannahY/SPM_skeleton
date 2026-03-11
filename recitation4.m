%% Recitation 4: ECM as a Nonlinear System, Linearization, and State-Space
%  Energy 295 — Spring 2026
%  -----------------------------------------------------------------------
%  Objectives:
%    1. Understand the first-order Equivalent Circuit Model (ECM)
%    2. See why the ECM is a nonlinear system
%    3. Formulate the discrete-time state equations
%    4. Use CasADi to define the model symbolically
%    5. Linearize the output equation using finite differences and CasADi AD
%  -----------------------------------------------------------------------
clear; clc; close all;
addpath('../../Midterm Project/SPM/casadi-3.7.2/')
import casadi.*
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultAxesFontSize', 14);

%% 0. Interpolating using MATLAB and CasADi
% =========================================================================
% For our simple example output function of y = x^2, we can build an interpolant in both
% MATLAB and CasADi to estimate the derivative at a point. The MATLAB derivate calculation
% uses a finite difference approximation, while CasADi's automatic differentiation provides
% "exact" derivatives based on the symbolic expression graph.

% y = x^2
dx = 0.01;
x = -2:dx:2;
y = x.^2;
x0 = 1;  % Point of interest

% MATLAB jacobian via finite difference approximation
matlab_interp = @(xq) interp1(x, y, xq, 'pchip');  % MATLAB interpolant
dydx_matlab = (matlab_interp(x0 + dx) - matlab_interp(x0 - dx)) / (2*dx);  % Finite difference approximation

% CasADi jacobian via symbolic expression
x = SX.sym('x');
y = Function('y', {x}, {x^2});
dydx_casadi = jacobian(y(x), x);
dydx_casadi = gradient(y(x), x);  % Alternatively, gradient can be used for scalar outputs
dydx_casadi_func = Function('dydx_casadi_func', {x}, {dydx_casadi});
dydx_casadi_value = full(dydx_casadi_func(x0));

fprintf('Derivative at x=%.2f:\n', x0);
fprintf('  MATLAB finite difference: %.4f\n', dydx_matlab);
fprintf('  CasADi automatic differentiation: %.4f\n', dydx_casadi_value);

%% 0b. Vector-Valued Jacobian Example
% =========================================================================
% Now extend the same idea to a vector-valued function f: R^2 -> R^2.
%
%   f(x1, x2) = [ x1^2 + x2   ]
%               [ sin(x1)*x2  ]
%
% The Jacobian is the 2x2 matrix of partial derivatives:
%
%   J = [ df1/dx1  df1/dx2 ]  =  [ 2*x1       1        ]
%       [ df2/dx1  df2/dx2 ]     [ cos(x1)*x2  sin(x1) ]

% Evaluation point
x0_vec = [1; 2];

% MATLAB Jacobian via finite differences
f_matlab = @(xv) [xv(1)^2 + xv(2); sin(xv(1))*xv(2)];  % Function handle
n = length(x0_vec);        % Number of inputs
m = length(f_matlab(x0_vec));  % Number of outputs
J_fd = zeros(m, n);
for j = 1:n
    e_j = zeros(n, 1);
    e_j(j) = dx;  % Perturb the j-th input by dx (defined in Section 0)
    J_fd(:, j) = (f_matlab(x0_vec + e_j) - f_matlab(x0_vec - e_j)) / (2*dx);
end

% MATLAB Jacobian via griddedInterpolant (multi-dim analog of interp1)
% Build a 2-D grid of input points and evaluate f at every grid node.
% The grid must extend beyond the evaluation point so that finite-difference
% perturbations (x0 +/- dx) stay within bounds (interpn returns NaN otherwise).
x1_grid = -3:dx:3;
x2_grid = -3:dx:3;
[X1, X2] = meshgrid(x1_grid, x2_grid);
F1_grid = X1.^2 + X2;             % f1 evaluated on the grid
F2_grid = sin(X1) .* X2;          % f2 evaluated on the grid

% Create a griddedInterpolant for each output component
F1_interp = griddedInterpolant({x1_grid, x2_grid}, F1_grid', 'cubic');
F2_interp = griddedInterpolant({x1_grid, x2_grid}, F2_grid', 'cubic');

% Finite differences on the interpolants
J_gridded = zeros(m, n);
for j = 1:n
    e_j = zeros(1, n);
    e_j(j) = dx;
    xp = x0_vec' + e_j;  xm = x0_vec' - e_j;
    J_gridded(1, j) = (F1_interp(xp) - F1_interp(xm)) / (2*dx);
    J_gridded(2, j) = (F2_interp(xp) - F2_interp(xm)) / (2*dx);
end

% MATLAB Jacobian via interpn (no object creation needed)
J_interpn = zeros(m, n);
for j = 1:n
    e_j = zeros(1, n);
    e_j(j) = dx;
    xp = x0_vec' + e_j;  xm = x0_vec' - e_j;
    J_interpn(1, j) = (interpn(x1_grid, x2_grid, F1_grid', xp(1), xp(2), 'cubic') ...
                      - interpn(x1_grid, x2_grid, F1_grid', xm(1), xm(2), 'cubic')) / (2*dx);
    J_interpn(2, j) = (interpn(x1_grid, x2_grid, F2_grid', xp(1), xp(2), 'cubic') ...
                      - interpn(x1_grid, x2_grid, F2_grid', xm(1), xm(2), 'cubic')) / (2*dx);
end

% CasADi Jacobian via automatic differentiation
x_vec = SX.sym('x_vec', 2, 1);

% Define the vector-valued function symbolically
f = [x_vec(1)^2 + x_vec(2);      ...
     sin(x_vec(1)) * x_vec(2)];

% Wrap into a CasADi Function
f_func = Function('f', {x_vec}, {f});

% Compute the full Jacobian matrix symbolically
J = jacobian(f, x_vec);
J_func = Function('J_func', {x_vec}, {J});
J_value = full(J_func(x0_vec));

% Hand-computed Jacobian for verification
J_expected = [2*x0_vec(1),              1;            ...
              cos(x0_vec(1))*x0_vec(2), sin(x0_vec(1))];

fprintf('Vector-valued Jacobian example:\n');
fprintf('  f(x) = [x1^2 + x2; sin(x1)*x2]\n');
fprintf('  Evaluated at x = [%.1f; %.1f]:\n', x0_vec(1), x0_vec(2));
fprintf('  Jacobian (MATLAB finite difference on function handle):\n');
disp(J_fd);
fprintf('  Jacobian (MATLAB griddedInterpolant + finite difference):\n');
disp(J_gridded);
fprintf('  Jacobian (MATLAB interpn + finite difference):\n');
disp(J_interpn);
fprintf('  Jacobian (CasADi):\n');
disp(J_value);
fprintf('  Jacobian (hand-computed):\n');
disp(J_expected);

%% 0c. Linearization of a Nonlinear ODE System
% =========================================================================
% Consider the third-order nonlinear ODE:
%
%   d^3x/dt^3 + 8 d^2x/dt^2 + 30 dx/dt + 3 x^4 = (3/5) u
%
% State-space form with x1 = x,  x2 = dx/dt,  x3 = d^2x/dt^2:
%
%   x1_dot = x2
%   x2_dot = x3
%   x3_dot = -3*x1^4 - 30*x2 - 8*x3 + (3/5)*u
%
% The system is nonlinear because of the x1^4 term.
%
% Steps:
%   1. Find the equilibrium point for u_e = 25
%   2. Linearize around the equilibrium using MATLAB finite differences
%   3. Linearize around the equilibrium using CasADi automatic differentiation
%   4. Compare with the hand-computed Jacobian

% --- Define the nonlinear state-space function ---
f = @(xv, uv) [xv(2);
                    xv(3);
                    -3*xv(1)^4 - 30*xv(2) - 8*xv(3) + (3/5)*uv];

% --- Find equilibrium point for u_e = 25 ---
% At equilibrium, x_dot = 0:
%   x2e = 0,  x3e = 0,  -3*x1e^4 + (3/5)*u_e = 0
%   => x1e^4 = (3/5)*u_e / 3 = u_e/5
%   => x1e = +/- (u_e/5)^(1/4)
u_e = 25;
x1e = (u_e / 5)^(1/4);   % Positive equilibrium

xe_pos = [x1e; 0; 0];

% --- Method 1: MATLAB finite-difference linearization ---
% The linearized model is:  delta_x_dot = F*delta_x + H*delta_u
% where F = df/dx  and  H = df/du  evaluated at (x_e, u_e).
dx_hw = 1e-7;
n_hw = 3;  % Number of states

% Compute F = df/dx via central finite differences
F_fd = zeros(n_hw, n_hw);
for j = 1:n_hw
    ej = zeros(n_hw, 1);
    ej(j) = dx_hw;
    F_fd(:, j) = (f(xe_pos + ej, u_e) - f(xe_pos - ej, u_e)) / (2*dx_hw);
end

% Compute H = df/du via central finite differences
H_fd = (f(xe_pos, u_e + dx_hw) - f(xe_pos, u_e - dx_hw)) / (2*dx_hw);

fprintf('MATLAB Finite-Difference Linearization\n');
fprintf('At x_e = [%.4f; 0; 0]:\n', x1e);
fprintf('  F = df/dx =\n'); disp(F_fd);
fprintf('  H = df/du =\n'); disp(H_fd);

% --- Method 2: CasADi automatic differentiation ---
x_hw1 = SX.sym('x_hw1', 3, 1);
u_hw1 = SX.sym('u_hw1');

% Define the state equations symbolically
f_hw1_sym = [x_hw1(2);
             x_hw1(3);
             -3*x_hw1(1)^4 - 30*x_hw1(2) - 8*x_hw1(3) + (3/5)*u_hw1];

% Compute the Jacobians symbolically
F_hw1_sym = jacobian(f_hw1_sym, x_hw1);   % 3x3: df/dx
H_hw1_sym = jacobian(f_hw1_sym, u_hw1);   % 3x1: df/du

% Wrap as callable CasADi Functions and evaluate at equilibrium
F_hw1_func = Function('F_hw1', {x_hw1, u_hw1}, {F_hw1_sym});
H_hw1_func = Function('H_hw1', {x_hw1, u_hw1}, {H_hw1_sym});

F_casadi_pos = full(F_hw1_func(xe_pos, u_e));
H_casadi_pos = full(H_hw1_func(xe_pos, u_e));

fprintf('CasADi\n');
fprintf('At x_e = [%.4f; 0; 0]:\n', x1e);
fprintf('  F =\n'); disp(F_casadi_pos);
fprintf('  H =\n'); disp(H_casadi_pos);

% --- Hand-computed Jacobian for verification ---
%   F = [  0          1     0  ]
%       [  0          0     1  ]
%       [ -12*x1e^3  -30   -8  ]
%
%   H = [0; 0; 3/5]
F_exact = [0, 1, 0; 0, 0, 1; -12*x1e^3, -30, -8];
H_exact = [0; 0; 3/5];

fprintf('Hand-Computed Jacobians\n');
fprintf('At x_e = [%.4f; 0; 0]:\n', x1e);
fprintf('  F =\n'); disp(F_exact);
fprintf('  H = [0; 0; 3/5] for both:\n'); disp(H_exact);

%% 1. Cell Parameters
% =========================================================================
% We consider a single lithium-ion cell described by a first-order ECM
%
%   Terminal voltage:  V_t = OCV(z) - V1 - R0*I
%   RC dynamics:       dV1/dt = -(1/(R1*C1))*V1 + (1/C1)*I
%   SoC dynamics:      dz/dt  = -I / Q
%
% where:
%   z     = state of charge (SoC), dimensionless [0,1]
%   I     = applied current (positive = discharge)  [A]
%   V1    = voltage across the RC pair              [V]
%   Q     = cell capacity                            [As]
%   OCV(z)= open-circuit voltage, a nonlinear function of SoC

Q_Ah = 4.8687;                 % Nominal capacity [Ah] (from C/20 test)
Q    = Q_Ah * 3600;           % Capacity in [As] (Coulombs)
R0   = 0.01;                  % Series (ohmic) resistance [Ohm]
R1   = 0.015;                 % Polarization resistance   [Ohm]
C1   = 2500;                  % Polarization capacitance  [F]
tau1 = R1 * C1;               % RC time constant [s]

%% 2. Open-Circuit Voltage from Experimental C/20 Discharge
% =========================================================================
% OCV(z) is the key nonlinearity in the ECM.  A C/20 discharge is slow
% enough that the terminal voltage ≈ OCV at each SoC.  We load the
% experimental C/20 data, compute SoC via Coulomb counting, and build
% a B-spline interpolant to obtain a smooth, differentiable OCV(SoC).
%
% Why interpolation instead of a polynomial fit?
%   - CasADi can differentiate through its interpolant automatically

% Load C/20 experimental data
data_path = 'data_co20.mat';
raw       = load(data_path);
I_co20    = raw.output.I;     % Current [A] (negative = discharge)
V_co20    = raw.output.V;     % Terminal voltage [V]
t_co20    = raw.output.t;     % Time [s]

% Compute SoC from Coulomb counting
% Current is negative during discharge; flip sign so I_dchg > 0
I_dchg   = -I_co20;
Q_total  = trapz(t_co20, I_dchg);           % Total charge [As]
soc_co20 = 1 - cumtrapz(t_co20, I_dchg) / Q_total;  % SoC: 1 -> 0

fprintf('C/20 data loaded: %d samples, Q = %.2f Ah\n', ...
        length(t_co20), Q_total/3600);

% Build a clean, uniformly-spaced SoC-vs-OCV lookup table
[soc_sorted, sort_idx] = sort(soc_co20);
V_sorted = V_co20(sort_idx);

% Remove duplicate SoC values (e.g., from rest periods)
[soc_unique, unique_idx] = unique(soc_sorted, 'stable');
V_unique = V_sorted(unique_idx);

% Resample onto a uniform grid using MATLAB pchip interpolation
N_grid   = 500;
soc_grid = linspace(soc_unique(1), soc_unique(end), N_grid);
V_grid   = interp1(soc_unique, V_unique, soc_grid, 'pchip');

fprintf('OCV grid: %d points, SoC range [%.4f, %.4f]\n', ...
        N_grid, soc_grid(1), soc_grid(end));

% --- CasADi B-spline interpolant ---
% CasADi's interpolant() builds a callable Function from gridded data.
% 'bspline' gives a smooth curve whose derivatives are well-defined
% everywhere — essential for computing Jacobians via automatic diff.
OCV_interp = casadi.interpolant('OCV', 'bspline', {soc_grid}, V_grid);

% --- CasADi symbolic OCV and its derivative ---
SOC_sym    = SX.sym('SoC');
OCV_sym  = OCV_interp(SOC_sym);              % OCV(SoC) as SX expression
dOCV_sym = jacobian(OCV_sym, SOC_sym);       % dOCV/dSOC via automatic diff

% Wrap as callable CasADi Functions
OCV_cas  = Function('OCV',  {SOC_sym}, {OCV_sym});
dOCV_cas = Function('dOCV', {SOC_sym}, {dOCV_sym});

% --- Use CasADi functions for plotting and comparison ---
% We wrap them in full() to convert CasADi DM results to MATLAB doubles
OCV  = @(SoC) full(OCV_cas(SoC));
dOCV = @(SoC) full(dOCV_cas(SoC));

fprintf('OCV B-spline interpolant created\n\n');

% --- Plot: experimental data vs interpolant ---
SOC_interpvec = linspace(soc_grid(1), soc_grid(end), 500);
figure('Name','OCV Curve','Position',[235,602,775,248]);

subplot(1,2,1);
plot(soc_co20, V_co20, 'k.', 'MarkerSize', 1); hold on;
plot(SOC_interpvec, OCV(SOC_interpvec), 'b-', 'LineWidth', 2);
xlabel('SoC [-]'); ylabel('OCV [V]');
title('Open-Circuit Voltage vs. SoC');
legend('C/20 data', 'Interpolant', 'Location', 'best');
grid on; xlim([0 1]);

subplot(1,2,2);
plot(SOC_interpvec, dOCV(SOC_interpvec), 'r-', 'LineWidth', 2);
xlabel('SoC [-]'); ylabel('dOCV/dSOC [V]');
title('Slope of OCV vs. SoC');
grid on; xlim([0 1]);

%% 3. Nonlinear State-Space Model: Discrete- Time Formulation
% =========================================================================
% For simulation and estimation we work in discrete time with step dt.

dt = 1;  % Time step [s]

% --- (a) Forward Euler step function ---
f_euler = @(xk, uk) xk + dt * [ -uk/Q; -(1/tau1)*xk(2) + (R1/tau1)*uk ];

% --- (b) Exact discrete-time matrices ---
% We use the augmented matrix exponential to handle the singular A.
%   A_cont = [  0        0      ]    B_cont = [ -1/Q      ]
%            [  0    -1/tau1    ]             [ R1/tau1   ]
A_cont = [0, 0; 0, -1/tau1];
B_cont = [-1/Q; R1/tau1];
n_x = size(A_cont, 1);
n_u = size(B_cont, 2);
M_aug = expm([A_cont, B_cont; zeros(n_u, n_x + n_u)] * dt);
A_bar = M_aug(1:n_x, 1:n_x);
B_bar = M_aug(1:n_x, n_x+1:end);

fprintf('=== Discrete-Time Matrices (dt = %d s) ===\n', dt);
fprintf('A_bar = expm(A*dt):\n'); disp(A_bar);
fprintf('B_bar (from augmented matrix exponential):\n'); disp(B_bar);

% --- (c) CasADi integrator: built in Section 4 after symbolic model ---
% (Requires f_sym, which is defined in the next section.)

fprintf('Forward-Euler and exact discrete-time models ready.\n\n');


%% 4. Define the Nonlinear Model Using CasADi
% =========================================================================
% CasADi lets us define the model SYMBOLICALLY.  We declare symbolic
% variables for the states and input, write the equations once, and then
% CasADi can:
%   (a) evaluate them numerically (like anonymous functions)
%   (b) differentiate them automatically (Jacobians for linearization)
%   (c) build ODE integrators for simulation

% --- Symbolic state and input ---
x_sym = SX.sym('x', 2, 1);   % x = [SoC; V1]
u_sym = SX.sym('u');          % u = I  (scalar)

% Unpack states for readability
SOC_s  = x_sym(1);            % SoC
V1_s   = x_sym(2);            % RC voltage

% --- Symbolic OCV via CasADi interpolant (built in Section 2) ---
OCV_of_SOC = OCV_interp(SOC_s);

% --- State equations:  dx/dt = f(x, u) ---
f_sym = [ -u_sym / Q;
          -(1/tau1)*V1_s + (R1/tau1)*u_sym ];

% --- Output equation:  y = g(x, u) ---
g_sym = OCV_of_SOC - V1_s - R0*u_sym;

% --- Standalone voltage function: voltage_func(x, I) -> V_t ---
% Build with a FRESH symbolic variable (independent of x_sym above).
% This is important: when you later want to differentiate voltage w.r.t.
% states, CasADi needs a pure symbolic input — not one that's already
% embedded in another expression graph.
x_v = SX.sym('x_v', 2, 1);
I_v = SX.sym('I_v');
% OCV via interpolant on the fresh variable x_v(1)
V_standalone = OCV_interp(x_v(1)) - x_v(2) - R0*I_v;
voltage_func = Function('voltage_func', {x_v, I_v}, {V_standalone});


%% 5. Linearization of the Output Equation
% =========================================================================
% The output (terminal voltage) is:
%   V_t = OCV(SoC) - V1 - R0*I  =  g(x, u)
%
% This is nonlinear because OCV(SoC) is a nonlinear function of SoC.
% To linearize around an operating point (x0, u0), we compute the Jacobian
% of g with respect to x = [SoC; V1] and u = I:
%
%   C = dg/dx = [dOCV/dSoC, -1]     (1x2)
%   D = dg/du = -R0                  (scalar)
%
% so that near the operating point:
%   V_t ≈ V_t0 + C*(x - x0) + D*(u - u0)
%
% We compute C two ways: (1) MATLAB finite differences, (2) CasADi AD.

% --- Operating point ---
z0     = 0.50;            % SoC = 50%
I0     = 0;               % No current
V1_0   = 0;               % Steady-state: V1 = 0 at equilibrium
x0     = [z0; V1_0];
u0     = I0;
y0     = full(voltage_func(x0, u0));

fprintf('=== Operating Point ===\n');
fprintf('  z0     = %.2f\n', z0);
fprintf('  V1_0   = %.2f V\n', V1_0);
fprintf('  I0     = %.2f A\n', I0);
fprintf('  V_t0   = OCV(z0) - V1_0 - R0*I0 = %.4f V\n\n', y0);

% Method 1
% Evaluate the nonlinear output g(x, u) = voltage_func(x, u)
% and perturb each input by dx to approximate the partial derivatives.
g_eval = @(xv, uv) full(voltage_func(xv, uv));

C_fd = zeros(1, 2);
dx = 1e-5;
for j = 1:2
    e_j = zeros(2, 1);
    e_j(j) = dx;
    C_fd(j) = (g_eval(x0 + e_j, u0) - g_eval(x0 - e_j, u0)) / (2*dx);
end
D_fd = (g_eval(x0, u0 + dx) - g_eval(x0, u0 - dx)) / (2*dx);

fprintf('MATLAB Finite-Difference Linearization of g(x,u)\n');
fprintf('  C = dg/dx = [%.4f, %.4f]\n', C_fd(1), C_fd(2));
fprintf('  D = dg/du = %.4f\n\n', D_fd);

% CasADi computes exact symbolic derivatives of g w.r.t. x and u.
C_sym = jacobian(g_sym, x_sym);   % 1x2: dg/dx
D_sym = jacobian(g_sym, u_sym);   % 1x1: dg/du

fprintf('CasADi Symbolic Jacobians of g(x,u)\n');
fprintf('  C_sym = '); disp(C_sym);
fprintf('  D_sym = '); disp(D_sym);

% Wrap as CasADi Functions and evaluate at operating point
C_func = Function('C', {x_sym, u_sym}, {C_sym});
D_func = Function('D', {x_sym, u_sym}, {D_sym});
C_casadi = double(full(C_func(x0, u0)));
D_casadi = double(full(D_func(x0, u0)));

fprintf('CasADi Evaluated at Operating Point\n');
fprintf('  C = dg/dx = [%.4f, %.4f]\n', C_casadi(1), C_casadi(2));
fprintf('  D = dg/du = %.4f\n\n', D_casadi);

% --- Comparison ---
fprintf('Comparison\n');
fprintf('  C(1) = dOCV/dSoC at SoC=%.2f:\n', z0);
fprintf('    MATLAB FD:  %.6f\n', C_fd(1));
fprintf('    CasADi AD:  %.6f\n', C_casadi(1));
fprintf('    dOCV func:  %.6f\n', dOCV(z0));
fprintf('  C(2) = -1 (dg/dV1):\n');
fprintf('    MATLAB FD:  %.6f\n', C_fd(2));
fprintf('    CasADi AD:  %.6f\n', C_casadi(2));
fprintf('    Exact:      %.6f\n', -1);
fprintf('  D = -R0:\n');
fprintf('    MATLAB FD:  %.6f\n', D_fd);
fprintf('    CasADi AD:  %.6f\n', D_casadi);
fprintf('    Exact:      %.6f\n', -R0);

%% 6. Comparison: Nonlinear vs. Linearized Output
% =========================================================================
% We simulate a constant-current discharge and compare three ways of
% computing the terminal voltage from the SAME state trajectory:

T_sim  = 1800;              % Simulate 30 minutes [s]
N_sim  = T_sim / dt;        % Number of time steps
t_vec  = (0:N_sim) * dt;    % Time vector [s]

% Initial conditions: start at 80% SoC, no polarization
z_init   = 0.80;
V1_init  = 0;
x_sim    = [z_init; V1_init];

% Input: constant 1C discharge
I_discharge = Q_Ah;           % 1C rate in Amps (positive = discharge)
figure('Position', [232,146,560,420]);
plot(I_discharge * ones(length(t_vec)), 'r-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Current [A]');
title('Input Current: 1C Discharge');
grid on;

x_traj    = zeros(2, N_sim + 1);     % State trajectory
V_nonlin  = zeros(1, N_sim + 1);     % Nonlinear output
V_fixed   = zeros(1, N_sim + 1);     % Fixed-Jacobian linearized output
V_updated = zeros(1, N_sim + 1);     % Updated-Jacobian linearized output

x_traj(:, 1) = x_sim;

% Fixed Linearization
OCV_at_z0  = OCV(z_init);            % OCV(z0)
dOCV_at_z0 = dOCV(z_init);           % dOCV/dz at z0

% Simulation loop
SoC_prev = z_init;

for k = 1:(N_sim + 1)
    SoC_k = x_traj(1, k);
    V1_k  = x_traj(2, k);

    % 1. Nonlinear output: use full OCV lookup
    V_nonlin(k) = OCV(SoC_k) - V1_k - R0*I_discharge;

    % 2. Fixed linearization: OCV ~ OCV(z0) + dOCV/dz|_{z0} * (SoC - z0)
    V_fixed(k) = OCV_at_z0 + dOCV_at_z0*(SoC_k - z_init) ...
                 - V1_k - R0*I_discharge;

    % 3. Updated linearization: OCV ~ OCV(SoC_{k-1}) + dOCV/dz|_{SoC_{k-1}} * (SoC_k - SoC_{k-1})
    if k == 1
        V_updated(k) = V_nonlin(k);   % First step: no previous point
    else
        V_updated(k) = OCV(SoC_prev) + dOCV(SoC_prev)*(SoC_k - SoC_prev) ...
                       - V1_k - R0*I_discharge;
    end
    SoC_prev = SoC_k;

    % Propagate state to next step (exact discretization from Section 3b)
    if k <= N_sim
        x_traj(:, k+1) = A_bar * x_traj(:, k) + B_bar * I_discharge;
    end
end

% --- Plot results ---
figure('Name','Linearized vs. Nonlinear Output','Position',[685,143,563,444]);

% Terminal voltage comparison
subplot(2,2,1);
plot(t_vec/60, V_nonlin,  'k-',  'LineWidth', 2); hold on;
plot(t_vec/60, V_fixed,   'r--', 'LineWidth', 1.5);
plot(t_vec/60, V_updated, 'b-.', 'LineWidth', 1.5);
xlabel('Time [min]'); ylabel('$V_t$ [V]');
title('Terminal Voltage: Nonlinear vs.\ Linearized');
legend('Nonlinear', 'Fixed Jacobian', 'Updated Jacobian', 'Location', 'best');
grid on;

% Linearization error
subplot(2,2,2);
plot(t_vec/60, (V_fixed - V_nonlin)*1000,   'r-', 'LineWidth', 1.5); hold on;
plot(t_vec/60, (V_updated - V_nonlin)*1000,  'b-', 'LineWidth', 1.5);
xlabel('Time [min]'); ylabel('Error [mV]');
title('Linearization Error vs.\ Nonlinear');
legend('Fixed Jacobian', 'Updated Jacobian', 'Location', 'best');
grid on;

% SoC trajectory
subplot(2,2,3);
plot(t_vec/60, x_traj(1,:), 'k-', 'LineWidth', 2);
xlabel('Time [min]'); ylabel('SoC [-]');
title('State of Charge During Discharge');
grid on;

% OCV approximation comparison
SoC_range = linspace(x_traj(1,end), z_init, 200);
OCV_true  = OCV(SoC_range);
OCV_lin   = OCV_at_z0 + dOCV_at_z0 * (SoC_range - z_init);

subplot(2,2,4);
plot(SoC_range, OCV_true, 'k-',  'LineWidth', 2); hold on;
plot(SoC_range, OCV_lin,  'r--', 'LineWidth', 1.5);
plot(z_init, OCV_at_z0, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('SoC [-]'); ylabel('OCV [V]');
title('OCV: True Curve vs.\ Fixed Tangent Line');
legend('OCV(SoC)', 'Lin. at $z_0$', 'Linearization point', 'Location', 'best');
grid on;
