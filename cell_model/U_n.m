%% ===========================================================================================
%% U_n - ANODE OPEN CIRCUIT POTENTIAL
%% ===========================================================================================
%
% DESCRIPTION:
%   Computes the open circuit potential (OCP) of the graphite anode as a function
%   of lithium stoichiometry (normalized concentration).
%
% INPUTS:
%   x       - Stoichiometry (dimensionless lithium concentration) [-]
%             Range: 0 (no lithium) to 1 (fully lithiated)
%   T_core  - Core temperature [K] (not used in current implementation)
%   param   - Parameter structure (not used in current implementation)
%
% OUTPUTS:
%   y       - Anode open circuit potential [V] vs. Li/Li+
function y = U_n(x, T_core, param)
    % Multi-term empirical fit for graphite OCP
    % Form: exponential + constant + multiple tanh terms

    y = 1.9793*exp(-39.3631*x) + ...
        0.2482 - ...
        0.0909*tanh(29.8538*(x-0.1234)) - ...
        0.04478*tanh(14.9159*(x-0.2769)) - ...
        0.0205*tanh(30.4444*(x-0.6103));
end
