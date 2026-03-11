%% ===========================================================================================
%% U_p - CATHODE OPEN CIRCUIT POTENTIAL
%% ===========================================================================================
%
% DESCRIPTION:
%   Computes the open circuit potential (OCP) of the NMC (Nickel Manganese Cobalt oxide)
%   cathode as a function of lithium stoichiometry.
% INPUTS:
%   x       - Stoichiometry (dimensionless lithium concentration) [-]
%             Range: 0 (fully lithiated) to 1 (delithiated)
%   T_core  - Core temperature [K] (not used in current implementation)
%   param   - Parameter structure (not used in current implementation)
%
% OUTPUTS:
%   y       - Cathode open circuit potential [V] vs. Li/Li+
%
% NOTES:
%   - Fitted from experimental data for NMC

function y = U_p(x, T_core, param)
    % Multi-term empirical fit for NMC OCP
    % Form: linear + multiple tanh terms
    
    y = -0.8090*x + 4.4875 - ...
        0.0428*tanh(18.5138*(x-0.5542)) - ...
        17.7326*tanh(15.7890*(x-0.3117)) + ...
        17.5842*tanh(15.9308*(x-0.3120));
end
