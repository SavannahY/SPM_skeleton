%% ===========================================================================================
%% eta_cathode - CATHODE OVERPOTENTIAL (BUTLER-VOLMER)
%% ===========================================================================================
%
% DESCRIPTION:
%   Computes the activation overpotential at the cathode using the Butler-Volmer equation.
%   Similar to anode but with opposite sign due to opposite reaction direction.
%
% PHYSICS:
%   Same Butler-Volmer form as anode, but:
%   - Reaction proceeds in opposite direction
%   - Different kinetic parameters (THETA(5), THETA(6))
%   - Different geometric parameters (Rsp, Lp, csp_max)
%
% INPUTS:
%   x          - Surface stoichiometry [-] (0 to 1)
%   T_core     - Core temperature [K]
%   input_crt  - Applied current [A] (positive = discharge)
%   param      - Parameter structure with cathode parameters
%
% OUTPUTS:
%   eta_p      - Cathode overpotential [V]
%                Positive during discharge (energy gain)
%                Negative during charge
%
function eta_p = eta_cathode(x, T_core, input_crt, param)
    % Thermal voltage: RT/F [V]
    V_thermal = param.Rg.*T_core./param.F;
    
    % Butler-Volmer overpotential (note negative sign)
    eta_p = 2.0.*V_thermal.*asinh( ...
        -(param.Rsp./(param.k0p.*param.epsp.*param.Lp.*param.Acell.*param.F.*param.csp_max.*sqrt(param.ce))./6.0).*(input_crt./sqrt(x.*(1-x))) ...
    );
end
