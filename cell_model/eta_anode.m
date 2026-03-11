%% ===========================================================================================
%% eta_anode - ANODE OVERPOTENTIAL (BUTLER-VOLMER)
%% ===========================================================================================
%
% DESCRIPTION:
%   Computes the activation overpotential at the anode using the Butler-Volmer equation.
%   Overpotential is the voltage loss due to slow kinetics of the lithium intercalation
%   reaction at the electrode-electrolyte interface.
%
% PHYSICS:
%   Butler-Volmer equation (linearized form):
%   η = (2RT/F) * asinh(i/(2*i₀*√(θ*(1-θ))))
%   
%   Where:
%   - i₀: Exchange current density (function of concentration)
%   - θ: Surface stoichiometry
%   - (1-θ): Available sites for lithium
%
% INPUTS:
%   x          - Surface stoichiometry [-] (0 to 1)
%   T_core     - Core temperature [K]
%   input_crt  - Applied current [A] (positive = discharge)
%   param      - Parameter structure with:
%                .Rsn, .Ln, .csn_max, .F, .Acell, .ce
%                .THETA(2), .THETA(3) - dimensionless kinetic parameters
%
% OUTPUTS:
%   eta_n      - Anode overpotential [V]
%                Positive during discharge (energy loss)
%                Negative during charge
%
% NOTES:
%   - Overpotential increases with current (faster = more loss)
%   - Overpotential increases near concentration limits (θ→0 or θ→1)
%   - The √(θ*(1-θ)) term comes from concentration dependence of exchange current
%
function eta_n = eta_anode(x, T_core, input_crt, param)
    % Thermal voltage: RT/F [V]
    % At 25°C: V_thermal ≈ 0.0257 V
    V_thermal = param.Rg.*T_core./param.F;
    
    % Butler-Volmer overpotential calculation
    eta_n = 2.0.*V_thermal.*asinh( ...
        (param.Rsn./(param.k0n.*param.epsn.*param.Ln.*param.Acell.*param.F.*param.csn_max.*sqrt(param.ce))./6.0).*(input_crt./sqrt(x.*(1-x))) ...
    );
end
