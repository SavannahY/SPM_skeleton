%% ===========================================================================================
%% soc_calculation - STATE OF CHARGE COMPUTATION
%% ===========================================================================================
%
% DESCRIPTION:
%   Converts volume-averaged stoichiometry to State of Charge (SOC) for both electrodes.
%   SOC is normalized to [0,1] based on the usable stoichiometry range.
%
% PHYSICS:
%   SOC is defined relative to the operating stoichiometry window:
%   
%   Anode:   SOC_n = (θ_n - θ_n,0%) / (θ_n,100% - θ_n,0%)
%   Cathode: SOC_p = (θ_p - θ_p,0%) / (θ_p,100% - θ_p,0%)
%   
%   Note the inverted formula for cathode (lithium moves out during charge)
%
% INPUTS:
%   theta_s_n  - Anode stoichiometry profile [dimensionless]
%   theta_s_p  - Cathode stoichiometry profile [dimensionless]
%   param      - Parameter structure with stoichiometry limits:
%                .THETA(7) = θ_n at 100% SOC
%                .THETA(8) = θ_p at 100% SOC
%                .THETA(9) = θ_n at 0% SOC
%                .THETA(10) = θ_p at 0% SOC
%
% OUTPUTS:
%   soc_bulk_n    - Anode bulk SOC [-] (0 = discharged, 1 = charged)
%   soc_bulk_p    - Cathode bulk SOC [-] (0 = discharged, 1 = charged)
%   theta_bulk_n  - Anode volume-averaged stoichiometry [-]
%   theta_bulk_p  - Cathode volume-averaged stoichiometry [-]
%
% NOTES:
%   - SOC = 1 means fully charged (anode full, cathode empty)
%   - SOC = 0 means fully discharged (anode empty, cathode full)
%   - Both electrodes should have similar SOC (slight mismatch possible)
function [soc_bulk_n, soc_bulk_p, theta_bulk_n, theta_bulk_p] = ...
    soc_calculation(theta_s_n, theta_s_p, param, method)

    % Get number of time points
    num_t = size(theta_s_n, 2);
    
    % Initialize output arrays
    theta_bulk_n = zeros(num_t, 1);
    theta_bulk_p = zeros(num_t, 1);

    % Compute volume average for each time point
    for k = 1:num_t
        [theta_bulk_n(k), theta_bulk_p(k)] = volume_average(...
            theta_s_n(:,k), theta_s_p(:,k), param, method);
    end

    % Ensure column vectors
    theta_bulk_n = theta_bulk_n(:);
    theta_bulk_p = theta_bulk_p(:);

    % Convert stoichiometry to SOC using operating window
    % Anode: lithium content increases with SOC
    soc_bulk_n = (theta_bulk_n - param.theta_n_0) / ...
                 (param.theta_n_100 - param.theta_n_0);
    
    % Cathode: lithium content decreases with SOC (inverted)
    soc_bulk_p = (theta_bulk_p - param.theta_p_0) / ...
                 (param.theta_p_100 - param.theta_p_0);
end
