function inventory = compute_lithium_inventory(all_data, method)
%COMPUTE_LITHIUM_INVENTORY Approximate total lithium inventory in the solid phase.
%
%   This is useful for Option 1 because FVM should preserve this quantity more
%   accurately than FDM at coarse spatial resolution.

    param = all_data.param;
    [csn_avg, csp_avg] = volume_average(all_data.csn, all_data.csp, param, method);

    % Inventory is proportional to electrode active-solid volume times the
    % average solid concentration. The proportionality constant is identical
    % across methods, so this is a fair conservation metric.
    inventory = param.Acell .* ( ...
        param.Ln .* param.epsn .* csn_avg + ...
        param.Lp .* param.epsp .* csp_avg );

    inventory = inventory(:);
end
