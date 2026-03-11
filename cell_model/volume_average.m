function [theta_s_n_ave, theta_s_p_ave] = volume_average(theta_s_n, theta_s_p, param, method)
%VOLUME_AVERAGE Compute volume-averaged concentration / stoichiometry
%
% Inputs:
%   theta_s_n, theta_s_p : radial profiles OR average states
%   param                : parameter structure
%   method               : 'FDM', 'FVM', 'PADE2', 'PADE3'
%
% Outputs:
%   theta_s_n_ave, theta_s_p_ave

    switch upper(method)

        case 'FDM'
            % radial nodes
            r_n = param.r_n(:);
            r_p = param.r_p(:);

            % trapezoidal integration in spherical coordinates
            num_t = size(theta_s_n,2);
            theta_s_n_ave = zeros(num_t,1);
            theta_s_p_ave = zeros(num_t,1);

            for k = 1:num_t
                theta_s_n_ave(k) = 3 * trapz(r_n, theta_s_n(:,k) .* (r_n.^2));
                theta_s_p_ave(k) = 3 * trapz(r_p, theta_s_p(:,k) .* (r_p.^2));
            end

        case 'FVM'
            % shell-centered control volumes
            r_n_edges = param.r_n(:);
            r_p_edges = param.r_p(:);

            r_n_mid = 0.5 * (r_n_edges(1:end-1) + r_n_edges(2:end));
            r_p_mid = 0.5 * (r_p_edges(1:end-1) + r_p_edges(2:end));

            dr_n = diff(r_n_edges);
            dr_p = diff(r_p_edges);

            w_n = 3 * (r_n_mid.^2) .* dr_n;
            w_p = 3 * (r_p_mid.^2) .* dr_p;

            theta_s_n_ave = (w_n.' * theta_s_n).';
            theta_s_p_ave = (w_p.' * theta_s_p).';

        case {'PADE2','PADE3','LLPADEECM'}
            % For Padé, theta_s_n / theta_s_p are already average quantities
            theta_s_n_ave = theta_s_n(:);
            theta_s_p_ave = theta_s_p(:);

        otherwise
            error('Unknown method: %s', method)
    end

    theta_s_n_ave = theta_s_n_ave(:);
    theta_s_p_ave = theta_s_p_ave(:);
end