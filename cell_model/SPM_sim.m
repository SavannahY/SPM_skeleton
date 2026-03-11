%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN SIMULATION WRAPPER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol, t_eval, param_out] = ...
            SPM_sim(param, t_data, I_data, T_amb, solver_opts)

    %% SETUP INPUT PARAMETERS
    param.t_data = t_data;
    param.I_data = I_data;

    %% TEMPERATURE
    param.T_amb = 273.15 + T_amb;
    param.V_thermal = param.Rg .* param.T_amb ./ param.F;

    %% INITIALIZE CONCENTRATION PROFILES

    if solver_opts.uniform_initial_cond

        % Get initial stoichiometry values
        [theta_s_initial, cs_n0_const, cs_p0_const] = conc_initial_sd(param, solver_opts);

        % Special initialization for Padé models
        if strcmpi(solver_opts.method,'PADE2') || strcmpi(solver_opts.method,'PADE3')

            switch upper(solver_opts.method)
                case 'PADE2'
                    order = 2;
                case 'PADE3'
                    order = 3;
                otherwise
                    error('Unknown Padé method')
            end

            % conc_initial_sd already returns concentrations in mol/m^3
            csn_avg0 = cs_n0_const;
            csp_avg0 = cs_p0_const;

            % Padé state structure:
            % [c_avg_n ; dev_states ; c_avg_p ; dev_states]
            theta_s_initial = [
                csn_avg0;
                zeros(order,1);
                csp_avg0;
                zeros(order,1)
            ];

        end

    else

        % Non-uniform initialization (kept unchanged)

        [~, cs_n0_const, cs_p0_const] = conc_initial_sd(param, solver_opts);

        [theta_s_n0, theta_s_p0] = volume_average(...
            (0.5.*(param.r_n(1:end-1) + param.r_n(2:end))).^2, ...
            (0.5.*(param.r_p(1:end-1) + param.r_p(2:end))).^2, param);

        theta_s_initial = [
            ((cs_n0_const-0.75)/(theta_s_n0-1)).* ...
                (0.5.*(param.r_n(1:end-1) + param.r_n(2:end))).^2 + ...
                (0.75 - ((cs_n0_const-0.75)/(theta_s_n0-1)));

            ((cs_p0_const-0.75)/(theta_s_p0-1)).* ...
                (0.5.*(param.r_p(1:end-1) + param.r_p(2:end))).^2 + ...
                (0.75 - ((cs_p0_const-0.75)/(theta_s_p0-1)))
        ];

    end
                 
    %% SELECT AND RUN ODE INTEGRATOR

    switch solver_opts.integrator
        
        case 'casadi'

            [sol, t_eval, param_out] = ...
                run_sim_CasADi(theta_s_initial(:), param, solver_opts);
                
        otherwise

            [sol, t_eval, param_out] = ...
                run_sim_MATLAB(theta_s_initial(:), param);

    end

    %% CHECK SOLUTION LENGTH

    if (length(t_eval) ~= length(param_out.t_data))
        disp("Something went wrong! Evaluation times not as requested!")
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%