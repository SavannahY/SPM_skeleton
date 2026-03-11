function n_states = utils_n_states(param)
%UTILS_N_STATES Compute number of states per electrode based on discretization.
%
%   n_states = UTILS_N_STATES(param) returns the number of state variables
%   per electrode based on the discretization method specified in param.method
%   and the number of radial nodes param.Nr.
%
%   Supported methods:
%       'FDM'   - Finite Difference Method
%       'FVM'   - Finite Volume Method
%       'PADE2' - 2nd order Nonlinear SPM-Padé model
%       'PADE3' - 3rd order Nonlinear SPM-Padé model
%       'LLPADEECM' - Local Linear Padé-ECM model
%
%   Inputs:
%       param.method
%       param.Nr
%
%   Output:
%       n_states - number of states per electrode

    switch upper(param.method)

        case 'FVM'
            % FVM uses cell-centered control volumes
            % results in Nr-1 states
            n_states = param.Nr - 1;

        case 'FDM'
            % FDM uses Nr nodes
            n_states = param.Nr;

        case 'PADE2'
            % Padé order 2:
            %   1 average concentration
            %   2 diffusion deviation states
            n_states = 3;

        case 'PADE3'
            % Padé order 3:
            %   1 average concentration
            %   3 diffusion deviation states
            n_states = 4;

        case 'LLPADEECM'
            % Local Linear Padé-ECM:
            %   3 anode deviation states + average concentration
            % per electrode-equivalent block
            n_states = 4;

        otherwise
            error('Unknown method: %s', param.method)

    end

end