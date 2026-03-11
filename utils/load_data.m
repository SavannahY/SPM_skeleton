function [t_data, I_data, V_data, solver_opts] = load_data(filename, dt, Tend)
    %LOAD_DATA Load experimental profile data from MAT-file
    %
    % INPUTS:
    %   filename - Path to MAT-file containing experimental data
    %   dt       - Desired time step for resampling [s]
    %   Tend     - End time for simulation [s]
    %
    % OUTPUTS:
    %   t_data      - Time vector [s]
    %   I_data      - Current vector [A]
    %   V_data      - Voltage vector [V]
    %   solver_opts - Structure with interpolant functions
    
    load(filename);
    % Create uniform time grid for interpolation
    t_data = linspace(0.0, output.t(end), ceil(output.t(end)./dt) + 1)';  % [s]  
    t_data = t_data(:);  % Ensure column vector
    t_data = t_data(t_data <= Tend);  % Truncate if Tend < profile length

    %% CREATE INTERPOLATION FUNCTIONS
    % Create interpolant for current profile
    % "previous" = zero-order hold (constant between samples)
    % "nearest" = extrapolation method for out-of-bounds queries
    solver_opts.I_interpolant = griddedInterpolant(output.t(:), -output.I(:), ...
                                                "previous", "nearest");

    % Create interpolant for voltage (for comparison/validation)
    solver_opts.V_interpolant = griddedInterpolant(output.t(:), output.V(:), ...
                                    "previous", "nearest");

    % Evaluate interpolants at simulation time points
    I_data = solver_opts.I_interpolant(t_data);  % Applied current [A]
    V_data = solver_opts.V_interpolant(t_data);            % Reference voltage [V]
end