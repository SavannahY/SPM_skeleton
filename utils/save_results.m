function save_results(output_dir, all_data, config)
    %SAVE_RESULTS Save simulation results to MAT-file
    %
    % INPUTS:
    %   output_dir - Output directory path [string]
    %   all_data   - Structure containing all simulation results with units
    %   config     - Configuration structure with discretization and solver info
    %
    % OUTPUTS:
    %   None (saves data to disk)
    if nargin < 2
        error('Configuration structure is required to determine output directory.');
    end
    % Build output file name
    % Generate descriptive filename
    matlab_outname = [output_dir '/matlab_sim_results_' ...
                        convertStringsToChars(config.method) ...
                        '_Nr_', num2str(config.Nr, '%.0f'), ...
                        '_' config.dchg_type_str, ...
                        '_int_', config.integrator, ...
                        '_outmat.mat'];

    % Save all data to MAT file
    save(matlab_outname, '-struct', 'all_data');

    % Display output location
    disp(['Matlab output: ' matlab_outname])
end