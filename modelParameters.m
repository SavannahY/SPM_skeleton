function param = modelParameters(Nr)

    % Initialize Model Parameters (SI Units)
    % Framework for geometric parameter identification
    % FOR ELECTROCHEMICAL MODELS: 
    % Cell Name: INR21700
    % Manufacturer: LG Chem
    % Rated Capacity = 4.85 A-h
    % Nominal Voltage = 3.63V 
    % Cathode Chemistry = NMC / Anode Chemistry = Graphite
    
    param.Rg = 8.314;       % Gas constant [J/(mol·K)]
    param.F = 96487;        % Faraday's constant [C/mol]
    param.alpha_cell = 0.5; % Anode/Cathode transfer coefficient [-]
    param.ce0 = 1000;       % Average electrolyte concentration [mol/m³]
    param.capacity = 4.85;  % Nominal cell capacity [A·h]
    
    param.c_n_max = 29583;  % Maximum anode concentration [mol/m³] (from Chen et al. 2020)
    param.c_p_max = 51765;  % Maximum cathode concentration [mol/m³] (from Chen et al. 2020)
    
    param.V_upper_lim = 4.3; % Upper voltage limit [V]
    param.V_lower_lim = 2.4; % Lower voltage limit [V]

    % known parameters
    param.Acell = 0.103675; % Cell active area [m²]
    param.Ln    = 85.2e-6;  % Anode thickness [m]
    param.Lp    = 75.6e-6;  % Cathode thickness [m]
    
    param.csn_max = 29583;  % Maximum anode solid concentration [mol/m³]
    param.csp_max = 51765;  % Maximum cathode solid concentration [mol/m³]

    param.ce = 1000;        % Electrolyte concentration [mol/m³]

    param.Rsn = 3.6e-6;     % Anode particle radius [m]
    param.Rsp = 6.25e-6;    % Cathode particle radius [m]

    param.theta_n_0 = 0.027935145151956345;   % Anode stoichiometry at 0% SOC [-]
    param.theta_p_0 = 0.9066422666349141;     % Cathode stoichiometry at 0% SOC [-]
    param.theta_n_100 = 0.9219907588763882;   % Anode stoichiometry at 100% SOC [-]
    param.theta_p_100 = 0.267963982283029;    % Cathode stoichiometry at 100% SOC [-]

    param.R0 = 0.02983;     % Ohmic resistance [Ω]
    
    % define kinetic parameters 
    param.Dsn  = 9.2299e-15;  % Anode solid diffusivity [m²/s]
    param.epsn = 0.7788;      % Anode active material volume fraction [-]
    param.k0n  = 5.2828e-06;  % Anode reaction rate constant [m^2.5/(mol^0.5·s)]
    param.Dsp  = 4.5900e-15;  % Cathode solid diffusivity [m²/s]
    param.epsp = 0.7022;      % Cathode active material volume fraction [-]
    param.k0p  = 4.2411e-06;  % Cathode reaction rate constant [m^2.5/(mol^0.5·s)]
    
    % Grid interval size for anode and cathode particle
    param.dr_n = 1.0./(Nr-1);    % Anode dimensionless radial step size [-]
    param.dr_p = 1.0./(Nr-1);    % Cathode dimensionless radial step size [-]
    
    % store grid
    param.r_n = linspace(0.0, 1.0, Nr)';  % Anode dimensionless radial grid [-]
    param.r_p = linspace(0.0, 1.0, Nr)';  % Cathode dimensionless radial grid [-]
    
    % store half-grid
    param.r_half_n = 0.5*(param.r_n(1:end-1) + param.r_n(2:end));  % Anode half-grid points [-]
    param.r_half_p = 0.5*(param.r_p(1:end-1) + param.r_p(2:end));  % Cathode half-grid points [-]

end
