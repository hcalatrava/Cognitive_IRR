function [config] = load_config(config_id)
% ------------------------------------------------------------------------------
% Cognitive Interference Resilient Radar (Cognitive_IRR) 
% Author: Helena Calatrava
% Affiliation: Northeastern University, Boston, United Sates
% Date: July 2023
%
% LOAD_CONFIG loads configuration for specified configuration ID
%   ID | Description
%    0 | Values from paper #1
%    1 | Values from paper #2
%
% References:
% 1 - Target Detection in Clutter Using Adaptive OFDM Radar
%     Authors: Satyabrata Sen, Arye Nehorai
% 2 - Adaptive OFDM Radar for Target Detection in Multipath Scenarios
%     Authors: Satyabrata Sen, Arye Nehorai
% ------------------------------------------------------------------------------

if config_id == 3
            % Plot config
    config.plot_color = ["blue", "red", "black"];
    config.plot_linestyle = ["-.", "--", "-"];

    % General config
    config.var_noise = 1; % std measurement noise e(t) (known)
    config.N = 15; % number of temporal measurements
    config.N_mc = 2e3;
    config.c = 3e08;
    config.gamma_vector = 0:0.01:10; % values of \gamma tested to compute the ROC curve
    
    % OFDM
    config.L = 1; % number of active subcarriers
    config.Deltaf = 20e6; % subcarrier spacing in Hz
    config.T = 1/config.Deltaf; % OFDM TX pulse duration in seconds
    config.B = config.L*config.Deltaf; % OFDM bandwidth B in Hz
    config.fc = 2e9; % carrier frequency in Hz

    %% 
    config.T_PRI = 20e-6;
    config.T = config.T_PRI/300; % from provided code
    config.Deltaf = 1/config.T;
    
    % Target
    config.SNR_predefined = -10; % target-to-noie ratio in dB
    config.target_velocity = [0.5; 2]; % target velocity in m/s
    config.varX = 1; % variance of the complex scattering coefficients (target)    

    config.target_dist_north = 0.2e3; % in meters
    config.target_dist_east = 5; 
    
    % ML estimation
    config.clairvoyant = 0;

    % Adaptive Waveform Design Module
    config.enable_awd = 1;

    % Adversarial effect on relative doppler shift beta
    config.add_doppler_jamming = 1;

    %
    config.known_velocity = 1;


elseif config_id == 0
        % Plot config
    config.plot_color = ["blue", "red", "black"];
    config.plot_linestyle = ["-.", "--", "-"];

    % General config
    config.var_noise = 1; % std measurement noise e(t) (known)
    config.N = 15; % number of temporal measurements
    config.N_mc = 2e3;
    config.c = 3e08;
    config.gamma_vector = 0:0.01:10; % values of \gamma tested to compute the ROC curve
    
    % OFDM
    config.L = 1; % number of active subcarriers
    config.Deltaf = 20e6; % subcarrier spacing in Hz
    config.T = 1/config.Deltaf; % OFDM TX pulse duration in seconds
    config.B = config.L*config.Deltaf; % OFDM bandwidth B in Hz
    config.fc = 2e9; % carrier frequency in Hz

    %% 
    config.T_PRI = 20e-6;
    config.T = config.T_PRI/300; % from provided code
    config.Deltaf = 1/config.T;
    
    % Target
    config.SNR_predefined = -10; % target-to-noie ratio in dB
    config.target_velocity = [0.5; 2]; % target velocity in m/s
    config.varX = 1; % variance of the complex scattering coefficients (target)    

    config.target_dist_north = 2e3; % in meters
    config.target_dist_east = 5; 
    
    % ML estimation
    config.clairvoyant = 1;

    % Adaptive Waveform Design Module
    config.enable_awd = 1;

    % 
    config.known_velocity = 1; % 0 if we want to maximize the GLRT w.r.t. the target velocity
elseif config_id == -1

    % Plot config
    config.plot_color = ["blue", "red", "black"];
    config.plot_linestyle = ["-.", "--", "-"];

    % General config
    config.var_noise = 1; % std measurement noise e(t) (known)
    config.N = 15; % number of temporal measurements
    config.N_mc = 10e3;
    config.c = 3e08;
    config.gamma_vector = 0:0.01:20; % values of \gamma tested to compute the ROC curve
    
    % OFDM
    config.L = 1; % number of active subcarriers
    config.Deltaf = 20e6; % subcarrier spacing in Hz
    config.T = 1/config.Deltaf; % OFDM TX pulse duration in seconds
    config.B = config.L*config.Deltaf; % OFDM bandwidth B in Hz
    config.fc = 2e9; % carrier frequency in Hz

    %% 
    config.T_PRI = 20e-6;
    config.T = config.T_PRI/300; % from provided code
    config.Deltaf = 1/config.T;
    
    % Target
    config.SNR_predefined = -10; % target-to-noie ratio in dB
    config.target_velocity = [10; 10]; % target velocity in m/s
    config.varX = 1; % variance of the complex scattering coefficients (target)    
    % Line-of-Sight (LOS) vector
    % target in range cell centered at config.target_dist_north meters North and
    % config.target_dist_east meters East w.r.t. the radar
    config.target_dist_north = 2e3; % in meters
    config.target_dist_east = 5; 
    
    % Clutter
    config.CNR = -2; % clutter-to-noise ratio in dB
    config.varSigma_c = 1; % variance for the generation of the clutter response covariance matrix Sigma_c

    % ML estimation
    config.clairvoyant = 1

    % Adaptive Waveform Design Module

    
elseif config_id == 1
    % General config
    config.N_mc = 20e3; % Number of Monte Carlo realizations
    config.N = 15; % Number of temporal measurements N
    config.gamma_vector = 0:0.01:10; % values of \gamma tested to compute the ROC curve
    config.c = 3e8;

    % OFDM
    config.L = 4; % number of active subcarriers
    config.Deltaf = 20e6; % subcarrier spacing in Hz
    config.T = 50e-9; % OFDM TX pulse duration in seconds
    config.B = 100e6; % OFDM bandwidth B in Hz
    config.fc = 1e9; % carrier frequency in Hz
    config.T_PRI = 20e-6; % pulse repetition interval
    
    % Target
    config.SNR = -15; % signal-to-noise ratio in dB
    config.target_velocity = [7.07, 7.07]'; % target velocity in m/s
    config.varX = 1; % variance of the complex scattering coefficients (target)
    % Line-of-Sight (LOS) vector
    % target in range cell centered at config.target_dist_north meters North and
    % config.target_dist_east meters East w.r.t. the radar
    config.target_dist_north = 2e3; % in meters
    config.target_dist_east = 5; % 
    
    % Adaptive Waveform Design (AWD)
    config.mu = 0.5;
    config.E_A = 1;
elseif config_id == 5%'dynamical'
        % Plot config
%     config.plot_color = ["blue", "red", "black"];
%     config.plot_linestyle = ["-.", "--", "-"];

    % General config
    config.var_noise = 1; % std measurement noise e(t) (known)
    config.N = 15; % number of temporal measurements
    config.N_mc = 2e3;
    config.c = 3e08;
    config.gamma_vector = 0:0.01:10; % values of \gamma tested to compute the ROC curve
    
    % OFDM
    config.L = 3; % number of active subcarriers
    config.Deltaf = 20e6; % subcarrier spacing in Hz
    config.T = 1/config.Deltaf; % OFDM TX pulse duration in seconds
    config.B = config.L*config.Deltaf; % OFDM bandwidth B in Hz
    config.fc = 2e9; % carrier frequency in Hz

    %% 
    config.T_PRI = 20e-6;
    config.T = config.T_PRI/300; % from provided code
    config.Deltaf = 1/config.T;
    
    % Target
    config.SNR_predefined = -10; % target-to-noie ratio in dB
    config.target_velocity = [7.07, 4]; % target velocity in m/s
    config.varX = 1; % variance of the complex scattering coefficients (target)    

    config.target_dist_north = 5e1; % in meters
    config.target_dist_east = 1e1; 
    
    % ML estimation
    config.clairvoyant = 0;

    % Adaptive Waveform Design Module
    config.enable_awd = 1;
elseif config_id == 6 %analysis_system' % main_sep26.m script
    % Plotting configuration
    config.plot_color = ["blue", "red", "black"];
    config.plot_linestyle = ["-.", "--", "-"];

    % General config
%     config.var_noise = 1; % std measurement noise e(t) (known)
    config.N_mc = 2e3;
    config.c = 3e08;
    config.gamma_vector = 0:0.01:20; % values of \gamma threshold tested to compute the ROC curve
    
    % OFDM
    config.L = 3; % number of active subcarriers
    config.Deltaf = 20e6; % subcarrier spacing in Hz
    config.T = 1/config.Deltaf; % OFDM TX pulse duration in seconds
    config.B = config.L*config.Deltaf; % OFDM bandwidth B in Hz
    config.fc = 2e9; % carrier frequency in Hz
    config.N = 15; % number of temporal measurements available per data chunk
    %
    config.T_PRI = 20e-6;
    config.T = config.T_PRI/300; % from provided code
    config.Deltaf = 1/config.T;
%     
    % Target
    config.SNR_predefined = -15; % target-to-noie ratio in dB
    config.target_velocity = [7.07, 7.07]'; % target velocity in m/s
    config.varX = 1; % variance of the complex scattering coefficients (target)    
    config.target_dist_north = 5e1; % distance of target from radar in the north direction
    config.target_dist_east = 1e1; % disance of target from radar in the east direction
    
    % Parameter estimation
    % Velocity is unknown so we maximize GLRT w.r.t. velocity
    % \tau_0 is assumed to be known so we calculate T_rangecell at start
    % covariance or factor \rho ???
    config.clairvoyant = 0; % 1 if target coefficients are known by the receiver (unrealistic)
    config.known_velocity = 0; % 0 if we want to maximize the GLRT w.r.t. the target velocity

    % Adaptive Waveform Design Module
    config.enable_awd = 0;
end % end if config_id

end