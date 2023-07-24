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

if config_id == 0
    % General config
    config.var_noise = 1; % std measurement noise e(t) (known)
    config.N = 15; % number of temporal measurements
    config.N_mc = 20e3;
    
    % OFDM
    config.L = 1; % number of active subcarriers
    config.Deltaf = 20e6; % subcarrier spacing in Hz
    config.T = 1/config.Deltaf; % OFDM TX pulse duration in seconds
    config.B = config.L*config.Deltaf; % OFDM bandwidth B in Hz
    config.fc = 8e9; % carrier frequency in Hz
    
    % Target
    config.TCR = 5; % target-to-noie ratio in dB
    config.target_velocity = [10, 10]'; % target velocity in m/s
    config.varX = 1; % variance of the complex scattering coefficients (target)    
    % Line-of-Sight (LOS) vector
    % target in range cell centered at config.target_dist_north meters North and
    % config.target_dist_east meters East w.r.t. the radar
    config.target_dist_north = 2e3; % in meters
    config.target_dist_east = 5; 
    
    % Clutter
    config.CNR = 10; % clutter-to-noise ratio in dB
    config.varSigma_c = 1; % variance for the generation of the clutter response covariance matrix Sigma_c

elseif config_id == 1
    % General config
    config.N_mc = 20e3; % Number of Monte Carlo realizations
    config.N = 15; % Number of temporal measurements N
    config.gamma_vector = 0:0.01:10; % values of \gamma tested to compute the ROC curve

    % OFDM
    config.L = 1; % number of active subcarriers
    config.Deltaf = 20e6; % subcarrier spacing in Hz
    config.T = 50e-9; % OFDM TX pulse duration in seconds
    config.B = 100e6; % OFDM bandwidth B in Hz
    config.fc = 1e9; % carrier frequency in Hz
    
    % Target
    config.SNR = -10; % signal-to-noise ratio in dB
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
end % end if config_id

end