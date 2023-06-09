function [config] = load_config
%LOAD_CONFIG loads the scenario configuration
    config.CNR = 10; % clutter-to-noise ratio in dB
    config.TCR = 5; % target-to-noie ratio in dB
    config.N = 15; % number of temporal measurements
    config.L = 5; % number of active subcarriers
    config.Deltaf = 20e6; % subcarrier spacing in Hz
    config.T = 1/config.Deltaf; % OFDM TX pulse duration in seconds
    config.B = config.L*config.Deltaf; % OFDM bandwidth B in Hz
    config.varX = 1; % variance of the complex scattering coefficients (target)
    config.varSigma_c = 1; % variance for the generation of the clutter response covariance matrix Sigma_c
    config.target_velocity = [10, 10]'; % target velocity in m/s
    config.target_angle = [10, 10]'; % TODO: check best value for this parameter
    config.eta_true = [10, 10]';
    config.fc = 8e9; % carrier frequency in Hz
    config.var_noise = 5; % std measurement noise e(t)
    config.gamma = 3; % gamma detection threshold
end