% ------------------------------------------------------------------------------
% Cognitive Interference Resilient Radar (Cognitive_IRR) 
% Helena Calatrava
%
% June 2023
% ------------------------------------------------------------------------------

%% Initialize Scenario
close all,
clear all,
config = load_config;

%% Build Measurements Model Y

[Y, A, Phi_t] = get_measurements(config, 'H0');

%% Generalized Likelihood Ratio Test (GLRT)

% MLE of coefficients X_t
X_t_hat = diag(diag(A\Y*Phi_t'*inv(Phi_t*Phi_t')));

% GLRT
GLRT = det(Y*Y')/det((Y-A*X_t_hat*Phi_t)*(Y-A*X_t_hat*Phi_t)');

%% Test Performance (P_D - Probability of Detection) 
detected = 0;
for idx_mc = 1:config.N_mc
    % Get observations
    [Y, A, Phi_t] = get_measurements(config, 'H1');

    % MLE of coefficients X_t
    X_t_hat = diag(diag(A\Y*Phi_t'*inv(Phi_t*Phi_t')));

    % GLRT
    GLRT(idx_mc) = det(Y*Y')/det((Y-A*X_t_hat*Phi_t)*(Y-A*X_t_hat*Phi_t)');
    if real(GLRT) > 1.8
        detected = detected + 1;
    end
end
disp(['Probability of detection (P_D) = ',num2str(round(detected/config.N_mc*100,2)), '%']);

%% Test Performance (P_FA - Probability of False Alarm) 
detected = 0;
for idx_mc = 1:config.N_mc
    % Get observations
    [Y, A, Phi_t] = get_measurements(config, 'H0');

    % MLE of coefficients X_t
    X_t_hat = diag(diag(A\Y*Phi_t'*inv(Phi_t*Phi_t')));

    % GLRT
    GLRT = det(Y*Y')/det((Y-A*X_t_hat*Phi_t)*(Y-A*X_t_hat*Phi_t)');
    if real(GLRT) > 1.8
        detected = detected + 1;
    end
end
disp(['Probability of detection (P_FA) = ',num2str(round(detected/config.N_mc*100,2)), '%']);

%% Test Performance (different gamma threshold values)
%for gamma = 

%% Adaptive Waveform Design (AWD)
