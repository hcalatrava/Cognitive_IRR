% ------------------------------------------------------------------------------
% Cognitive Interference Resilient Radar (Cognitive_IRR) 
% Author: Helena Calatrava
% Affiliation: Northeastern University, Boston, United Sates
% Date: July 2023
%
% MAIN calculates and plots ROC curves for different configurations
%   - different number of OFDM subcarriers
%   - different values of SNR
% ------------------------------------------------------------------------------

%% Initialize scenario
close all,
clear all,
clc,
config = load_config(1); % Load configuration

%% Test for different number of subcarriers L
% Results are plotted on the bottom subfigure

% Configuration
rng('default') % Initialize seed
% Get LOS vector
LOS_vector = get_LOS_vector(config);

% Create figure
figure()
subplot(2,1,2)

% Loop through different number of subcarriers L
L_vector = [2,4,6]; % values of L for which a ROC curve is computed
gamma_vector = config.gamma_vector % values of \gamma tested to compute the ROC curve
for idx_L = 1:length(L_vector)

    % Generate OFDM measurements
    config.L = L_vector(idx_L); % update number of subcarriers
    % Prior probabilities (we generate measurements following hypothesis 0
    % and hypothesis 1)
    p_H0 = 0.5; p_H1 = 1 - p_H0;
    labels = rand(config.N_mc,1) <= p_H1; % 1 when H1 and 0 when H0
    
    % OFDM weights
    % These weights can be updated in the Adaptive Waveform Design (AWD) module
    A = eye(config.L);
    
    % Doppler information
    Phi_c = ones(config.L, config.N); % LxN matrix containing the Doppler information (clutter)
    % Phi_c is an all-ones matrix given that the clutter is considered to be static
    Phi_t = get_doppler_mat(config, LOS_vector); % LxN matrix containing the Doppler information through the parameter config.eta (target)
    
    % Initialize metrics to compute the ROC curve
    Ppp = zeros(length(gamma_vector), 1); % probability of positive detector and positive label p(D=1|H1)
    Ppn = zeros(length(gamma_vector), 1); % p(D=1|H0)
    Pnp = zeros(length(gamma_vector), 1); % p(D=0|H1)
    Pe = zeros(length(gamma_vector), 1); % probability of error (Pe)
    Nl0 = [length(find(labels==0))]; % number of samples from class 0
    Nl1 = [length(find(labels==1))]; % number of samples from class 1
    
    % For each Monte Carlo realization, generate OFDM measurements
    % according to the hypothesis indicated by the label, and calculate
    % GLRT.
    GLRT = zeros(config.N_mc, 1);
    for idx_mc = 1:config.N_mc
        if labels(idx_mc, 1) == 1
            hypothesis = 'H1';
        else % labels(1, idx_mc) == 1
            hypothesis = 'H0';
        end
        Y = get_measurements_vmulti(config, hypothesis, A, Phi_t, Phi_c);

        % MLE of target coefficients X_t
        X_t_hat = diag(diag(A\Y*Phi_t'*inv(Phi_t*Phi_t')));

        % GLRT
        GLRT(idx_mc, 1)= det(Y*Y')/det((Y-A*X_t_hat*Phi_t)*(Y-A*X_t_hat*Phi_t)');
    end
    
    % For each \gamma, calculate the probability of falsa alarm (FPR) and
    % the probability of detection (TPR)
    for gamma_idx = 1:length(gamma_vector)

        gamma = gamma_vector(gamma_idx);
        detected = real(GLRT)>gamma;

        % For each value of gamma, compute the FPR (false positive rate) and
        % the TPR (true positive rate)
        % Number of True/False decisions for each True/False label
        Nd0l1 = [length(find(detected==0 & labels==1))];
        Nd1l0 = [length(find(detected==1 & labels==0))];
        Nd1l1 = [length(find(detected==1 & labels==1))];
        % False negative: P(D=0 | L=1;gamma)
        Pnp(gamma_idx, 1) = Nd0l1/Nl1;
        % False positive: P(D=1 | L=0;gamma)
        Ppn(gamma_idx, 1) = Nd1l0/Nl0;
        % True positive: P(D=1 | L=1;gamma)
        Ppp(gamma_idx, 1) = Nd1l1/Nl1;
        % Probability of error: P(error; gamma)
        Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
        disp(['Results have been computed for gamma=', num2str(gamma), ' with TPR=', num2str(Ppp(gamma_idx, 1))])
    end

    % Plot ROC for this configuration
    semilogx(Ppn(:,1), Ppp(:,1), 'displayname', ['L = ', num2str(L_vector(idx_L))])
    hold on,
end

legend('show', 'location', 'southeast')
grid on
title('ROC')
axis([0 1 0 1])
xlim([0.005, 1]) % same x axis limit as in paper #2
xlabel('FPR P(D=1 | H0)')
ylabel('TPR P(D=1 | H1)')

%% Test for different values of SNR
% Results are plotted on the top subfigure

% Configuration
rng('default') % Initialize seed
% Get LOS vector
LOS_vector = get_LOS_vector(config);
config.L = 4;

% Create figure
subplot(2,1,1)

% Loop through different values of SNR
SNR_vector = [-5, -10, -15]; % in dB
gamma_vector = config.gamma_vector % values of \gamma tested to compute the ROC curve
for idx_SNR = 1:length(SNR_vector)
    
    % Generate OFDM measurements
    config.SNR = SNR_vector(idx_SNR);
    % Prior probabilities (we generate measurements following hypothesis 0
    % and hypothesis 1)
    p_H0 = 0.5; p_H1 = 1 - p_H0;
    labels = rand(config.N_mc,1) <= p_H1; % 1 when H1 and 0 when H0
    
    % OFDM weights
    % These weights can be updated in the Adaptive Waveform Design (AWD) module
    A = eye(config.L);
    
    % Doppler information
    Phi_c = ones(config.L, config.N); % LxN matrix containing the Doppler information (clutter)
    % Phi_c is an all-ones matrix given that the clutter is considered to be static
    Phi_t = get_doppler_mat(config, LOS_vector); % LxN matrix containing the Doppler information through the parameter config.eta (target)
    
    % Initialize metrics to compute the ROC curve
    Ppp = zeros(length(gamma_vector), 1); % probability of positive detector and positive label p(D=1|H1)
    Ppn = zeros(length(gamma_vector), 1); % p(D=1|H0)
    Pnp = zeros(length(gamma_vector), 1); % p(D=0|H1)
    Pe = zeros(length(gamma_vector), 1); % probability of error (Pe)
    Nl0 = [length(find(labels==0))]; % number of samples from class 0
    Nl1 = [length(find(labels==1))]; % number of samples from class 1
    
    % For each Monte Carlo realization, generate OFDM measurements
    % according to the hypothesis indicated by the label, and calculate
    % GLRT.
    GLRT = zeros(config.N_mc, 1);
    for idx_mc = 1:config.N_mc
        if labels(idx_mc, 1) == 1
            hypothesis = 'H1';
        else % labels(1, idx_mc) == 1
            hypothesis = 'H0';
        end
        Y = get_measurements_vmulti(config, hypothesis, A, Phi_t, Phi_c);

        % MLE of target coefficients X_t
        X_t_hat = diag(diag(A\Y*Phi_t'*inv(Phi_t*Phi_t')));

        % GLRT
        GLRT(idx_mc, 1)= det(Y*Y')/det((Y-A*X_t_hat*Phi_t)*(Y-A*X_t_hat*Phi_t)');
    end
    
    % For each \gamma, calculate the probability of falsa alarm (FPR) and
    % the probability of detection (TPR)
    for gamma_idx = 1:length(gamma_vector)

        gamma = gamma_vector(gamma_idx);
        detected = real(GLRT)>gamma;

        % For each value of gamma, compute the FPR (false positive rate) and
        % the TPR (true positive rate)
        % Number of True/False decisions for each True/False label
        Nd0l1 = [length(find(detected==0 & labels==1))];
        Nd1l0 = [length(find(detected==1 & labels==0))];
        Nd1l1 = [length(find(detected==1 & labels==1))];
        % False negative: P(D=0 | L=1;gamma)
        Pnp(gamma_idx, 1) = Nd0l1/Nl1;
        % False positive: P(D=1 | L=0;gamma)
        Ppn(gamma_idx, 1) = Nd1l0/Nl0;
        % True positive: P(D=1 | L=1;gamma)
        Ppp(gamma_idx, 1) = Nd1l1/Nl1;
        % Probability of error: P(error; gamma)
        Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
        disp(['Results have been computed for gamma=', num2str(gamma), ' with TPR=', num2str(Ppp(gamma_idx, 1))])
    end

    % Plot ROC curve for this configuration
    semilogx(Ppn(:,1), Ppp(:,1), 'displayname', ['SNR = ', num2str(SNR_vector(idx_SNR))])
    hold on,
end
legend('show', 'location', 'southeast')
grid on
title('ROC')
axis([0 1 0 1])
xlim([0.005, 1]) % same x axis limit as in paper #2
xlabel('FPR P(D=1 | H0)')
ylabel('TPR P(D=1 | H1)')
