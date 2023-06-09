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

%% Test Performance (different gamma threshold values)
rng('default')

% Prior probabilities
p_H0 = 0.2; p_H1 = 1 - p_H0;
labels = rand(config.N_mc,1) <= p_H1; % 1 when H1 and 0 when H0

% OFDM weights
% These weights can be updated in the Adaptive Waveform Design (AWD) module
A = eye(config.L); % LxL complex diagonal matrix ontaining the transmitted weights

% Doppler information
Phi_c = ones(config.L, config.N); % LxN matrix containing the Doppler information (clutter)
% Phi_c is an all-ones matrix given that the clutter is considered to be static
Phi_t = get_doppler_mat(config); % LxN matrix containing the Doppler information through the parameter config.eta (target)

gamma_vector = 0:0.1:10;

Ppp = zeros(length(gamma_vector), 1);
Ppn = zeros(length(gamma_vector), 1);
Pnp = zeros(length(gamma_vector), 1);
Pe = zeros(length(gamma_vector), 1);
Nl0 = [length(find(labels==0))]; % number of samples from class 0
Nl1 = [length(find(labels==1))]; % number of samples from class 1

GLRT = zeros(config.N_mc, 1);
% Get measurements for each label (as many as Monte Carlo realizations)
for idx_mc = 1:config.N_mc
    if labels(idx_mc, 1) == 1
        hypothesis = 'H1';
    else % labels(1, idx_mc) == 1
        hypothesis = 'H0';
    end
    Y = get_measurements(config, hypothesis, A, Phi_t, Phi_c);
    % MLE of coefficients X_t
    X_t_hat = diag(diag(A\Y*Phi_t'*inv(Phi_t*Phi_t')));
    % GLRT
    GLRT(idx_mc, 1)= det(Y*Y')/det((Y-A*X_t_hat*Phi_t)*(Y-A*X_t_hat*Phi_t)');
end

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

% Plot ROC
figure()
title('ROC')
axis([0 1 0 1])
plot(Ppn(:,1), Ppp(:,1))
xlabel('P(D=1 | H0)')
ylabel('P(D=1 | H1)')
legend('ROC')

%% Test for different number of subcarriers
rng('default')

figure()

L_vector = [1, 3, 5, 10]; % vector with different number of subcarriers
gamma_vector = 0:0.01:10;

for idx_L = 1:length(L_vector)
    config.L = L_vector(idx_L);

    % Prior probabilities
    p_H0 = 0.2; p_H1 = 1 - p_H0;
    labels = rand(config.N_mc,1) <= p_H1; % 1 when H1 and 0 when H0
    
    % OFDM weights
    % These weights can be updated in the Adaptive Waveform Design (AWD) module
    A = eye(config.L); % LxL complex diagonal matrix ontaining the transmitted weights
    
    % Doppler information
    Phi_c = ones(config.L, config.N); % LxN matrix containing the Doppler information (clutter)
    % Phi_c is an all-ones matrix given that the clutter is considered to be static
    Phi_t = get_doppler_mat(config); % LxN matrix containing the Doppler information through the parameter config.eta (target)
    
    Ppp = zeros(length(gamma_vector), 1);
    Ppn = zeros(length(gamma_vector), 1);
    Pnp = zeros(length(gamma_vector), 1);
    Pe = zeros(length(gamma_vector), 1);
    Nl0 = [length(find(labels==0))]; % number of samples from class 0
    Nl1 = [length(find(labels==1))]; % number of samples from class 1
    
    GLRT = zeros(config.N_mc, 1);
    % Get measurements for each label (as many as Monte Carlo realizations)
    for idx_mc = 1:config.N_mc
        if labels(idx_mc, 1) == 1
            hypothesis = 'H1';
        else % labels(1, idx_mc) == 1
            hypothesis = 'H0';
        end
        Y = get_measurements(config, hypothesis, A, Phi_t, Phi_c);
        % MLE of coefficients X_t
        X_t_hat = diag(diag(A\Y*Phi_t'*inv(Phi_t*Phi_t')));
        % GLRT
        GLRT(idx_mc, 1)= det(Y*Y')/det((Y-A*X_t_hat*Phi_t)*(Y-A*X_t_hat*Phi_t)');
    end
    
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

    % Plot ROC
    plot(Ppn(:,1), Ppp(:,1), 'displayname', ['L = ', num2str(L_vector(idx_L))])
    hold on,
end
legend('show', 'location', 'southeast')
grid on
title('ROC')
axis([0 1 0 1])
xlabel('FPR P(D=1 | H0)')
ylabel('TPR P(D=1 | H1)')


%% Adaptive Waveform Design (AWD)
