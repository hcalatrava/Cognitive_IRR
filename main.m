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
%% Initialize scenario (paper 0) - ADVERSARIAL EFFECTS
% close all,
% clear all,
clc,
config = load_config(0); % Load configuration

%%% Test for different number of subcarriers L
% Results are plotted on the bottom subfigure

% Configuration
rng('default') % Initialize seed
% Get LOS vector
LOS_vector = get_LOS_vector(config);
% r = sqrt(config.target_dist_north^2 + config.target_dist_east^2); % distance between target and radar (meters)
% LOS_vector = [config.target_dist_east, config.target_dist_north]/r;

% Create figure
figure()
% subplot(2,1,1)

plot_strings = ["Safe Radar", "Radar under soft attack", "Radar under strong attack"];
plot_strings = ["Safe Radar", "Radar under soft attack", "Protected radar under soft attack"];

% Loop through different number of subcarriers L
config.L = 3;
config.jamming_coeff_vector = [0, 1, 1]; % values of L for which a ROC curve is computed
gamma_vector = config.gamma_vector; % values of \gamma tested to compute the ROC curve
for idx_jamming = 1:length(config.jamming_coeff_vector)

    % Generate OFDM measurements
    config.jamming_coeff = config.jamming_coeff_vector(idx_jamming); % update number of subcarriers
%     config.Deltaf = config.B/(config.L);
%     config.T = 1/config.Deltaf;

    % OFDM weights
    % These weights can be updated in the Adaptive Waveform Design (AWD) module
    A = eye(config.L)/sqrt(config.L);
    A = eye(config.L);
    if idx_jamming ==3
        A = [0.1 0 0; 0 2 0; 0 0 2];
    end
    
    % Doppler information
    beta = 2*config.target_velocity'*LOS_vector'/physconst('LightSpeed');
    omega_dopp = 2*pi*beta*(config.fc + (0:(config.L-1))'*config.Deltaf);
    doppler_mat = exp(1i*omega_dopp*(1:config.N));
    %Phi_t = get_doppler_mat(config, LOS_vector); % LxN matrix containing the Doppler information through the parameter config.eta (target)
    Phi_t = doppler_mat;

    % Doppler information (JAMMER)
    beta = 2*config.target_velocity'*LOS_vector'/physconst('LightSpeed');
    omega_dopp = 2*pi*2.4e9;
    doppler_mat = exp(1i*omega_dopp*ones(config.L,config.N));
    %Phi_t = get_doppler_mat(config, LOS_vector); % LxN matrix containing the Doppler information through the parameter config.eta (target)
    Phi_jammer = doppler_mat;
    Phi_jammer(2:3,:) = zeros(2,config.N);
%     Phi_jammer = zeros(config.L,config.N);

    % Initialize metrics to compute the ROC curve
    Ppp = zeros(length(gamma_vector), 1); % probability of positive detector and positive label p(D=1|H1)
    Ppn = zeros(length(gamma_vector), 1); % p(D=1|H0)
    
    % For each Monte Carlo realization, generate OFDM measurements
    % according to the hypothesis indicated by the label, and calculate
    % GLRT.
    GLRT_H0 = zeros(config.N_mc, 1);
    GLRT_H1 = zeros(config.N_mc, 1);
    for idx_mc = 1:config.N_mc
        
        % Generate target coefficients
        X_H0 = zeros(config.L); % all zeros
        X_H1 = diag(sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1))); % LxL complex diagonal matrix repressenting the scattering coefficients (target)
        
        % Get Sigma_c
        sum_n = 0; % TODO: clean this loop
        for n = 1:config.N
            sum_n = sum_n + (A*X_H1*Phi_t(:,n))'*(A*X_H1*Phi_t(:,n));
        end
        tcr_lin = 10^(config.TCR/10);
        alpha = 1/config.N*sum_n/tcr_lin;
        R_sqrt = (randn(config.L) + 1i*randn(config.L))*(1/sqrt(2));
        R = R_sqrt'*R_sqrt;
        k = sqrt(alpha/trace(A*R*A'));
        Sigma_c = k^2*R;

        % Get clutter coefficients and doppler info
        XcPhic = get_noise_mat(config, Sigma_c);

        % Get noise matrix
        E = zeros(config.L,config.N);
        cnr_lin = 10^(config.CNR/10);
        var_noise = trace(A*Sigma_c*A')/config.L/cnr_lin;
        for i=1:config.N
            % We generate N vectors e(t), each one satisfying:
            % E[e(t)*e^H(t)] = Sigma_c
            e = sqrt(var_noise/2)*(randn(config.L,1)+1i*randn(config.L,1));
            E(:,i) = e;
        end

        % Get OFDM measurments
        Y_H0 = A*XcPhic +E;
        Y_H1 = A*X_H1*Phi_t + A*XcPhic + E;
        Y_H1 = A*X_H1*Phi_t + config.jamming_coeff*A*X_H1*Phi_jammer+ A*XcPhic + E;

        % MLE of target coefficients X under hypothesis H1
        X_hat = diag(diag(A\Y_H1*Phi_t'*inv(Phi_t*Phi_t')));
%         X_hat = X_H1;
%         X_hat = X_H1/10;
%         Pi_Phi = Phi_t'*inv(Phi_t*Phi_t')*Phi_t;
%         Pi_Phi_dag = eye(config.N) - Pi_Phi;
%         G = Y_H1*Pi_Phi_dag*Y_H1';
%         G = G - diag(G) + real(diag(G));
%         S = chol(G)*A;
%         Gamma = kron(S',chol(G)*A)
%         X_hat = inv(Gamma'*Gamma)*Gamma'*(chol(G)*Y_H1*Pi_Phi)

        % GLRT
        GLRT_H0(idx_mc, 1)= det(Y_H0*Y_H0')/det((Y_H0-A*X_hat*Phi_t)*(Y_H0-A*X_hat*Phi_t)');
        GLRT_H1(idx_mc, 1)= det(Y_H1*Y_H1')/det((Y_H1-A*X_hat*Phi_t)*(Y_H1-A*X_hat*Phi_t)');
    end
    
%     % For each \gamma, calculate the probability of falsa alarm (FPR) and
%     % the probability of detection (TPR)
%     for gamma_idx = 1:length(gamma_vector)
% 
%         gamma = gamma_vector(gamma_idx);
%         detected = real(GLRT)>gamma;
% 
%         % For each value of gamma, compute the FPR (false positive rate) and
%         % the TPR (true positive rate)
%         % Number of True/False decisions for each True/False label
%         Nd0l1 = [length(find(detected==0 & labels==1))];
%         Nd1l0 = [length(find(detected==1 & labels==0))];
%         Nd1l1 = [length(find(detected==1 & labels==1))];
%         % False negative: P(D=0 | L=1;gamma)
%         Pnp(gamma_idx, 1) = Nd0l1/Nl1;
%         % False positive: P(D=1 | L=0;gamma)
%         Ppn(gamma_idx, 1) = Nd1l0/Nl0;
%         % True positive: P(D=1 | L=1;gamma)
%         Ppp(gamma_idx, 1) = Nd1l1/Nl1;
%         % Probability of error: P(error; gamma)
%         Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
%         disp(['Results have been computed for gamma=', num2str(gamma), ' with TPR=', num2str(Ppp(gamma_idx, 1))])
%     end

    % For each \gamma, calculate the probability of falsa alarm (FPR) and
    % the probability of detection (TPR)
    for gamma_idx = 1:length(gamma_vector)

        gamma = gamma_vector(gamma_idx);
        detected_H0 = real(GLRT_H0)>gamma;
        detected_H1 = real(GLRT_H1)>gamma;

        % For each value of gamma, compute the FPR (false positive rate) and
        % the TPR (true positive rate)
        % False positive: P(D=1 | L=0;gamma)
        Ppn(gamma_idx, 1) = sum(detected_H0)/config.N_mc;
        % True positive: P(D=1 | L=1;gamma)
        Ppp(gamma_idx, 1) = sum(detected_H1)/config.N_mc;
        % Probability of error: P(error; gamma)
        %Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
        disp(['Results have been computed for gamma=', num2str(gamma), ' with TPR=', num2str(Ppp(gamma_idx, 1))])
    end


    % Plot ROC for this configuration
    semilogx(Ppn(:,1), Ppp(:,1), 'displayname', plot_strings(idx_jamming), ...
        'linewidth', 1.5, 'color', config.plot_color(idx_jamming), ...
        'linestyle', config.plot_linestyle(idx_jamming));
    hold on,
end

legend('show', 'location', 'southeast')
grid on
title('OFDM L = 3; ISM band 2.4 GHz is attacked', 'fontsize', 20)
axis([0 1 0 1])
xlim([0.005, 1]) % same x axis limit as in paper #2
xlabel('FPR P(D=1 | H0)')
ylabel('TPR P(D=1 | H1)')

%% ADAPTIVE WAVEFORM DESIGN - REPRODUCE FIGURE 2
close all,
clear all,
clc,
config = load_config(0); % Load configuration

%%% Test for different number of subcarriers L
% Results are plotted on the bottom subfigure

% Configuration
rng('default') % Initialize seed
% Get LOS vector
LOS_vector = get_LOS_vector(config);
% r = sqrt(config.target_dist_north^2 + config.target_dist_east^2); % distance between target and radar (meters)
% LOS_vector = [config.target_dist_east, config.target_dist_north]/r;

% Create figure
figure()
subplot(2,1,1)

% Loop through different number of subcarriers L
config.L = 3;
config.TCR = -10;
config.N_dwells = 2;
config.N_mc = 5e2;
gamma_vector = config.gamma_vector; % values of \gamma tested to compute the ROC curve


% Generate OFDM measurements
%     config.TCR = TCR_vector(idx_TCR); % update number of subcarriers
%     config.Deltaf = config.B/(config.L);
%     config.T = 1/config.Deltaf;

% OFDM weights
% These weights can be updated in the Adaptive Waveform Design (AWD) module
A = eye(config.L)/sqrt(config.L);
A = eye(config.L);

% Doppler information
beta = 2*config.target_velocity'*LOS_vector/physconst('LightSpeed');
omega_dopp = 2*pi*beta*(config.fc + (0:(config.L-1))'*config.Deltaf);
doppler_mat = exp(1i*omega_dopp*(1:config.N));
%Phi_t = get_doppler_mat(config, LOS_vector); % LxN matrix containing the Doppler information through the parameter config.eta (target)
Phi_t = doppler_mat;

% Initialize metrics to compute the ROC curve
Ppp = zeros(length(gamma_vector), 1); % probability of positive detector and positive label p(D=1|H1)
Ppn = zeros(length(gamma_vector), 1); % p(D=1|H0)

% For each Monte Carlo realization, generate OFDM measurements
% according to the hypothesis indicated by the label, and calculate
% GLRT.
GLRT_H0 = zeros(config.N_mc, 1);
GLRT_H1 = zeros(config.N_mc, 1);

config.N_dwells_vector = [1,2];


for idx_Ndwells = 1:length(config.N_dwells_vector)
    config.N_dwells = config.N_dwells_vector(idx_Ndwells);
    for idx_mc = 1:config.N_mc
        % OFDM weights
        % These weights can be updated in the Adaptive Waveform Design (AWD) module
        A = eye(config.L)/sqrt(config.L);
        A = eye(config.L)*0.1;
        for idx_dwell = 1:config.N_dwells
           
            % Generate target coefficients
            X_H0 = zeros(config.L); % all zeros
            X_H1 = diag(sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1))); % LxL complex diagonal matrix repressenting the scattering coefficients (target)

            % Get Sigma_c
            sum_n = 0; % TODO: clean this loop
            for n = 1:config.N
                sum_n = sum_n + (A*X_H1*Phi_t(:,n))'*(A*X_H1*Phi_t(:,n));
            end
            tcr_lin = 10^(config.TCR/10);
            alpha = 1/config.N*sum_n/tcr_lin;
            R_sqrt = (randn(config.L) + 1i*randn(config.L))*(1/sqrt(2));
            R = R_sqrt'*R_sqrt;
            k = sqrt(alpha/trace(A*R*A'));
            Sigma_c = k^2*R;
        
            % Get clutter coefficients and doppler info
            XcPhic = get_noise_mat(config, Sigma_c);
            XcPhic(1,:) = XcPhic(1,:)*1e5;
        
            % Get noise matrix
            E = zeros(config.L,config.N);
            cnr_lin = 10^(config.CNR/10);
            var_noise = trace(A*Sigma_c*A')/config.L/cnr_lin;
            for i=1:config.N
                % We generate N vectors e(t), each one satisfying:
                % E[e(t)*e^H(t)] = Sigma_c
                e = sqrt(var_noise/2)*(randn(config.L,1)+1i*randn(config.L,1));
                E(:,i) = e;
            end
        
            % Get OFDM measurments
            Y_H0 = A*XcPhic +E;
            Y_H1 = A*X_H1*Phi_t + A*XcPhic + E;
        
            % MLE of target coefficients X under hypothesis H1
            X_hat = diag(diag(A\Y_H1*Phi_t'*inv(Phi_t*Phi_t')));
            X_hat = X_H1;

%             X_H1(1,1) = X_H1(1,1)*1e3;

            var_noise = var_noise*1e4;
            if idx_dwell<config.N_dwells
                % ***** Second dwell N/2
                config.mu = 0;
                config.E_A = 1000;
                fun = @(A_f)(-1*trace(inv(A_f*Sigma_c*A_f'+var_noise*eye(config.L))*A_f*X_H1*Phi_t*Phi_t'*X_H1'*A_f'));     
                A0 = eye(config.L);
                A_min = fminsearch(fun,A0);
                A = diag(diag(A_min));
                A = eye(config.L)/sqrt(config.L);
            end
            
        end % end dwell loop
        % GLRT
        GLRT_H0(idx_mc, 1)= det(Y_H0*Y_H0')/det((Y_H0-A*X_hat*Phi_t)*(Y_H0-A*X_hat*Phi_t)');
        GLRT_H1(idx_mc, 1)= det(Y_H1*Y_H1')/det((Y_H1-A*X_hat*Phi_t)*(Y_H1-A*X_hat*Phi_t)');
%         disp(idx_mc)
    end % end monte carlo loop
    
    % For each \gamma, calculate the probability of falsa alarm (FPR) and
    % the probability of detection (TPR)
    for gamma_idx = 1:length(gamma_vector)
    
        gamma = gamma_vector(gamma_idx);
        detected_H0 = real(GLRT_H0)>gamma;
        detected_H1 = real(GLRT_H1)>gamma;
    
        % For each value of gamma, compute the FPR (false positive rate) and
        % the TPR (true positive rate)
        % False positive: P(D=1 | L=0;gamma)
        Ppn(gamma_idx, 1) = sum(detected_H0)/config.N_mc;
        % True positive: P(D=1 | L=1;gamma)
        Ppp(gamma_idx, 1) = sum(detected_H1)/config.N_mc;
        % Probability of error: P(error; gamma)
        %Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
        disp(['Results have been computed for gamma=', num2str(gamma), ' with TPR=', num2str(Ppp(gamma_idx, 1))])
    end
    
    
    % Plot ROC for this configuration
    semilogx(Ppn(:,1), Ppp(:,1), 'displayname', [num2str(config.N_dwells),' dwells'], ...
        'linewidth', 1.5, 'color', config.plot_color(idx_Ndwells), ...
        'linestyle', config.plot_linestyle(idx_Ndwells));
    hold on,
end

legend('show', 'location', 'southeast')
grid on
title('ROC')
axis([0 1 0 1])
xlim([0.005, 1]) % same x axis limit as in paper #2
xlabel('FPR P(D=1 | H0)')
ylabel('TPR P(D=1 | H1)')

%% Initialize scenario (paper 0) - FIGURE 1
close all,
clear all,
clc,
config = load_config(0); % Load configuration

%%% Test for different number of subcarriers L
% Results are plotted on the bottom subfigure

% Configuration
rng('default') % Initialize seed
% Get LOS vector
LOS_vector = get_LOS_vector(config);
% r = sqrt(config.target_dist_north^2 + config.target_dist_east^2); % distance between target and radar (meters)
% LOS_vector = [config.target_dist_east, config.target_dist_north]/r;

% Create figure
figure()
subplot(2,1,1)

% Loop through different number of subcarriers L
config.L = 3;
TCR_vector = [-15, -10, -5]; % values of L for which a ROC curve is computed
gamma_vector = config.gamma_vector; % values of \gamma tested to compute the ROC curve
for idx_TCR = 1:length(TCR_vector)

    % Generate OFDM measurements
    config.TCR = TCR_vector(idx_TCR); % update number of subcarriers
%     config.Deltaf = config.B/(config.L);
%     config.T = 1/config.Deltaf;

    % OFDM weights
    % These weights can be updated in the Adaptive Waveform Design (AWD) module
    A = eye(config.L)/sqrt(config.L);
    A = eye(config.L);
    
    % Doppler information
    beta = 2*config.target_velocity'*LOS_vector/physconst('LightSpeed');
    omega_dopp = 2*pi*beta*(config.fc + (0:(config.L-1))'*config.Deltaf);
    doppler_mat = exp(1i*omega_dopp*(1:config.N));
    %Phi_t = get_doppler_mat(config, LOS_vector); % LxN matrix containing the Doppler information through the parameter config.eta (target)
    Phi_t = doppler_mat;

    % Initialize metrics to compute the ROC curve
    Ppp = zeros(length(gamma_vector), 1); % probability of positive detector and positive label p(D=1|H1)
    Ppn = zeros(length(gamma_vector), 1); % p(D=1|H0)
    
    % For each Monte Carlo realization, generate OFDM measurements
    % according to the hypothesis indicated by the label, and calculate
    % GLRT.
    GLRT_H0 = zeros(config.N_mc, 1);
    GLRT_H1 = zeros(config.N_mc, 1);
    for idx_mc = 1:config.N_mc
        
        % Generate target coefficients
        X_H0 = zeros(config.L); % all zeros
        X_H1 = diag(sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1))); % LxL complex diagonal matrix repressenting the scattering coefficients (target)
        
        % Get Sigma_c
        sum_n = 0; % TODO: clean this loop
        for n = 1:config.N
            sum_n = sum_n + (A*X_H1*Phi_t(:,n))'*(A*X_H1*Phi_t(:,n));
        end
        tcr_lin = 10^(config.TCR/10);
        alpha = 1/config.N*sum_n/tcr_lin;
        R_sqrt = (randn(config.L) + 1i*randn(config.L))*(1/sqrt(2));
        R = R_sqrt'*R_sqrt;
        k = sqrt(alpha/trace(A*R*A'));
        Sigma_c = k^2*R;

        % Get clutter coefficients and doppler info
        XcPhic = get_noise_mat(config, Sigma_c);

        % Get noise matrix
        E = zeros(config.L,config.N);
        cnr_lin = 10^(config.CNR/10);
        var_noise = trace(A*Sigma_c*A')/config.L/cnr_lin;
        for i=1:config.N
            % We generate N vectors e(t), each one satisfying:
            % E[e(t)*e^H(t)] = Sigma_c
            e = sqrt(var_noise/2)*(randn(config.L,1)+1i*randn(config.L,1));
            E(:,i) = e;
        end

        % Get OFDM measurments
        Y_H0 = A*XcPhic +E;
        Y_H1 = A*X_H1*Phi_t + A*XcPhic + E;

        % MLE of target coefficients X under hypothesis H1
        X_hat = diag(diag(A\Y_H1*Phi_t'*inv(Phi_t*Phi_t')));
        X_hat = X_H1;
%         X_hat = X_H1/10;
%         Pi_Phi = Phi_t'*inv(Phi_t*Phi_t')*Phi_t;
%         Pi_Phi_dag = eye(config.N) - Pi_Phi;
%         G = Y_H1*Pi_Phi_dag*Y_H1';
%         G = G - diag(G) + real(diag(G));
%         S = chol(G)*A;
%         Gamma = kron(S',chol(G)*A)
%         X_hat = inv(Gamma'*Gamma)*Gamma'*(chol(G)*Y_H1*Pi_Phi)

        % GLRT
        GLRT_H0(idx_mc, 1)= det(Y_H0*Y_H0')/det((Y_H0-A*X_hat*Phi_t)*(Y_H0-A*X_hat*Phi_t)');
        GLRT_H1(idx_mc, 1)= det(Y_H1*Y_H1')/det((Y_H1-A*X_hat*Phi_t)*(Y_H1-A*X_hat*Phi_t)');
    end
    
%     % For each \gamma, calculate the probability of falsa alarm (FPR) and
%     % the probability of detection (TPR)
%     for gamma_idx = 1:length(gamma_vector)
% 
%         gamma = gamma_vector(gamma_idx);
%         detected = real(GLRT)>gamma;
% 
%         % For each value of gamma, compute the FPR (false positive rate) and
%         % the TPR (true positive rate)
%         % Number of True/False decisions for each True/False label
%         Nd0l1 = [length(find(detected==0 & labels==1))];
%         Nd1l0 = [length(find(detected==1 & labels==0))];
%         Nd1l1 = [length(find(detected==1 & labels==1))];
%         % False negative: P(D=0 | L=1;gamma)
%         Pnp(gamma_idx, 1) = Nd0l1/Nl1;
%         % False positive: P(D=1 | L=0;gamma)
%         Ppn(gamma_idx, 1) = Nd1l0/Nl0;
%         % True positive: P(D=1 | L=1;gamma)
%         Ppp(gamma_idx, 1) = Nd1l1/Nl1;
%         % Probability of error: P(error; gamma)
%         Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
%         disp(['Results have been computed for gamma=', num2str(gamma), ' with TPR=', num2str(Ppp(gamma_idx, 1))])
%     end

    % For each \gamma, calculate the probability of falsa alarm (FPR) and
    % the probability of detection (TPR)
    for gamma_idx = 1:length(gamma_vector)

        gamma = gamma_vector(gamma_idx);
        detected_H0 = real(GLRT_H0)>gamma;
        detected_H1 = real(GLRT_H1)>gamma;

        % For each value of gamma, compute the FPR (false positive rate) and
        % the TPR (true positive rate)
        % False positive: P(D=1 | L=0;gamma)
        Ppn(gamma_idx, 1) = sum(detected_H0)/config.N_mc;
        % True positive: P(D=1 | L=1;gamma)
        Ppp(gamma_idx, 1) = sum(detected_H1)/config.N_mc;
        % Probability of error: P(error; gamma)
        %Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
        disp(['Results have been computed for gamma=', num2str(gamma), ' with TPR=', num2str(Ppp(gamma_idx, 1))])
    end


    % Plot ROC for this configuration
    semilogx(Ppn(:,1), Ppp(:,1), 'displayname', ['TCR = ', num2str(TCR_vector(idx_TCR))], ...
        'linewidth', 1.5, 'color', config.plot_color(idx_TCR), ...
        'linestyle', config.plot_linestyle(idx_TCR));
    hold on,
end

legend('show', 'location', 'southeast')
grid on
title('ROC')
axis([0 1 0 1])
xlim([0.005, 1]) % same x axis limit as in paper #2
xlabel('FPR P(D=1 | H0)')
ylabel('TPR P(D=1 | H1)')


% %% Initialize scenario (paper 0)
% % close all,
% clear all,
% clc,
config = load_config(0); % Load configuration

% Test for different number of subcarriers L
% Results are plotted on the bottom subfigure

% Configuration
rng('default') % Initialize seed
% Get LOS vector
LOS_vector = get_LOS_vector(config);
% r = sqrt(config.target_dist_north^2 + config.target_dist_east^2); % distance between target and radar (meters)
% LOS_vector = [config.target_dist_east, config.target_dist_north]/r;

% Create figure
% figure()
subplot(2,1,2)

% Loop through different number of subcarriers L
L_vector = [1, 3, 5]; % values of L for which a ROC curve is computed
gamma_vector = config.gamma_vector; % values of \gamma tested to compute the ROC curve
for idx_L = 1:length(L_vector)

    % Generate OFDM measurements
    config.L = L_vector(idx_L); % update number of subcarriers
    config.Deltaf = config.B/(config.L);
    config.T = 1/config.Deltaf;

    % OFDM weights
    % These weights can be updated in the Adaptive Waveform Design (AWD) module
    A = eye(config.L)/sqrt(config.L);
    A = eye(config.L);
    
    % Doppler information
    beta = 2*config.target_velocity'*LOS_vector/physconst('LightSpeed');
    omega_dopp = 2*pi*beta*(config.fc + (0:(config.L-1))'*config.Deltaf);
    doppler_mat = exp(1i*omega_dopp*(1:config.N));
    %Phi_t = get_doppler_mat(config, LOS_vector); % LxN matrix containing the Doppler information through the parameter config.eta (target)
    Phi_t = doppler_mat;

    % Initialize metrics to compute the ROC curve
    Ppp = zeros(length(gamma_vector), 1); % probability of positive detector and positive label p(D=1|H1)
    Ppn = zeros(length(gamma_vector), 1); % p(D=1|H0)
    
    % For each Monte Carlo realization, generate OFDM measurements
    % according to the hypothesis indicated by the label, and calculate
    % GLRT.
    GLRT_H0 = zeros(config.N_mc, 1);
    GLRT_H1 = zeros(config.N_mc, 1);
    for idx_mc = 1:config.N_mc
        
        % Generate target coefficients
        X_H0 = zeros(config.L); % all zeros
        X_H1 = diag(sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1))); % LxL complex diagonal matrix repressenting the scattering coefficients (target)
        
        % Get Sigma_c
        sum_n = 0; % TODO: clean this loop
        for n = 1:config.N
            sum_n = sum_n + (A*X_H1*Phi_t(:,n))'*(A*X_H1*Phi_t(:,n));
        end
        tcr_lin = 10^(config.TCR/10);
        alpha = 1/config.N*sum_n/tcr_lin;
        R_sqrt = (randn(config.L) + 1i*randn(config.L))*(1/sqrt(2));
        R = R_sqrt'*R_sqrt;
        k = sqrt(alpha/trace(A*R*A'));
        Sigma_c = k^2*R;

        % Get clutter coefficients and doppler info
        XcPhic = get_noise_mat(config, Sigma_c);

        % Get noise matrix
        E = zeros(config.L,config.N);
        cnr_lin = 10^(config.CNR/10);
        var_noise = trace(A*Sigma_c*A')/config.L/cnr_lin;
        for i=1:config.N
            % We generate N vectors e(t), each one satisfying:
            % E[e(t)*e^H(t)] = Sigma_c
            e = sqrt(var_noise/2)*(randn(config.L,1)+1i*randn(config.L,1));
            E(:,i) = e;
        end

        % Get OFDM measurments
        Y_H0 = A*XcPhic +E;
        Y_H1 = A*X_H1*Phi_t + A*XcPhic + E;

        % MLE of target coefficients X under hypothesis H1
        X_hat = diag(diag(A\Y_H1*Phi_t'*inv(Phi_t*Phi_t')));
        X_hat = X_H1;
%         X_hat = X_H1/10;
%         Pi_Phi = Phi_t'*inv(Phi_t*Phi_t')*Phi_t;
%         Pi_Phi_dag = eye(config.N) - Pi_Phi;
%         G = Y_H1*Pi_Phi_dag*Y_H1';
%         G = G - diag(G) + real(diag(G));
%         S = chol(G)*A;
%         Gamma = kron(S',chol(G)*A)
%         X_hat = inv(Gamma'*Gamma)*Gamma'*(chol(G)*Y_H1*Pi_Phi)

        % GLRT
        GLRT_H0(idx_mc, 1)= det(Y_H0*Y_H0')/det((Y_H0-A*X_hat*Phi_t)*(Y_H0-A*X_hat*Phi_t)');
        GLRT_H1(idx_mc, 1)= det(Y_H1*Y_H1')/det((Y_H1-A*X_hat*Phi_t)*(Y_H1-A*X_hat*Phi_t)');
    end
    
%     % For each \gamma, calculate the probability of falsa alarm (FPR) and
%     % the probability of detection (TPR)
%     for gamma_idx = 1:length(gamma_vector)
% 
%         gamma = gamma_vector(gamma_idx);
%         detected = real(GLRT)>gamma;
% 
%         % For each value of gamma, compute the FPR (false positive rate) and
%         % the TPR (true positive rate)
%         % Number of True/False decisions for each True/False label
%         Nd0l1 = [length(find(detected==0 & labels==1))];
%         Nd1l0 = [length(find(detected==1 & labels==0))];
%         Nd1l1 = [length(find(detected==1 & labels==1))];
%         % False negative: P(D=0 | L=1;gamma)
%         Pnp(gamma_idx, 1) = Nd0l1/Nl1;
%         % False positive: P(D=1 | L=0;gamma)
%         Ppn(gamma_idx, 1) = Nd1l0/Nl0;
%         % True positive: P(D=1 | L=1;gamma)
%         Ppp(gamma_idx, 1) = Nd1l1/Nl1;
%         % Probability of error: P(error; gamma)
%         Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
%         disp(['Results have been computed for gamma=', num2str(gamma), ' with TPR=', num2str(Ppp(gamma_idx, 1))])
%     end

    % For each \gamma, calculate the probability of falsa alarm (FPR) and
    % the probability of detection (TPR)
    for gamma_idx = 1:length(gamma_vector)

        gamma = gamma_vector(gamma_idx);
        detected_H0 = real(GLRT_H0)>gamma;
        detected_H1 = real(GLRT_H1)>gamma;

        % For each value of gamma, compute the FPR (false positive rate) and
        % the TPR (true positive rate)
        % False positive: P(D=1 | L=0;gamma)
        Ppn(gamma_idx, 1) = sum(detected_H0)/config.N_mc;
        % True positive: P(D=1 | L=1;gamma)
        Ppp(gamma_idx, 1) = sum(detected_H1)/config.N_mc;
        % Probability of error: P(error; gamma)
        %Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
        disp(['Results have been computed for gamma=', num2str(gamma), ' with TPR=', num2str(Ppp(gamma_idx, 1))])
    end


    % Plot ROC for this configuration
    semilogx(Ppn(:,1), Ppp(:,1), 'displayname', ['L = ', num2str(L_vector(idx_L))], ...
        'linewidth', 1.5, 'color', config.plot_color(idx_L), ...
        'linestyle', config.plot_linestyle(idx_L));
    hold on,
end

legend('show', 'location', 'southeast')
grid on
title('ROC')
axis([0 1 0 1])
xlim([0.005, 1]) % same x axis limit as in paper #2
xlabel('FPR P(D=1 | H0)')
ylabel('TPR P(D=1 | H1)')

%% Initialize scenario
% close all,
clear all,
clc,
config = load_config(1); % Load configuration

%%% Test for different number of subcarriers L
% Results are plotted on the bottom subfigure

% Configuration
rng('default') % Initialize seed
% Get LOS vector
LOS_vector = get_LOS_vector(config);

% Create figure
figure()
subplot(2,1,2)

% Loop through different number of subcarriers L
L_vector = [2, 4, 6]; % values of L for which a ROC curve is computed
gamma_vector = config.gamma_vector; % values of \gamma tested to compute the ROC curve
for idx_L = 1:length(L_vector)

    % Generate OFDM measurements
    config.L = L_vector(idx_L); % update number of subcarriers
    config.Deltaf = config.B/(config.L+1);
    config.T = 1/config.Deltaf;

    % OFDM weights
    % These weights can be updated in the Adaptive Waveform Design (AWD) module
    A = eye(config.L)/sqrt(config.L);
%     A = eye(config.L);
    
    % Doppler information
    Phi_t = get_doppler_mat(config, LOS_vector); % LxN matrix containing the Doppler information through the parameter config.eta (target)
    
    % Initialize metrics to compute the ROC curve
    Ppp = zeros(length(gamma_vector), 1); % probability of positive detector and positive label p(D=1|H1)
    Ppn = zeros(length(gamma_vector), 1); % p(D=1|H0)
    
    % For each Monte Carlo realization, generate OFDM measurements
    % according to the hypothesis indicated by the label, and calculate
    % GLRT.
    GLRT_H0 = zeros(config.N_mc, 1);
    GLRT_H1 = zeros(config.N_mc, 1);
    for idx_mc = 1:config.N_mc
        
        % Generate target coefficients
        X_H0 = zeros(config.L); % all zeros
        X_H1 = diag(sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1))); % LxL complex diagonal matrix repressenting the scattering coefficients (target)
        
        % Get Sigma_c
        sum_n = 0; % TODO: clean this loop
        for n = 1:config.N
            sum_n = sum_n + (A*X_H1*Phi_t(:,n))'*(A*X_H1*Phi_t(:,n));
        end
        snr_lin = 10^(config.SNR/10);
        alpha = 1/config.N*sum_n/snr_lin;
        R_sqrt = (randn(config.L) + 1i*randn(config.L))*(1/sqrt(2));
        R = R_sqrt'*R_sqrt;
        k = sqrt(alpha/trace(R));
        Sigma_c = k^2*R;

        % Get noise matrix E
        E = get_noise_mat(config, Sigma_c);

        % Get OFDM measurments
        Y_H0 = E;
        Y_H1 = A*X_H1*Phi_t + E;

        % MLE of target coefficients X under hypothesis H1
        X_hat = diag(diag(A\Y_H1*Phi_t'*inv(Phi_t*Phi_t')));
%         X_hat = X_H1/10;
%         Pi_Phi = Phi_t'*inv(Phi_t*Phi_t')*Phi_t;
%         Pi_Phi_dag = eye(config.N) - Pi_Phi;
%         G = Y_H1*Pi_Phi_dag*Y_H1';
%         G = G - diag(G) + real(diag(G));
%         S = chol(G)*A;
%         Gamma = kron(S',chol(G)*A)
%         X_hat = inv(Gamma'*Gamma)*Gamma'*(chol(G)*Y_H1*Pi_Phi)

        % GLRT
        GLRT_H0(idx_mc, 1)= det(Y_H0*Y_H0')/det((Y_H0-A*X_hat*Phi_t)*(Y_H0-A*X_hat*Phi_t)');
        GLRT_H1(idx_mc, 1)= det(Y_H1*Y_H1')/det((Y_H1-A*X_hat*Phi_t)*(Y_H1-A*X_hat*Phi_t)');
    end
    
%     % For each \gamma, calculate the probability of falsa alarm (FPR) and
%     % the probability of detection (TPR)
%     for gamma_idx = 1:length(gamma_vector)
% 
%         gamma = gamma_vector(gamma_idx);
%         detected = real(GLRT)>gamma;
% 
%         % For each value of gamma, compute the FPR (false positive rate) and
%         % the TPR (true positive rate)
%         % Number of True/False decisions for each True/False label
%         Nd0l1 = [length(find(detected==0 & labels==1))];
%         Nd1l0 = [length(find(detected==1 & labels==0))];
%         Nd1l1 = [length(find(detected==1 & labels==1))];
%         % False negative: P(D=0 | L=1;gamma)
%         Pnp(gamma_idx, 1) = Nd0l1/Nl1;
%         % False positive: P(D=1 | L=0;gamma)
%         Ppn(gamma_idx, 1) = Nd1l0/Nl0;
%         % True positive: P(D=1 | L=1;gamma)
%         Ppp(gamma_idx, 1) = Nd1l1/Nl1;
%         % Probability of error: P(error; gamma)
%         Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
%         disp(['Results have been computed for gamma=', num2str(gamma), ' with TPR=', num2str(Ppp(gamma_idx, 1))])
%     end

    % For each \gamma, calculate the probability of falsa alarm (FPR) and
    % the probability of detection (TPR)
    for gamma_idx = 1:length(gamma_vector)

        gamma = gamma_vector(gamma_idx);
        detected_H0 = real(GLRT_H0)>gamma;
        detected_H1 = real(GLRT_H1)>gamma;

        % For each value of gamma, compute the FPR (false positive rate) and
        % the TPR (true positive rate)
        % False positive: P(D=1 | L=0;gamma)
        Ppn(gamma_idx, 1) = sum(detected_H0)/config.N_mc;
        % True positive: P(D=1 | L=1;gamma)
        Ppp(gamma_idx, 1) = sum(detected_H1)/config.N_mc;
        % Probability of error: P(error; gamma)
        %Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
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
gamma_vector = config.gamma_vector; % values of \gamma tested to compute the ROC curve
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
