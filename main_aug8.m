% %% Initialize scenario (paper 0) - FIGURE 1
% close all,
% clear all,
% clc,
% config = load_config(0); % Load configuration
% format long e;
% 
% 
% %%% Test for different number of subcarriers L
% % Results are plotted on the bottom subfigure
% 
% % Configuration
% rng('default') % Initialize seed
% % Get LOS vector
% LOS_vector = get_LOS_vector(config);
% % r = sqrt(config.target_dist_north^2 + config.target_dist_east^2); % distance between target and radar (meters)
% % LOS_vector = [config.target_dist_east, config.target_dist_north]/r;
% 
% % Create figure
% figure()
% subplot(2,1,1)
% 
% % Loop through different number of subcarriers L
% config.L = 3;
% TCR_vector = [-15, -10, -5]; % values of L for which a ROC curve is computed
% gamma_vector = config.gamma_vector; % values of \gamma tested to compute the ROC curve
% for idx_TCR = 1:length(TCR_vector)
% 
%     % Generate OFDM measurements
%     config.TCR = TCR_vector(idx_TCR); % update number of subcarriers
% %     config.Deltaf = config.B/(config.L);
% %     config.T = 1/config.Deltaf;
% 
%     % OFDM weights
%     % These weights can be updated in the Adaptive Waveform Design (AWD) module
%     A = eye(config.L)/sqrt(config.L);
%     A = eye(config.L);
%     
%     % Doppler information
%     beta = 2*config.target_velocity'*LOS_vector'/config.c;
%     omega_dopp = 2*pi*beta*(config.fc + (0:(config.L-1))'*config.Deltaf);
%     doppler_mat = exp(1i*omega_dopp*(1:config.N));
%     %Phi_t = get_doppler_mat(config, LOS_vector); % LxN matrix containing the Doppler information through the parameter config.eta (target)
%     Phi_t = doppler_mat;
% 
%     % Initialize metrics to compute the ROC curve
%     Ppp = zeros(length(gamma_vector), 1); % probability of positive detector and positive label p(D=1|H1)
%     Ppn = zeros(length(gamma_vector), 1); % p(D=1|H0)
%     
%     % For each Monte Carlo realization, generate OFDM measurements
%     % according to the hypothesis indicated by the label, and calculate
%     % GLRT.
%     GLRT_H0 = zeros(config.N_mc, 1);
%     GLRT_H1 = zeros(config.N_mc, 1);
%     for idx_mc = 1:config.N_mc
%         
%         % Generate target coefficients
%         X_H0 = zeros(config.L); % all zeros
%         X_H1 = diag(sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1))); % LxL complex diagonal matrix repressenting the scattering coefficients (target)
%         
%         % Get Sigma_c
%         sum_n = 0; % TODO: clean this loop
%         for n = 1:config.N
%             sum_n = sum_n + (A*X_H1*Phi_t(:,n))'*(A*X_H1*Phi_t(:,n));
%         end
%         tcr_lin = 10^(config.TCR/10);
%         alpha = 1/config.N*sum_n/tcr_lin;
%         R_sqrt = (randn(config.L) + 1i*randn(config.L))*(1/sqrt(2));
%         R = R_sqrt'*R_sqrt;
%         k = sqrt(alpha/trace(A*R*A'));
%         Sigma_c = k^2*R;
% 
%         % Get clutter coefficients and doppler info
%         XcPhic = get_noise_mat(config, Sigma_c);
% 
%         % Get noise matrix
%         E = zeros(config.L,config.N);
%         cnr_lin = 10^(config.CNR/10);
%         var_noise = trace(A*Sigma_c*A')/config.L/cnr_lin;
%         for i=1:config.N
%             % We generate N vectors e(t), each one satisfying:
%             % E[e(t)*e^H(t)] = Sigma_c
%             e = sqrt(var_noise/2)*(randn(config.L,1)+1i*randn(config.L,1));
%             E(:,i) = e;
%         end
% 
%         % Get OFDM measurments
%         Y_H0 = A*XcPhic +E;
%         Y_H1 = A*X_H1*Phi_t + A*XcPhic + E;
% 
%         % MLE of target coefficients X under hypothesis H1
%         X_hat = diag(diag(A\Y_H1*Phi_t'*inv(Phi_t*Phi_t')));
% %         X_hat = X_H1;
% %         X_hat = X_H1/10;
% %         Pi_Phi = Phi_t'*inv(Phi_t*Phi_t')*Phi_t;
% %         Pi_Phi_dag = eye(config.N) - Pi_Phi;
% %         G = Y_H1*Pi_Phi_dag*Y_H1';
% %         G = G - diag(G) + real(diag(G));
% %         S = chol(G)*A;
% %         Gamma = kron(S',chol(G)*A)
% %         X_hat = inv(Gamma'*Gamma)*Gamma'*(chol(G)*Y_H1*Pi_Phi)
% 
%         % GLRT
%         GLRT_H0(idx_mc, 1)= det(Y_H0*Y_H0')/det((Y_H0-A*X_hat*Phi_t)*(Y_H0-A*X_hat*Phi_t)');
%         GLRT_H1(idx_mc, 1)= det(Y_H1*Y_H1')/det((Y_H1-A*X_hat*Phi_t)*(Y_H1-A*X_hat*Phi_t)');
%     end
%     
% %     % For each \gamma, calculate the probability of falsa alarm (FPR) and
% %     % the probability of detection (TPR)
% %     for gamma_idx = 1:length(gamma_vector)
% % 
% %         gamma = gamma_vector(gamma_idx);
% %         detected = real(GLRT)>gamma;
% % 
% %         % For each value of gamma, compute the FPR (false positive rate) and
% %         % the TPR (true positive rate)
% %         % Number of True/False decisions for each True/False label
% %         Nd0l1 = [length(find(detected==0 & labels==1))];
% %         Nd1l0 = [length(find(detected==1 & labels==0))];
% %         Nd1l1 = [length(find(detected==1 & labels==1))];
% %         % False negative: P(D=0 | L=1;gamma)
% %         Pnp(gamma_idx, 1) = Nd0l1/Nl1;
% %         % False positive: P(D=1 | L=0;gamma)
% %         Ppn(gamma_idx, 1) = Nd1l0/Nl0;
% %         % True positive: P(D=1 | L=1;gamma)
% %         Ppp(gamma_idx, 1) = Nd1l1/Nl1;
% %         % Probability of error: P(error; gamma)
% %         Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
% %         disp(['Results have been computed for gamma=', num2str(gamma), ' with TPR=', num2str(Ppp(gamma_idx, 1))])
% %     end
% 
%     % For each \gamma, calculate the probability of falsa alarm (FPR) and
%     % the probability of detection (TPR)
%     for gamma_idx = 1:length(gamma_vector)
% 
%         gamma = gamma_vector(gamma_idx);
%         detected_H0 = real(GLRT_H0)>gamma;
%         detected_H1 = real(GLRT_H1)>gamma;
% 
%         % For each value of gamma, compute the FPR (false positive rate) and
%         % the TPR (true positive rate)
%         % False positive: P(D=1 | L=0;gamma)
%         Ppn(gamma_idx, 1) = sum(detected_H0)/config.N_mc;
%         % True positive: P(D=1 | L=1;gamma)
%         Ppp(gamma_idx, 1) = sum(detected_H1)/config.N_mc;
%         % Probability of error: P(error; gamma)
%         %Pe(gamma_idx, 1) = Ppn(end)*p_H0 + Pnp(end)*p_H1;
%         disp(['Results have been computed for gamma=', num2str(gamma), ' with TPR=', num2str(Ppp(gamma_idx, 1))])
%     end
% 
% 
%     % Plot ROC for this configuration
%     semilogx(Ppn(:,1), Ppp(:,1), 'displayname', ['TCR = ', num2str(TCR_vector(idx_TCR))], ...
%         'linewidth', 1.5, 'color', config.plot_color(idx_TCR), ...
%         'linestyle', config.plot_linestyle(idx_TCR));
%     hold on,
% end
% 
% legend('show', 'location', 'southeast')
% grid on
% title('ROC')
% axis([0 1 0 1])
% xlim([0.005, 1]) % same x axis limit as in paper #2
% xlabel('FPR P(D=1 | H0)')
% ylabel('TPR P(D=1 | H1)')


%% Initialize scenario (paper 0)
close all,
clear all,
clc,
config = load_config(0); % Load configuration

format long e;

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

vini = round(rand(1)*10);
GLR_H0_L5 = zeros(config.N_mc, 1)
GLR_H1_L5 = zeros(config.N_mc, 1)

% Loop through different number of subcarriers L
L_vector = [2,4,6]; % values of L for which a ROC curve is computed
gamma_vector = config.gamma_vector; % values of \gamma tested to compute the ROC curve
for idx_L = 1:length(L_vector)

    % Generate OFDM measurements
    config.L = L_vector(idx_L); % update number of subcarriers
%     config.Deltaf = config.B/(config.L);
%     config.T = 1/config.Deltaf;

    % OFDM weights
    % These weights can be updated in the Adaptive Waveform Design (AWD) module
    A = eye(config.L)/sqrt(config.L);
%     A = eye(config.L);
    
    % Doppler information
    T_rangecell = sqrt(config.target_dist_east^2 + config.target_dist_north^2)*2/config.c;
%     T_rangecell = 0;
    beta = 2*config.target_velocity'*LOS_vector/config.c;
    omega_dopp = 2*pi*beta*(config.fc + (1:(config.L))'*config.Deltaf);
    doppler_mat = exp(1i*omega_dopp*((0:config.N-1)*config.T_PRI + T_rangecell)); %
%     doppler_mat = exp(1i*omega_dopp*((1:config.N)*config.T_PRI + T_rangecell));
    %Phi_t = get_doppler_mat(config, LOS_vector); % LxN matrix containing the Doppler information through the parameter config.eta (target)
    Phi_t = doppler_mat;

% %         % Doppler information
%     omega_dopp = 2*pi*beta*(config.fc + (0:(config.L-1))'*config.Deltaf);
%     doppler_mat = exp(1i*omega_dopp*(1:config.N));
% %     %Phi_t = get_doppler_mat(config, LOS_vector); % LxN matrix containing the Doppler information through the parameter config.eta (target)
% %     Phi_t = doppler_mat;
% 
% %     beta = 2*config.target_velocity'*LOS_vector'/config.c;
%     omega_dopp = 2*pi*beta*(config.fc + (0:(config.L-1))'*config.Deltaf);
%     doppler_mat = exp(1i*omega_dopp*(1:config.N));
%     %Phi_t = get_doppler_mat(config, LOS_vector); % LxN matrix containing the Doppler information through the parameter config.eta (target)
%     Phi_t = doppler_mat

    % Generate target coefficients
    X_H0 = zeros(config.L); % all zeros
    X_H1 = X_H0;
    for l = 1:config.L
        X_H1(l, l) = sqrt(config.varX/2)*(randn(1,1) + 1i*randn(1,1));
        disp(X_H1(l,l))
    end
%         X_H1 = diag(sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1))); % LxL complex diagonal matrix repressenting the scattering coefficients (target)
    
    % Scale Sigma_c
    signal_power = 0; % TODO: clean this loop
    for n = 1:config.N
        signal_power = signal_power + (A*X_H1*Phi_t(:,n))'*(A*X_H1*Phi_t(:,n));
    end
    signal_power = signal_power/config.N;
    tcr_lin = 10^(config.TCR/10);
%     R_sqrt = (randn(config.L) + 1i*randn(config.L))*(1/sqrt(2));
    R_sqrt = eye(config.L);
    noise_power = trace(R_sqrt'*R_sqrt);
    R = R_sqrt'*R_sqrt;
    SNR = signal_power/noise_power;
%         k = sqrt(alpha/trace(A*R*A'));
    scaling_factor = sqrt(SNR/tcr_lin);
    Sigma_c_sqrt = scaling_factor*R_sqrt;
    Sigma_c = Sigma_c_sqrt'*Sigma_c_sqrt;
    

    % Get clutter coefficients and doppler info
%         XcPhic = get_noise_mat(config, Sigma_c);

    % Initialize metrics to compute the ROC curve
    Ppp = zeros(length(gamma_vector), 1); % probability of positive detector and positive label p(D=1|H1)
    Ppn = zeros(length(gamma_vector), 1); % p(D=1|H0)
    

    % For each Monte Carlo realization, generate OFDM measurements
    % according to the hypothesis indicated by the label, and calculate
    % GLRT.
    GLRT_H0 = zeros(config.N_mc, 1);
    GLRT_H1 = zeros(config.N_mc, 1);
    for idx_mc = 1:config.N_mc

        E = zeros(config.L, config.N);
        for i=1:config.N
            % We generate N vectors e(t), each one satisfying:
            % E[e(t)*e^H(t)] = Sigma_c
            e = chol(Sigma_c)'*(randn(config.L,1)+1i*randn(config.L,1));
            E(:,i) = e;
%             disp((randn(config.L,1)+1i*randn(config.L,1)))
        end

        % Get OFDM measurments
        Y_H0 = E;
        Y_H1 = A*X_H1*Phi_t + E;

        % MLE of target coefficients X under hypothesis H1
%         X_hat = diag(diag(A\Y_H1*Phi_t'*inv(Phi_t*Phi_t')));
        X_hat = X_H1;
        X_hat = inv(A) * (Y_H1*Phi_t') * pinv(Phi_t*Phi_t');
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

        % From colleague code
        bigR = Y_H0;
        [v, GLR_max] = fminsearch(@(v) GLR_function(v, config.N, config.L, config.c, config.fc, config.Deltaf, LOS_vector, LOS_vector, LOS_vector, config.T_PRI, T_rangecell, bigR, A, eye(config.L)), [vini;vini]);
        GLR_H0_L5(idx_mc,1) = -GLR_max;
        bigR = Y_H1;
        [v, GLR_max] = fminsearch(@(v) GLR_function(v, config.N, config.L, config.c, config.fc, config.Deltaf, LOS_vector, LOS_vector, LOS_vector, config.T_PRI, T_rangecell, bigR, A, eye(config.L)), [vini;vini]);
        GLR_H1_L5(idx_mc,1) = -GLR_max;
        disp(idx_mc)
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
%     semilogx(Ppn(:,1), Ppp(:,1), 'displayname', ['L = ', num2str(L_vector(idx_L))], ...
%         'linewidth', 1.5, 'color', config.plot_color(idx_L), ...
%         'linestyle', config.plot_linestyle(idx_L));
    semilogx(Ppn(:,1), Ppp(:,1), 'displayname', ['L = ', num2str(L_vector(idx_L))], ...
    'linewidth', 1.5);
    hold on,
end

legend('show', 'location', 'southeast')
grid on
title('ROC')
axis([0 1 0 1])
xlim([0.005, 1]) % same x axis limit as in paper #2
xlabel('FPR P(D=1 | H0)')
ylabel('TPR P(D=1 | H1)')



%% Plot their figure
sorted_GLR_H0_L1 = sort(GLR_H0_L1, 1, 'ascend');
sorted_GLR_H1_L1 = sort(GLR_H1_L1, 1, 'ascend');
clear GLR_H0_L1
clear GLR_H1_L1

sorted_GLR_H0_L3 = sort(GLR_H0_L3, 1, 'ascend');
sorted_GLR_H1_L3 = sort(GLR_H1_L3, 1, 'ascend');
clear GLR_H0_L3
clear GLR_H1_L3

sorted_GLR_H0_L5 = sort(GLR_H0_L5, 1, 'ascend');
sorted_GLR_H1_L5 = sort(GLR_H1_L5, 1, 'ascend');
clear GLR_H0_L5
clear GLR_H1_L5

pFA = [];
pD_L1 = [];
pD_L3 = [];
pD_L5 = [];
for i = 1 : 15 : max_MC - 15

    pFA = [pFA, (max_MC - i)/max_MC];
    pD_L1 = [pD_L1, sum(sorted_GLR_H1_L1>=sorted_GLR_H0_L1(i,1))/max_MC];
    pD_L3 = [pD_L3, sum(sorted_GLR_H1_L3>=sorted_GLR_H0_L3(i,1))/max_MC];
    pD_L5 = [pD_L5, sum(sorted_GLR_H1_L5>=sorted_GLR_H0_L5(i,1))/max_MC];

end
clear sorted_GLR_H0_L1;
clear sorted_GLR_H1_L1;
clear sorted_GLR_H0_L3;
clear sorted_GLR_H1_L3;
clear sorted_GLR_H0_L5;
clear sorted_GLR_H1_L5;

figure;
plot(pFA, pD_L1, 'b'); hold on;
plot(pFA, pD_L3, 'r'); hold on;
plot(pFA, pD_L5, 'k');