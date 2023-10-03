function [A, X_H1, Phi_t, Sigma_c] = get_OFDM_model(config, LOS_vector, T_rangecell)

    % [A matrix]
    % OFDM weights
    % These weights may be updated in the Adaptive Waveform Design (AWD) module
    A = eye(config.L)/sqrt(config.L); % set like this for good AWD performance
%     A = eye(config.L);
    
    % [Phi_t matrix]
    % Doppler information 
    % exponential of 2\pi\f_l\beta (nT_{PRI} + \tau_0)
    beta = 2*config.target_velocity'*LOS_vector/config.c; % paper 1, paragraph after equation (5)
    omega_dopp = 2*pi*beta*(config.fc + (1:config.L)'*config.Deltaf);
    t = ((0:config.N-1)*config.T_PRI + T_rangecell);
    doppler_mat = exp(1i*omega_dopp*t); %
    Phi_t = doppler_mat;

    % Trying to correct Phi_t matrix
%     beta = 2*config.target_velocity'*LOS_vector/config.c; % paper 1, paragraph after equation (5)
%     omega_dopp = 2*pi*beta*(config.fc + (1:config.L)'*config.Deltaf);
%     t = ((0:config.N-1)*config.T_PRI);
%     omega_delay = 2*pi*(config.fc + (1:config.L)'*config.Deltaf)*T_rangecell*ones(1,config.N);
%     doppler_mat = exp(1i*(omega_dopp*t - omega_delay)); %
%     Phi_t = doppler_mat;

    % [X_H1 (target coefficients)]
    % Generate target coefficients
    X_H1 = diag(sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1))); % LxL complex diagonal matrix repressenting the scattering coefficients (target)
    
    % [Sigma_c matrix]
    % Scale Sigma_c
    signal_power = 0; % TODO: clean this loop
    for n = 1:config.N
        signal_power = signal_power + (A*X_H1*Phi_t(:,n))'*(A*X_H1*Phi_t(:,n));
    end
    signal_power = signal_power/config.N;
    SNR_predefined_lin = 10^(config.SNR_predefined/10);
%     R_sqrt = (randn(config.L) + 1i*randn(config.L))*(1/sqrt(2));
    R_sqrt = eye(config.L);
    noise_power = trace(R_sqrt'*R_sqrt);
    R = R_sqrt'*R_sqrt;
    SNR = signal_power/noise_power;
%         k = sqrt(alpha/trace(A*R*A'));
    scaling_factor = sqrt(SNR/SNR_predefined_lin);
    Sigma_c_sqrt = scaling_factor*R_sqrt;
    Sigma_c = Sigma_c_sqrt'*Sigma_c_sqrt;
    
    % [X_H (clutter coefficients)]
    % Get clutter coefficients and doppler info
%         XcPhic = get_noise_mat(config, Sigma_c);