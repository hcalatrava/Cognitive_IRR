function Sigma_c = get_Sigma_c(config, Phi_t, A, X_H1)
% [Sigma_c matrix]
    % Scale Sigma_c
    signal_power = 0; % TODO: clean this loop
    for n = 1:config.N
        signal_power = signal_power + (A*X_H1*Phi_t(:,n))'*(A*X_H1*Phi_t(:,n));
    end
    signal_power = signal_power/config.N;
    SNR_predefined_lin = 10^(config.SNR_predefined/10);
    if config.sigma_random == 0
        R_sqrt = eye(config.L);
    else
        R_sqrt = (randn(config.L) + 1i*randn(config.L))*(1/sqrt(2));
    end
    noise_power = trace(R_sqrt'*R_sqrt);
    R = R_sqrt'*R_sqrt;
    SNR = signal_power/noise_power;
%         k = sqrt(alpha/trace(A*R*A'));
    scaling_factor = sqrt(SNR/SNR_predefined_lin);
    Sigma_c_sqrt = scaling_factor*R_sqrt;
    Sigma_c = Sigma_c_sqrt'*Sigma_c_sqrt;
    