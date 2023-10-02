function [A, X_H1, Phi_t_int, Sigma_c, Phi_t_true] = get_OFDM_model(config, LOS_vector, T_rangecell)

    % OFDM weights
    % These weights can be updated in the Adaptive Waveform Design (AWD) module
    A = eye(config.L)/sqrt(config.L); % set like this for good AWD performance
%     A = eye(config.L);

    % Generate target coefficients
    X_H1 = diag(sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1))); % LxL complex diagonal matrix repressenting the scattering coefficients (target)

    % ------------------------------------------------------
    % -------- INTERFERENCE MODULE: DISTURB REL. DOPP. SHIFT 
    % ------------------------------------------------------
    %     beta = 2*config.target_velocity'*LOS_vector/config.c; % paper 1, paragraph after equation (5)
    % Doppler information (as in colleague code)
    % exponential of 2\pi\f_l\beta (nT_{PRI} + \tau_0)
    original_vel = config.target_velocity;
    original_vel_x = config.target_velocity(1)
    beta = 2*config.target_velocity'*LOS_vector/config.c;
    beta_n_vector = [];
    Phi_t = zeros(config.L, config.N);
    target_velocity_vector = []
    for n = 1:config.N
        t = (n-1)*config.T_PRI + T_rangecell;
        config.target_velocity = original_vel + original_vel_x*[unifrnd(-100,100);unifrnd(-100,100)]
        LOS_vector = get_LOS_vector(config);
        beta_n = 2*config.target_velocity'*LOS_vector/config.c;
       % it is good to add the if here so that this random process is not
       % done if wanting to sanity check non interfence case
        if config.add_doppler_jamming
            adversarial_noise = unifrnd(-1000,1000)*0.5;
        else
            adversarial_noise = 0;
        end
        target_velocity_vector(n) = config.target_velocity(1);
        beta_n = beta * (1 + adversarial_noise);
        beta_n_vector(n) = beta_n;
        for l = 1:config.L
            fl = config.fc + l*config.Deltaf;
            omega_dopp = 2*pi*beta_n*fl;
            Phi_t(l,n) = exp(1i*omega_dopp*t);
        end
    end
    Phi_t_int = Phi_t;
    config.target_velocity = original_vel;

    %     % Plot adversarial effect on rel. doppler shift
%     figure()
%     scatter(0:(config.N-1), target_velocity_vector, 'ro')
%     hold on,
%     yline(original_vel_x,'-','v_x','LineWidth',3, 'fontsize', 25)
%     grid on
%     title('Adversarial effect on relative doppler shift v_x', 'fontsize', 15)
% 

% %     beta = 2*config.target_velocity'*LOS_vector/config.c; % paper 1, paragraph after equation (5)
%     % Doppler information (as in colleague code)
%     % exponential of 2\pi\f_l\beta (nT_{PRI} + \tau_0)
%     beta = 2*config.target_velocity'*LOS_vector/config.c;
%     beta_n_vector = [];
%     Phi_t = zeros(config.L, config.N);
%     for n = 1:config.N
%         t = (n-1)*config.T_PRI + T_rangecell;
%        % it is good to add the if here so that this random process is not
%        % done if wanting to sanity check non interfence case
%         if config.add_doppler_jamming
%             adversarial_noise = unifrnd(-1000,1000)*0.5;
%         else
%             adversarial_noise = 0;
%         end
%         beta_n = beta * (1 + adversarial_noise);
%         beta_n_vector(n) = beta_n;
%         for l = 1:config.L
%             fl = config.fc + l*config.Deltaf;
%             omega_dopp = 2*pi*beta_n*fl;
%             Phi_t(l,n) = exp(1i*omega_dopp*t);
%         end
%     end
%     Phi_t_int = Phi_t;

%     % Plot adversarial effect on rel. doppler shift
%     figure()
%     scatter(0:(config.N-1), beta_n_vector, 'ro')
%     hold on,
%     yline(beta,'-','\beta_T','LineWidth',3, 'fontsize', 25)
%     yline(beta*(3/2))
%     yline(beta*0.5)
%     grid on
%     title('Adversarial effect on relative doppler shift \beta', 'fontsize', 15)


%     omega_dopp = 2*pi*beta*(config.fc + (1:(config.L))'*config.Deltaf);
%     t = ((0:config.N-1)*config.T_PRI + T_rangecell);
%     doppler_mat = exp(1i*omega_dopp*t); %
%     Phi_t = doppler_mat;
    
    % Doppler information (as in colleague code)
    % exponential of 2\pi\f_l\beta (nT_{PRI} + \tau_0)
    beta = 2*config.target_velocity'*LOS_vector/config.c; % paper 1, paragraph after equation (5)
    omega_dopp = 2*pi*beta*(config.fc + (1:(config.L))'*config.Deltaf);
    t = ((0:config.N-1)*config.T_PRI + T_rangecell);
    doppler_mat = exp(1i*omega_dopp*t); %
    Phi_t_true = doppler_mat;
    
    % Scale Sigma_c
    signal_power = 0; % TODO: clean this loop
    for n = 1:config.N
        signal_power = signal_power + (A*X_H1*Phi_t(:,n))'*(A*X_H1*Phi_t(:,n));
    end
    signal_power = signal_power/config.N;
    SNR_predefined_lin = 10^(config.SNR_predefined/10);
    R_sqrt = (randn(config.L) + 1i*randn(config.L))*(1/sqrt(2));
    R_sqrt = eye(config.L);
    noise_power = trace(R_sqrt'*R_sqrt);
    R = R_sqrt'*R_sqrt;
    SNR = signal_power/noise_power;
%         k = sqrt(alpha/trace(A*R*A'));
    scaling_factor = sqrt(SNR/SNR_predefined_lin);
    Sigma_c_sqrt = scaling_factor*R_sqrt;
    Sigma_c = Sigma_c_sqrt'*Sigma_c_sqrt;

    % Get clutter coefficients and doppler info
%         XcPhic = get_noise_mat(config, Sigma_c);