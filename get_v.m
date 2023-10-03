function diff = get_v(v, config, LOS_vector, Y, X, A_nonopt, T_rangecell)
    
%     v = config.target_velocity
    % ESTIMATION OF [Phi_t matrix]
    % Doppler information 
    % exponential of 2\pi\f_l\beta (nT_{PRI} + \tau_0)
    beta = 2*v'*LOS_vector/config.c; % paper 1, paragraph after equation (5)
    omega_dopp = 2*pi*beta*(config.fc + (1:config.L)'*config.Deltaf);
    t = ((0:config.N-1)*config.T_PRI + T_rangecell);
    doppler_mat = exp(1i*omega_dopp*t); %
    Phi_t = doppler_mat;

    % Trying to correct Phi_t matrix
%     v = v_aux
%     beta = 2*v'*LOS_vector/config.c; % paper 1, paragraph after equation (5)
%     omega_dopp = 2*pi*beta*(config.fc + (1:config.L)'*config.Deltaf)
%     t = ((0:config.N-1)*config.T_PRI);
%     omega_delay = 2*pi*(config.fc + (1:config.L)'*config.Deltaf)*T_rangecell*ones(1,config.N);
%     doppler_mat = exp(1i*(omega_dopp*t - omega_delay)); %
%     Phi_t = doppler_mat;
% 
%     v = config.target_velocity
%     beta = 2*v'*LOS_vector/config.c; % paper 1, paragraph after equation (5)
%     omega_dopp = 2*pi*beta*(config.fc + (1:config.L)'*config.Deltaf)
%     t = ((0:config.N-1)*config.T_PRI);
%     omega_delay = 2*pi*(config.fc + (1:config.L)'*config.Deltaf)*T_rangecell*ones(1,config.N);
%     doppler_mat = exp(1i*(omega_dopp*t - omega_delay)); %
%     Phi_t = doppler_mat;

    if config.clairvoyant
        X_hat = X;
    else
        X_hat = inv(A_nonopt) * (Y*Phi_t') * pinv(Phi_t*Phi_t');
    end

%     GLRT_nu = (1/config.N)*Y*Y';
%     GLRT_de = (1/config.N)*(Y-A_nonopt*X_hat*Phi_t)*(Y-A_nonopt*X_hat*Phi_t)';
%     GLRT_max = -real(det(GLRT_nu))/real(det(GLRT_de));

    diff = abs(real(mean(mean((Y-A_nonopt*X_hat*Phi_t)))));
%     disp('FMINSEARCH')
%     disp('velocity:')
%     disp(v)
%     disp('function result:')
%     disp(GLRT_max)

end