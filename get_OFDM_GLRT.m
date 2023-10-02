function GLRT_max = get_OFDM_GLRT(v, config, LOS_vector, Y, A_nonopt, T_rangecell)

    % ESTIMATION OF [Phi_t matrix]
    % Doppler information 
    % exponential of 2\pi\f_l\beta (nT_{PRI} + \tau_0)
    beta = 2*v'*LOS_vector/config.c; % paper 1, paragraph after equation (5)
    omega_dopp = 2*pi*beta*(config.fc + (1:config.L)'*config.Deltaf);
    t = ((0:config.N-1)*config.T_PRI + T_rangecell);
    doppler_mat = exp(1i*omega_dopp*t); %
    Phi_t = doppler_mat;

    X_hat = inv(A_nonopt) * (Y*Phi_t') * pinv(Phi_t*Phi_t');

    GLRT_nu = det(Y*Y');
    GLRT_de = det((Y-A_nonopt*X_hat*Phi_t)*(Y-A_nonopt*X_hat*Phi_t)');
    GLRT_max = -real(GLRT_nu)/real(GLRT_de);

end