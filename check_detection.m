function detected = check_detection(config, target_position, uav_position, velocity_uav, LOS_vector) 

    aux_north = config.target_dist_north;
    aux_east = config.target_dist_east;
    aux_velocity = config.target_velocity;

    % RECEIVED MEASUREMENTS
    % The T_rangecell changes for each UAV, because depending on the
    % UAV position the delay is different
    position_north = uav_position(1);
    position_east = uav_position(2);
    T_rangecell = sqrt(position_east^2 + position_north^2)*2/config.c;
    config.target_dist_north = position_north;
    config.target_dist_east = position_east;
    config.target_velocity = velocity_uav';
    LOS_vector = get_LOS_vector(config);
    [A_nonopt, X_H1, Phi_t, Sigma_c] = get_OFDM_model(config, LOS_vector, T_rangecell);
    [Y_H0, Y_H1] = build_OFDM_meas(config, Sigma_c, A_nonopt, X_H1, Phi_t);
    Y_H1 = A_nonopt*config.X_H1*Phi_t;

    % 'ASSUMED MEASUREMENTS' (i.e., measurements computed with assumed
    % data)
    config.target_dist_north = aux_north;
    config.target_dist_east = aux_east;
    config.target_velocity = aux_velocity;
    T_rangecell = sqrt(config.target_dist_east^2 + config.target_dist_north^2)*2/config.c; % roundtrip time
    LOS_vector = get_LOS_vector(config);
    [A_nonopt, X_H1, Phi_t, Sigma_c] = get_OFDM_model(config, LOS_vector, T_rangecell);
    [Y_H0_correct, Y_H1_correct] = build_OFDM_meas(config, Sigma_c, A_nonopt, X_H1, Phi_t);
    Y_H1_correct = A_nonopt*config.X_H1*Phi_t;

    % COMPUTE GLRT
    if config.clairvoyant
        X_hat = X_H1;
    else 
        X_hat = inv(A_nonopt) * (Y_H1*Phi_t') * pinv(Phi_t*Phi_t'); % expression after eq 8 paper 0
    %             X_hat = inv(A_nonopt)*Y_H1*Phi_t'*inv(Phi_t*Phi_t');
    %             X_hat = diag(diag(X_hat));
    end
    % GLRT (non-optimal)
    GLRT_H1 = det(Y_H1*Y_H1')/det((Y_H1-A_nonopt*X_hat*Phi_t)*(Y_H1-A_nonopt*X_hat*Phi_t)')
    if real(GLRT_H1)>3e3
        detected=1;
    else
        detected = 0;
    end
    
end
