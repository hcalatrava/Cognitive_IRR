function varargout = get_ROC_from_config(config, LOS_vector, T_rangecell)   

    % Given the input configuration, obtain the OFDM model parameters
    [A_nonopt, X_H1, Phi_t, Sigma_c] = get_OFDM_model(config, LOS_vector, T_rangecell);

    % Initialize GLRT metrics
    GLRT_H0 = zeros(config.N_mc, 1);
    GLRT_H1 = zeros(config.N_mc, 1);

    % Initialize GLRT metrics for AWD, if enabled
    if config.enable_awd
        GLRT_H0_opt = zeros(config.N_mc, 1);
        GLRT_H1_opt = zeros(config.N_mc, 1);
    end

    % For each Monte Carlo realization, generate OFDM measurements
    finalA = zeros(size(A_nonopt,1), size(A_nonopt,2), config.N_mc);
    vini = round(rand(1)*10);
    for idx_mc = 1:config.N_mc

        [Y_H0, Y_H1] = build_OFDM_meas(config, Sigma_c, A_nonopt, X_H1, Phi_t);

        % If config.known_velocity is 0, we maximize the GLRT w.r.t. the
        % target velocity. If config.known velocity is 1, we assume the
        % target velocity to be known
        if config.known_velocity 
            % MLE of target coefficients X under hypothesis H1
            if config.clairvoyant
                X_hat = X_H1;
            else 
                X_hat = inv(A_nonopt) * (Y_H1*Phi_t') * pinv(Phi_t*Phi_t'); % expression after eq 8 paper 0
    %             X_hat = inv(A_nonopt)*Y_H1*Phi_t'*inv(Phi_t*Phi_t');
    %             X_hat = diag(diag(X_hat));
            end
            
            % GLRT (non-optimal)
            GLRT_H0(idx_mc, 1)= det(Y_H0*Y_H0')/det((Y_H0-A_nonopt*X_hat*Phi_t)*(Y_H0-A_nonopt*X_hat*Phi_t)');
            GLRT_H1(idx_mc, 1)= det(Y_H1*Y_H1')/det((Y_H1-A_nonopt*X_hat*Phi_t)*(Y_H1-A_nonopt*X_hat*Phi_t)');
        else
            [v, GLR_max] = fminsearch(@(v) get_OFDM_GLRT(v, config, LOS_vector, Y_H0, A_nonopt, T_rangecell), [vini;vini]);
            GLRT_H0(idx_mc, 1) = -GLR_max;
            [v, GLR_max] = fminsearch(@(v) get_OFDM_GLRT(v, config, LOS_vector, Y_H1, A_nonopt, T_rangecell), [vini;vini]);
            GLRT_H1(idx_mc, 1) = -GLR_max;
        end
    
        % Adaptive Waveform Design (AWD) module:
        if config.enable_awd
            A_ini = ones(1,config.L);
%             A_ini = diag(A_nonopt);
            sgm_sqr = 1; % set to 1 (best!) or sqrt(Sigma_c(1,1)) for good AWD performance
            [A_opt, ~] = fminsearch(@(A) opt_waveform(A, config.L, X_H1, Phi_t, chol(Sigma_c), sgm_sqr), A_ini);
            A_opt = diag(A_opt);
            finalA(:,:,idx_mc) = A_opt;  

            % GLRT (optimal)
            X_hat_opt = X_hat;
%             X_hat_opt = inv(A_opt)*Y_H1*Phi_t'*inv(Phi_t*Phi_t');
            GLRT_H0_opt(idx_mc, 1)= det(Y_H0*Y_H0')/det((Y_H0-A_opt*X_hat_opt*Phi_t)*(Y_H0-A_opt*X_hat_opt*Phi_t)');
            GLRT_H1_opt(idx_mc, 1)= det(Y_H1*Y_H1')/det((Y_H1-A_opt*X_hat_opt*Phi_t)*(Y_H1-A_opt*X_hat_opt*Phi_t)');
        end

        if mod(idx_mc, 50) == 0
            disp(['Monte carlo it ', num2str(idx_mc)])
        end
    end

%     opt_A = mean(finalA); % we average A matrices from all MC realizations

    [Ppn, Ppp] = get_ROC_from_GLRTs(config, GLRT_H0, GLRT_H1);
    varargout{1} = Ppn;
    varargout{2} = Ppp;

    if config.enable_awd
        [Ppn_opt, Ppp_opt] = get_ROC_from_GLRTs(config, GLRT_H0_opt, GLRT_H1_opt);
        varargout{3} = Ppn_opt;
        varargout{4} = Ppp_opt;
    end
   