function varargout = get_ROC_from_config(config, LOS_vector, T_rangecell)   

    % Given the input configuration, obtain the OFDM model parameters
    [A_nonopt, X_H1, Phi_t, Sigma_c] = get_OFDM_model(config, LOS_vector, T_rangecell);

    if config.coefficients_random == 0
        X_H1 = X_H1;
    end
        
    if config.sigma_random == 0
        Sigma_c = get_Sigma_c(config, Phi_t, A_nonopt, X_H1);
    end

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
    for idx_mc = 1:config.N_mc
        if config.coefficients_random
            X_H1 = get_X_H1(config);
        end
        if config.sigma_random
            Sigma_c = get_Sigma_c(config, Phi_t, A_nonopt, X_H1);
        end

        [Y_H0, Y_H1] = build_OFDM_meas(config, Sigma_c, A_nonopt, X_H1, Phi_t);

        % If config.known_velocity is 0, we maximize the GLRT w.r.t. the
        % target velocity. If config.known velocity is 1, we assume the
        % target velocity to be known
        if config.known_velocity 
            % MLE of target coefficients X under hypothesis H1
            if config.clairvoyant
                X_hat_H1 = X_H1;
                X_hat_H0 = X_hat_H1;
            else 
                X_hat_H1 = inv(A_nonopt) * (Y_H1*Phi_t') * pinv(Phi_t*Phi_t'); % expression after eq 8 paper 0
                X_hat_H0 = inv(A_nonopt) * (Y_H0*Phi_t') * pinv(Phi_t*Phi_t'); % expression after eq 8 paper 0
%                 X_hat_H0 = X_hat_H1;
                %             X_hat = inv(A_nonopt)*Y_H1*Phi_t'*inv(Phi_t*Phi_t');
    %             X_hat = diag(diag(X_hat));
            end
            
            % GLRT (non-optimal)
            GLRT_H0(idx_mc, 1)= real(det(Y_H0*Y_H0'))/real(det((Y_H0-A_nonopt*X_hat_H0*Phi_t)*(Y_H0-A_nonopt*X_hat_H0*Phi_t)'));
            GLRT_H1(idx_mc, 1)= real(det(Y_H1*Y_H1'))/real(det((Y_H1-A_nonopt*X_hat_H1*Phi_t)*(Y_H1-A_nonopt*X_hat_H1*Phi_t)'));
        else
            if config.debug_GLRT
                disp('debug')
                GLRT_evaluation_module(config, LOS_vector, Y_H1, X_H1, A_nonopt, T_rangecell);
            end
            options = optimset('Display','iter','PlotFcns',@optimplotfval);
%             options = optimset('MaxIter', 30000, 'MaxFunEvals', 1000000);
            options = [];
            n_GLRT = config.n_GLRT;
            GLR_max_H1 = zeros(n_GLRT,1);
            GLR_max_H0 = zeros(n_GLRT,1);
            for GLR_iter = 1:n_GLRT
                vini = round(rand(1)*10);
    %     vini = 10;
                [v, GLR_max] = fminsearch(@(v) get_OFDM_GLRT(v, config, LOS_vector, Y_H1, X_H1, A_nonopt, T_rangecell), [vini;vini], options);
    %             [v, diff] = fminsearch(@(v) get_v(v, config, LOS_vector, Y_H1, X_H1, A_nonopt, T_rangecell), [vini;vini], options);
                v_H1 = v;
                GLR_max_H1(GLR_iter) = -GLR_max;     
                [v, GLR_max] = fminsearch(@(v) get_OFDM_GLRT(v, config, LOS_vector, Y_H0, X_H1, A_nonopt, T_rangecell), [vini;vini], options);            
                v_H0 = v;
                GLR_max_H0(GLR_iter) = -GLR_max;
                if config.debug_estimated_v
                    disp('DEBUG ESTIMATED VELOCITY')
                    disp('hypothesis 1')
                    disp(v_H1)
                    disp('hypothesis 0')
                    disp(v_H0)
                end
            end
            GLRT_H1(idx_mc, 1) = mean(GLR_max_H1);
            GLRT_H0(idx_mc, 1) = mean(GLR_max_H0);
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

        if mod(idx_mc, 200) == 0
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

    if config.betatest_ROC_calculation
        sorted_GLR_H0 = sort(GLRT_H0, 1, 'ascend');
        sorted_GLR_H1 = sort(GLRT_H1, 1, 'ascend');
        pFA = [];
        pD = [];
        for i = 1 : 15 : config.N_mc - 15
            pFA = [pFA, (config.N_mc - i)/config.N_mc];
            pD = [pD, sum(sorted_GLR_H1>=sorted_GLR_H0(i,1))/config.N_mc];
        end
        varargout{1} = pFA;
        varargout{2} = pD; 
    end
   