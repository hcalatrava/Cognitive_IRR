function [Y_H0, Y_H1] = build_OFDM_meas(config, Sigma_c, A, X_H1, Phi_t)

    % Build noie mat E
    E = zeros(config.L, config.N);
    
    for i=1:config.N
        % We generate N vectors e(t), each one satisfying:
        % E[e(t)*e^H(t)] = Sigma_c
        e = chol(Sigma_c)'*(randn(config.L,1)+1i*randn(config.L,1));
        E(:,i) = e;
    end
    
    % Build OFDM meas
    Y_H0 = E;
    Y_H1 = A*X_H1*Phi_t + E;