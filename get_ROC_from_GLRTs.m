function [Ppn, Ppp] = get_ROC_from_GLRTs(config, GLRT_H0, GLRT_H1)       

    % Initialize metrics to compute the ROC curve
    Ppp = zeros(length(config.gamma_vector), 1); % probability of positive detector and positive label p(D=1|H1)
    Ppn = zeros(length(config.gamma_vector), 1); % p(D=1|H0)
        
    % For each \gamma, calculate the probability of falsa alarm (FPR) and
    % the probability of detection (TPR)
    for gamma_idx = 1:length(config.gamma_vector)
    
        gamma = config.gamma_vector(gamma_idx);
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