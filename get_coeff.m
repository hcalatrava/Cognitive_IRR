function [coefficients, Sigma_c] = get_coeff(source, config)
%GET_COEFF returns the Lx1 complex vector representing the scattering
%coefficients of the target or clutter

if strcmp(source,'target')
    % Generate target coefficients to build X_t
    coefficients = sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1));

else % strcmp(source,'clutter')
    % Generate clutter coefficients to build X_c

    % Step 1: Generate covariance matrix Sigma_c
    Sigma_c = sqrt(config.varSigma_c/2)*(randn(config.L, config.L) + 1i*randn(config.L, config.L));

    % Step 2: Calculate covariance matrix Gamma_c
    Gamma_c = sqrt(1/2)*[real(Sigma_c) -imag(Sigma_c); imag(Sigma_c) real(Sigma_c)]; % the sqrt() is included in this calculation so that the variance remains the same after computing A'*A, which is necessaryh so that the matrix is SDP
    % We need to make Gamma_c symmetric
%     for i = 1:(config.L*2)
%         for j = 1:(config.L*2)
%             if i ~= j
%                 Gamma_c(i,j) = Gamma_c(j,i);
%             end
%         end
%     end
    
    Gamma_c = Gamma_c'*Gamma_c;
%     try chol(Gamma_c)
%     disp('Matrix is symmetric positive definite.')
%     catch ME
%     disp('Matrix is not symmetric positive definite')
%     end

    % Step 3: Generate vector [X Y], where X contains the real part of the
    % coefficients and Y contains the imaginary part of the coefficients
    Z = mvnrnd(zeros(config.L*2,1), Gamma_c);

    % Step 4: Build the coefficients vector
    coefficients = Z(1:config.L) + 1i*Z(config.L+1:end);
end