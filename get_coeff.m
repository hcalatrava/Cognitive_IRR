function [coefficients, Sigma_c] = get_coeff(source, config, A, alpha)
%GET_COEFF returns the Lx1 complex vector representing the scattering
%coefficients of the target or clutter

if strcmp(source,'target')
    % Generate target coefficients to build X_t
    coefficients = sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1));

else % strcmp(source,'clutter')
    % Generate clutter coefficients to build X_c

    % Sigma_c^(1/2) = L is generated from CN(0,1)
    R_sqrt = (randn(config.L) + 1i*randn(config.L))*(1/sqrt(2));
    R = R_sqrt'*R_sqrt;
%     R = R/trace(R);

    % calculate Sigma_c scale factor assuming A is the identity matrix
%     A = eye(config.L);
    k = sqrt(alpha/trace(A*R*A'));
    k = sqrt(alpha/trace(R));

    L = chol(R);
    coefficients = k*L'*(randn(config.L,1) + 1i*randn(config.L,1))/sqrt(2);

    Sigma_c = k^2*R;
%     Sigma_c = ones(config.L, config.L);
    
    
end

%     % Generate Hermitian matrix
%     R = randn(config.L) + 1i*randn(config.L);
%     R = R_aux'*R_aux;
%     for i = 1:config.L
%         R(i,i) = R(i,i)/config.L*2;
%     end
%     R = zeros(config.L, config.L);
%     for i = 1:config.L
%         for j = 1:config.L
%             R(i,j) = (randn(1) + 1i*randn(1))*sqrt(1/2);
% %             R(j,i) = R(i,j);
% %             if i==j
% %                 R(i,j) = abs(real(R(i,j)));
% %             end
%         end
%     end
%     R = (R + R')/2;
%     for i = 1:config.L
%         for j = 1:config.L
%             if i==j
%                 R(i,j) = real(R(i,j));
%             end
%         end
%     end