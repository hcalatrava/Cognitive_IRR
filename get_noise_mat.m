function noise_mat = get_noise_mat(config, Sigma_c)
% ------------------------------------------------------------------------------
% Cognitive Interference Resilient Radar (Cognitive_IRR) 
% Author: Helena Calatrava
% Affiliation: Northeastern University, Boston, United Sates
% Date: July 2023
%
% GET_NOISE_MAT returns noise matrix E from the statistical model in (2)
%
% References:
% 1 - Target Detection in Clutter Using Adaptive OFDM Radar
%     Authors: Satyabrata Sen, Arye Nehorai
% 2 - Adaptive OFDM Radar for Target Detection in Multipath Scenarios
%     Authors: Satyabrata Sen, Arye Nehorai
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Statistical Model in (2)
% ------------------------------------------------------------------------------
% Option 2 (meeting June 18th)
E = zeros(config.L, config.N);
for i=1:config.N
    % We generate N vectors e(t), each one satisfying:
    % E[e(t)*e^H(t)] = Sigma_c
    e = chol(Sigma_c)'*sqrt(1/2)*(randn(config.L,1)+1i*randn(config.L,1));
    E(:,i) = e;
end
noise_mat = E;

% % Option 1 (deprecated)
% Gamma_c = 1/2*[real(Sigma_c), -imag(Sigma_c); imag(Sigma_c), real(Sigma_c)];
% % From Equation (13)
% E = zeros(config.L, config.N);
% mu = zeros(config.L*2,1);
% for i=1:config.N
%     % We generate N vectors e(t), each one satisfying:
%     % E[e(t)*e^H(t)] = Sigma_c
%     e_real_imag = mvnrnd(mu,Gamma_c);
%     E(:,i) = e_real_imag(1:config.L) + 1i*e_real_imag(config.L+1:end);
% end

end

% % [DEPRECATED]:
% % ------------------------------------------------------------------------------
% % Statistical Model in (1), Equation (13)
% % Noise matrix is obtained from CNR (clutter-to-noise ratio)
% % ------------------------------------------------------------------------------
% cnr_lin = 10^(config.CNR/20);
% var_noise = trace(A*Sigma_c*A')/config.L/cnr_lin;
% noise_mat = sqrt(var_noise/2)*(randn(config.L,config.N) + 1i*randn(config.L,config.N));