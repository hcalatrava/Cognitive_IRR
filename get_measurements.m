function [Y, Sigma_c, var_noise, X_t] = get_measurements(config, hypothesis, A, Phi_t, Phi_c)
% ------------------------------------------------------------------------------
% Cognitive Interference Resilient Radar (Cognitive_IRR) 
% Author: Helena Calatrava
% Affiliation: Northeastern University, Boston, United Sates
% Date: July 2023
%
% GET_MEASUREMENTS returns measurements Y as in statistical model from (1)
%
% References:
% 1 - Target Detection in Clutter Using Adaptive OFDM Radar
%     Authors: Satyabrata Sen, Arye Nehorai
% 2 - Adaptive OFDM Radar for Target Detection in Multipath Scenarios
%     Authors: Satyabrata Sen, Arye Nehorai
% ------------------------------------------------------------------------------

% Complex coefficients
X_t_H1 = diag(get_coeff('target', config, A, 0)); % LxL complex diagonal matrix repressenting the scattering coefficients (target)
X_t_H0 = zeros(config.L);

% Complex coefficients
if strcmp(hypothesis, 'H1')
    X_t = X_t_H1;
else % strcmp(hypothesis, 'H0')
    X_t = X_t_H0;
end

% Get scaling factor for Sigma_c
sum_n = 0;
for n = 1:config.N
    sum_n = sum_n + (A*X_t_H1*Phi_t(:,n))'*(A*X_t_H1*Phi_t(:,n));
end
tcr_lin = 10^(config.TCR/20);
alpha = 1/config.N*sum_n/tcr_lin;

% Get clutter coefficients
[clutter_coeffs, Sigma_c] = get_coeff('clutter', config, A, alpha);

% Scale Sigma_c
%Sigma_c_scaled = scale_sigma_c(X_t, );

% Build X_c with the computed clutter coefficients
X_c = diag(clutter_coeffs); % LxL complex diagonal matrix repressenting the scattering coefficients (clutter)

% Measurements model noise
% e(t) accounts for the measurement noise and co-channel interference (CCI)
[E, var_noise] = get_noise_mat(config, Sigma_c, A); % LxN measurements model noise matrix

% Build the measurements
Y = A*X_t*Phi_t + A*X_c*Phi_c + E;
end