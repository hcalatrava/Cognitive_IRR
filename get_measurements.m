function Y = get_measurements(config, hypothesis, A, Phi_t, Phi_c)
%GET_MEASUREMENTS returns measurements Y for the specified hypothesis

% Complex coefficients
if strcmp(hypothesis, 'H1')
    X_t = diag(get_coeff('target', config)); % LxL complex diagonal matrix repressenting the scattering coefficients (target)
else % strcmp(hypothesis, 'H0')
    X_t = zeros(config.L);
end
[clutter_coeffs, Sigma_c] = get_coeff('clutter', config);
X_c = diag(clutter_coeffs); % LxL complex diagonal matrix repressenting the scattering coefficients (clutter)

% Measurements model noise
% e(t) accounts for the measurement noise and co-channel interference (CCI)
E = get_noise_mat(config, Sigma_c, A); % LxN measurements model noise matrix

% Build the measurements
Y = A*X_t*Phi_t + A*X_c*Phi_c + E;
end