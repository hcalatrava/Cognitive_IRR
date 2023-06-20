% ------------------------------------------------------------------------------
% Test for debugging purposes wrt function GET_COEFF 
% Helena Calatrava
%
% June 2023
% ------------------------------------------------------------------------------

load_config;
n_mc = 20e4;
coeff = zeros(config.L, n_mc);
R = zeros(config.L, config.L, n_mc);
R_aux = zeros(config.L, config.L, n_mc);
coeff_covar = zeros(config.L, config.L, n_mc);
for i_mc = 1:n_mc
    [coefficients, Sigma_c, R_aux_o] = get_coeff('clutter', config);
    coeff(:,i_mc) = coefficients;
    R(:,:,i_mc) = Sigma_c;
    R_aux(:,:,i_mc) = R_aux_o;
    coeff_covar(:,:,i_mc) = coefficients*coefficients';
end

% We test that Aij is zero man with variance sigma_A^2 for all i,j (WORKS)
var(real(R_aux(1,1,:))) % variance real part = variance imag part = total variance divided by 2
mean(real(R_aux),3)
mean(real(R_aux(1,1,:)))

% We test whether Rij is zero mean with variance sigma_R^2 for all i,j (DOES NOT
% WORK)
mean(real(R),3)
mean(real(R(1,1,:)))

% We check whether mean of coefficients is 0 (WORKS)
mean(coeff,2)
% and whether covariance of coefficients is R (WORKS with condition:)
mean(coeff_covar,3) % this only works when R is constant for all realizations
% it has been checked and it works, so coefficients are actually
% distributed according to covariance R

