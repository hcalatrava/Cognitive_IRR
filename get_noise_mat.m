function [noise_mat] = get_noise_mat(config)
%GET_NOISE_MAT returns the measurements model noise matrix

%Sigma_E = kron(eye(config.N), config.sigma_e^2*eye(config.L));
noise_mat = sqrt(config.var_noise/2)*(randn(config.L,config.N) + 1i*randn(config.L,config.N));
end