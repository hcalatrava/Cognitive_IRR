function [coefficients] = get_coeff(source, config)
%GET_COEFF returns the Lx1 complex vector representing the scattering
%coefficients of the target or clutter
if strcmp(source,'target')
    coefficients = sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1));
else % strcmp(source,'clutter')
    coefficients = sqrt(config.varSigma_c/2)*(randn(config.L,1) + 1i*randn(config.L,1));
end