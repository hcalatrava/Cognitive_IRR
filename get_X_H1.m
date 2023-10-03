function X_H1 = get_X_H1(config)
    % [X_H1 (target coefficients)]
    % Generate target coefficients
    X_H1 = diag(sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1))); % LxL complex diagonal matrix repressenting the scattering coefficients (target)
    