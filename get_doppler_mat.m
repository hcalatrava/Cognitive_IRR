function doppler_mat = get_doppler_mat(config)
%GET_DOPPLER_MAT returns the LxN matrix containing the Doppler information 
% of the target through the parameter config.eta
beta = config.target_velocity'*config.target_angle/physconst('LightSpeed');
omega_dopp = 2*pi*beta*(config.fc + (0:(config.L-1))'*config.Deltaf);
doppler_mat = exp(1i*omega_dopp*(0:config.N-1));
end

