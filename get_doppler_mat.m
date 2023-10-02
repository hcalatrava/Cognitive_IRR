function doppler_mat = get_doppler_mat(config, LOS_vector)
%GET_DOPPLER_MAT returns the LxN matrix containing the Doppler information 
% of the target through the parameter config.eta
config.T_PRI = 20e-6;
r = sqrt(config.target_dist_north^2 + config.target_dist_east^2); % distance between target and radar (meters)
LOS_vector2 = [config.target_dist_east, config.target_dist_north]/r;
% LOS_vector2=LOS_vector;
beta = 2*config.target_velocity'*LOS_vector2'/physconst('LightSpeed');
tau0 = 2*r/physconst('LightSpeed'); % roundtrip delay between target and radar (seconds)
omega_dopp = 2*pi*(config.fc + (0:config.L-1)'*config.Deltaf)*(beta*config.T_PRI*(0:config.N-1)-tau0*ones(1,config.N));
omega_dopp = 2*pi*(config.fc + (0:config.L-1)'*config.Deltaf)*(beta*config.T_PRI*(0:config.N-1)-tau0*ones(1,config.N)) + 2*pi*config.Deltaf*(0:config.L-1)'*(0:config.N-1);
doppler_mat = exp(1i*omega_dopp);
inv(doppler_mat*doppler_mat')

% beta = 2*config.target_velocity'*LOS_vector'/physconst('LightSpeed');
% omega_dopp = 2*pi*beta*(config.fc + (0:(config.L-1))'*config.Deltaf);
% doppler_mat = exp(1i*omega_dopp*(1:config.N));
% 
% % Following statistical model in (2)
% r = sqrt(config.target_dist_north^2 + config.target_dist_east^2); % distance between target and radar (meters)
% tau0 = 2*r/physconst('LightSpeed'); % roundtrip delay between target and radar (seconds)
% fl = (config.fc + (0:(config.L-1))'*config.Deltaf); % subcarrier center frequency
% omega_dopp = 2*pi*fl*(beta*config.T_PRI*(0:config.N-1) - tau0); % eq. (9) from paper (2)
% doppler_mat_1 = exp(1i*omega_dopp);
% 
% % Second option for line of sight vector
% %tau0=0;
% config.T_PRI = 0.5;
% LOS_vector_2 = [config.target_dist_east, config.target_dist_north]/r;
% beta_2 = 2*config.target_velocity'*LOS_vector_2'/physconst('LightSpeed');
% omega_dopp = 2*pi*fl*(beta_2*config.T_PRI*(0:config.N-1)) - 2*pi*config.fc*tau0 + 2*pi*config.Deltaf*config.T_PRI*(0:(config.L-1))'*(0:config.N-1);
% doppler_mat_2 = exp(1i*omega_dopp);
% inv(doppler_mat_2*doppler_mat_2')
% 
% doppler_mat = doppler_mat_2;
end

