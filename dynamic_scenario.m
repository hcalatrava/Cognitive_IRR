% SCENARIO WITH TRAJECTORIES (September 2023) - Assumed known velocity
%% Configuration
close all,
clear all,
clc,
config = load_config(5); % Load configuration

N_uav = 6;
N_t = 20; % total number of seconds in the simulation
% Position at which spoofer joins the system
config.scenario.target_dist_north = config.target_dist_north; % in meters
config.scenario.target_dist_east = config.target_dist_east; 
config.X_H1 = diag(sqrt(config.varX/2)*(randn(config.L,1) + 1i*randn(config.L,1))); % LxL complex diagonal matrix repressenting the scattering coefficients (target)


% Different UAVs velocities (we assume first one is the correct trajectory)
config.scenario.uav_velocity(1,:) = config.target_velocity'; % target velocity in m/s
config.scenario.uav_velocity(2,:) = [-50, 4]'; % target velocity in m/s
config.scenario.uav_velocity(3,:) = [1, 4]'; % target velocity in m/s
config.scenario.uav_velocity(4,:) = [7.07, 1]'; % target velocity in m/s
config.scenario.uav_velocity(5,:) = [-2, 4]'; % target velocity in m/s
config.scenario.uav_velocity(6,:) = [7.07, -3]'; % target velocity in m/s
config.scenario.colors = ["g", "r", "b", "m", "y", "c"]

position = zeros(N_uav,N_t,2); % 2D north-east coords

for uav = 1:N_uav
    position(uav,:,1) = config.scenario.target_dist_north + (0:N_t-1)*config.scenario.uav_velocity(uav,1); % position w.r.t. the target in the north direction
    position(uav,:,2) = config.scenario.target_dist_east + (0:N_t-1)*config.scenario.uav_velocity(uav,2); % position w.r.t the target in the east direction
end

%% Draw scenario while detecting
f = figure()
f.Position = [100 100 700 600];
min_x = min(position(:,:,1),[],'all');
min_y = min(position(:,:,2),[],'all');
xlim([min(min_x,0), max(position(:,:,1),[],'all')])
ylim([min(min_y,0), max(position(:,:,2),[],'all')])
hold on,
config.target_velocity = config.target_velocity'; % constant velocity model
for t = 1:N_t

    % For each time instant, we assume the velocity and DOA of the target
    % to be known, so we calculate the detections with the velocity and LOS
    % vector given by the true position and velocity of the target (i.e.,
    % UAV 1)
    config.target_dist_north = position(1,t,1);
    config.target_dist_east = position(1,t,2);
    LOS_vector = get_LOS_vector(config);

    for uav = 1:N_uav
        plot(position(uav,1:t,1), position(uav,1:t,2), 'marker', '*', 'color', config.scenario.colors(uav), 'markersize', 20, 'linestyle', 'none')
        hold on,

        detected = check_detection(config, position(1,t,:), position(uav,t,:), config.scenario.uav_velocity(uav,:), LOS_vector);
        if detected
            plot(position(uav,t,1), position(uav,t,2), 'marker', 'o','color', config.scenario.colors(uav), 'markersize', 20, 'linestyle', 'none', 'HandleVisibility','off')
            hold on,
        end
    end
    pause(0.1)
end
plot(0, 0, 'color', 'black', 'marker', '.', 'markersize', 30)
text(0.5, 1, 'RADAR', 'fontsize', 18)
hold on,
plot(0, 0, 'color', 'black', 'marker', 'o', 'markersize', 30)
text(0.5, 1, 'RADAR', 'fontsize', 18)
plot(position(uav,1,1), position(uav,1,2), 'color', 'black', 'marker', '*', 'markersize', 20)
plot(position(uav,1,1), position(uav,1,2), 'color', 'black', 'marker', 'o', 'markersize', 20)
text(position(uav,1,1)+2, position(uav,1,2)-2, 'Spoofing Starts', 'fontsize', 18)
grid('on')
xlabel('Dist. North direction [m]','fontsize', 17)
ylabel('Dist. East direction [m]', 'fontsize', 17)
title('2D Target under Spoofer Scenario', 'fontsize', 20)

string_legend = ["Target"];
for uav = 1:(N_uav-1)
    string_legend = [string_legend, ['Spoofer ', num2str(uav)]];
end     
legend(string_legend, 'fontsize', 15, 'location', 'bestoutside')


