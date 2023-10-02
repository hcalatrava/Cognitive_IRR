function LOS_vector = get_LOS_vector(config)
%GET_LOS_VECTOR returns the 2D LOS vector

u1 = [config.target_dist_east; config.target_dist_north];
u1 = u1./norm(u1);
LOS_vector = u1;

% wgs84 = wgs84Ellipsoid('kilometer');
% 
% % Random ECEF position from https://www.mathworks.com/help/map/ref/ecef2ned.html
% radar_x = 1345.660;
% radar_y = -4350.891;
% radar_z = 4452.314;
% radar_pos_ECEF = [radar_x, radar_y, radar_z];
% 
% % Specify the origin of the local NED system with the geodetic coordinates lat0, lon0, and h0
% lat0 = 44.532;
% lon0 = -72.782;
% h0 = 1.699;
% 
% % Compute radar NED coordinates 
% [radar_xNorth, radar_yEast, radar_zDown] = ecef2ned(radar_x,radar_y,radar_z,lat0,lon0,h0,wgs84)
% 
% % Compute target NED coordinates considering the provided distance values
% % (given in NED coordinates)
% target_xNorth = radar_xNorth+config.target_dist_north/1e3;
% target_yEast = radar_yEast+config.target_dist_east/1e3;
% target_zDown = radar_zDown;
% 
% % Get target ECEF coordinates
% [target_x, target_y, target_z] = ned2ecef(target_xNorth, target_yEast, target_zDown, lat0,lon0,h0,wgs84);
% 
% % Only in 2D
% % target_pos_ECEF = [target_x, target_y, target_z];
% target_pos_ECEF = [target_x, target_y];
% radar_pos_ECEF = radar_pos_ECEF(1:2);
% 
% % Compute LOS vector
% LOS_vector = (radar_pos_ECEF - target_pos_ECEF)/norm((radar_pos_ECEF - target_pos_ECEF));
% LOS_vector = LOS_vector(1:2);