%% Load data
DATA_PATH = fullfile("T37 8 14 paccar 5 2 split.csv");
opts = detectImportOptions(DATA_PATH);
data = readtable(DATA_PATH, opts);

%% Speed from wheel speeds
Front_conversion_factor = 1.24482467/(13.25*60);
Rear_conversion_factor  = 1.24482467/(13.25*60);

FR_v = -data.FRrpm * Front_conversion_factor;
FL_v =  data.FLrpm * Front_conversion_factor;
RR_v = -data.RRrpm * Rear_conversion_factor;
RL_v =  data.RLrpm * Rear_conversion_factor;

average_vel = (FR_v + FL_v) / 2;

%% Signals
time   = data.Time;
steer  = data.Steer;      % deg
ay     = data.Ay;
ax = data.Ax;
yawr   = -data.YawR;      % deg/s
L      = 1.535;           % m

%% Convert to radians
steer_rad = deg2rad(steer);
yaw_rad   = deg2rad(yawr);

%% Steady-state mask + low speed cutoff
yaw_accel = gradient(yaw_rad, time);   % rad/s^2
vmin = 2;                              % m/s, tune this
mask = abs(yaw_accel) < 0.2 & abs(average_vel) > vmin;

%% Ackermann-equivalent steer from measured yaw rate
delta_ack = atan2(L .* yaw_rad, average_vel);   % rad

%% Steering error
steer_error = steer_rad - delta_ack;   % >0 usually means understeer

%% Expected yaw rate from steer
r_ack = (average_vel ./ L) .* tan(steer_rad);    % rad/s
yaw_error = abs(yaw_rad) - abs(r_ack);                     % <0 usually means understeer

%% Plots
figure()
scatter(average_vel(mask), yaw_error(mask), 5, steer(mask), 'filled')
xlabel('Velocity')
ylabel('Yaw error (rad/s)')
colorbar

figure()
scatter(steer(mask), yawr(mask), 5, average_vel(mask), 'filled')
xlabel('Steer (deg)')
ylabel('Yaw rate (deg/s)')
colorbar

figure()
hold on
plot(time, yawr)
plot(time, steer)
plot(time, average_vel)
legend('Yaw rate','Steer','Velocity')


steer_rad = deg2rad(steer);
yaw_rad   = deg2rad(yawr);
v         = abs(average_vel);
ay_ms2    = abs(ay);

L = 1.535;
g = 9.81;

% steady-state mask
yaw_accel = gradient(yaw_rad, time);

mask = abs(yaw_accel) < 0.2 & v > 2;

% geometric steer from yaw rate
delta_geo = atan2(L .* yaw_rad, v);   % rad

% understeer gradient, rad/g
K = (steer_rad - delta_geo) ./ (ay_ms2 / g);

% use only steady-state points
K_ss = K(mask);

% average understeer gradient
K_mean = mean(K_ss, 'omitnan');

fprintf('Understeer gradient = %.4f rad/g = %.4f deg/g\n', ...
    K_mean, rad2deg(K_mean));