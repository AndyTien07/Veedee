%% Load data
DATA_PATH = fullfile("Tag 9 v2.csv");
opts = detectImportOptions(DATA_PATH);
data = readtable(DATA_PATH, opts);

%% Speed from wheel speeds
Front_conversion_factor = 1.24482467/(13.25*60);
Rear_conversion_factor  = 1.24482467/(13.25*60);

FR_v = -data.FRrpm * Front_conversion_factor;
FL_v =  data.FLrpm * Front_conversion_factor;
RR_v = -data.RRrpm * Rear_conversion_factor;
RL_v =  data.RLrpm * Rear_conversion_factor;

average_vel = (FR_v + FL_v + RR_v + RL_v) / 4;

%% Signals
time   = data.Time;
ay     = data.Ay;
