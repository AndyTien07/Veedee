clear; clc; close all;
%Channels to export and rename: 
% aX, aY, aZ, Time, 
% YawR, PitchR, RollR
%Fl,Fr,Rr,Rl
%% 1. Configuration
folder='Paccar 6-6';
DATA_PATH = fullfile(folder, "Tag 38.csv");
CAR_PATH= fullfile(folder, 'Car.csv');
opts = detectImportOptions(DATA_PATH);
data = readtable(DATA_PATH, opts);
opts = detectImportOptions(CAR_PATH);
car=readtable(CAR_PATH,opts);
CUTOFF_FREQ = 20; %for lowpass 

Start_Time=4*100;
End_Time=10*100;

%% 2. Data Import
time = data.("Time");
dt=0.01
[b, a] = butter(2, CUTOFF_FREQ / ((1/(time(2)-time(1)))/ 2), 'low');

d2r=1; %change to actual value if log is in degrees 

YawR=filtfilt(b, a,data.("YawR"))*(1/d2r);
PitchR=filtfilt(b, a,data.("PitchR"))*(1/d2r);
RollR=filtfilt(b, a, data.("RollR"))*(1/d2r);

ax = filtfilt(b, a, data.("Ax"));
%ax=data.("Ax");
ay =filtfilt(b, a, data.("Ay"));
%ay=data.("Ay");
az =filtfilt(b, a, data.("Az"));
%az=data.("Az");

window=11;
zero_start=8*100;
zero_end=12*100;
%FL_pos = filtfilt(b, a, (data.FL * car.lvdt_conversion)); 
FL_pos = data.FL * car.lvdt_conversion; 
FL_pos = FL_pos - mean(FL_pos(zero_start:zero_end));

%FR_pos = filtfilt(b, a, data.FR * car.lvdt_conversion); 
FR_pos = data.FR * car.lvdt_conversion; 
FR_pos = FR_pos - mean(FR_pos(zero_start:zero_end));

%RL_pos = filtfilt(b, a, (data.RL *car.lvdt_conversion)); 
RL_pos = data.RL * car.lvdt_conversion; 
RL_pos = RL_pos - mean(RL_pos(zero_start:zero_end));


%RR_pos = filtfilt(b, a, (data.RR * car.lvdt_conversion)); 
RR_pos = data.RR * car.lvdt_conversion; 
RR_pos = RR_pos - mean(RR_pos(zero_start:zero_end));

%% car/wheel speed
Front_conversion_factor=1.24482467/(13.25*60);
Rear_conversion_factor=1.24482467/(13.88*60);

FR_v=-data.FRrpm*Front_conversion_factor;
FL_v=data.FLrpm*Front_conversion_factor;
RR_v=-data.RRrpm*Rear_conversion_factor;
RL_v=data.RLrpm*Rear_conversion_factor;
average_vel=FR_v+FL_v+RR_v+RL_v;
average_vel=average_vel/4;

%% Pitch Roll and Yaw
%DO NOT TOUCH THIS IT IS TUNED FOR BILATERAL FILTER
alpha = 0.05; %DO NOT TOUCH THIS IT IS TUNED
% DO NOT TOUCH THIS IT IS TUNED 0.05 good

YawPos = zeros(length(time), 1);
PitchPos = zeros(length(time), 1);
RollPos = zeros(length(time), 1);

for i = 2:length(time)
    % Roll Reference (Average of front and rear roll angles)
    % Changed atand to atan
    RollLvdtF = atand(((FR_pos(i) - FL_pos(i))*car.front_mr) / car.front_track);
    rollLvdtR = atand((RR_pos(i) - RL_pos(i)) / car.rear_track);
    roll_ref = (RollLvdtF+rollLvdtR)/2;
    
    % Pitch Reference
    avg_f = (FR_pos(i) + FL_pos(i))*car.front_mr / 2;
    avg_r = (RL_pos(i) + RR_pos(i)) / 2;
    % Changed atand to atan
    pitch_ref = atand((avg_r-avg_f) / car.wheelbase);
    
    %% Apply Complementary Filter
    % Roll Fusion (All terms now in radians or radians/sec * sec)
    RollPos(i) = alpha * (RollPos(i-1) + RollR(i) * dt) + (1 - alpha) * roll_ref;
    
    % Pitch Fusion
    PitchPos(i) = alpha * (PitchPos(i-1) + PitchR(i) * dt) + (1 - alpha) * pitch_ref;
    
    % Yaw Integration (Pure Dead Reckoning)
    YawPos(i) = YawPos(i-1) + YawR(i) * dt;
end

PitchR2 = [0; diff(PitchPos) ./ diff(time)]; % Recommended: divide by dt if time step is explicitly needed
RollR2  = [0; diff(RollPos) ./ diff(time)]; 

%% Accelerometer rectification using Pitch Roll/Yaw
ax = ax - car.imuz_offset * PitchR*d2r + car.imuy_offset * YawR*d2r;
ay = ay - car.imux_offset * YawR*d2r   + car.imuz_offset * RollR*d2r;
az = az - car.imuy_offset * RollR*d2r + car.imux_offset * PitchR*d2r;

g=9.81;
gs=sqrt(ax.^2+ay.^2)/g;
latg=ay/g;
longg=ax/g;

%% Calculating Tire Loads

%Pull from Optimum G LT calculator
C_front = 5.2; 
C_rear  = 4.10;
%Aero
cl=2.7;
area=1.18;
rho=1.204;

FLVel = gradient(FL_pos, time);
DampFL = C_front .* FLVel;
FRVel = gradient(FR_pos, time);
DampFR = C_front .* FRVel;
RLVel = gradient(RL_pos, time);
DampRL = C_rear .* RLVel;
RRVel = gradient(RR_pos, time);
DampRR = C_rear .* RRVel;

front_mass=car.FL_corner_weight+car.FR_corner_weight;
%front_mass=front_mass/g;
rear_mass=car.RL_corner_weight+car.RR_corner_weight;
%rear_mass=rear_mass/g;
front_sprung_mass=car.sprung*car.weight_dist_f;
rear_sprung_mass=car.sprung*(1-car.weight_dist_f);

long_transfer = (car.mass * ax * car.cg_height) / car.wheelbase;
lat_transfer_f = (ay * front_mass * car.cg_height) / car.front_track;
lat_transfer_r = (ay * rear_mass * car.cg_height) / car.rear_track;
tire_bump = az*car.mass;

ElasticWTf=car.sprung*0.52*ay*0.05/car.front_track;
ElasticWTr=car.sprung*0.48*ay*0.07/car.rear_track;

%% Calculate the total acceleration load
% Pre-allocate
num_pts = length(time);

PotFlLoad = zeros(num_pts, 1);
PotFrLoad = zeros(num_pts, 1);
PotRlLoad = zeros(num_pts, 1);
PotRrLoad = zeros(num_pts, 1);
PotTotalLoad= zeros(num_pts, 1);

AccelFlLoad = zeros(num_pts, 1);
AccelFrLoad = zeros(num_pts, 1);
AccelRlLoad = zeros(num_pts, 1);
AccelRrLoad = zeros(num_pts, 1);
AccelTotalLoad= zeros(num_pts, 1);

for i = 2:num_pts
    PotFlLoad(i) = (FL_pos(i) * car.front_rate * car.front_mr +DampFL(i)) + (car.FL_corner_weight* g);;
    PotFrLoad(i) = (FR_pos(i) * car.front_rate * car.front_mr +DampFR(i)) + (car.FR_corner_weight* g);
    PotRlLoad(i) = (RL_pos(i) * car.rear_rate * car.rear_mr +DampRL(i)) + (car.RL_corner_weight* g );
    PotRrLoad(i) = (RR_pos(i) * car.rear_rate * car.rear_mr +DampRR(i)) + (car.RR_corner_weight* g );
    PotTotalLoad(i) = max(PotFlLoad(i),0)+max(PotFrLoad(i),0)+max(PotRlLoad(i),0)+max(PotRrLoad(i),0);
    
    AccelFrLoad(i)=(front_mass * g / 2) + tire_bump(i)/4 - (long_transfer(i) / 2)+ lat_transfer_f(i);
    AccelFlLoad(i)=(front_mass * g / 2) + tire_bump(i)/4 - (long_transfer(i) / 2)- lat_transfer_f(i);
    AccelRrLoad(i)=(rear_mass * g / 2) + tire_bump(i)/4 + (long_transfer(i) / 2) + lat_transfer_r(i);
    AccelRlLoad(i)=(rear_mass * g / 2) + tire_bump(i)/4 + (long_transfer(i) / 2) - lat_transfer_r(i);
    AccelTotalLoad(i) = AccelFlLoad(i) + AccelFrLoad(i) + AccelRlLoad(i) + AccelRrLoad(i)+average_vel(i)*average_vel(i)*cl*area*rho;

    
end
FrontLT=PotFrLoad-PotFlLoad;
TotalLT=(PotFrLoad+PotRrLoad)-(PotRlLoad+PotFlLoad);
LongLT=(PotFrLoad+PotFlLoad)-(PotRlLoad+PotRrLoad);
max_load=max([PotRrLoad, PotRlLoad, PotFrLoad, PotFlLoad],[],2);
max_vel=max([abs(RLVel), abs(RRVel), abs(FRVel), abs(FLVel)],[],2);
mu = (gs* g * car.mass) ./AccelTotalLoad;

%% 4. Visualization
% call any visualization function you want to use
loadSensitivityTest
