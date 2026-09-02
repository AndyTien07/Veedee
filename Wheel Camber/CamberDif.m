clear; clc; close all;
%Channels to export and rename: 
% aX, aY, aZ, Time, 
% YawR, PitchR, RollR
%Fl,Fr,Rr,Rl
%% 1. Configuration
folder='sad';
DATA_PATH = fullfile(folder, "Pingus.csv");
CAR_PATH= fullfile(folder, 'Car.csv');
opts = detectImportOptions(DATA_PATH);
data = readtable(DATA_PATH, opts);
opts = detectImportOptions(CAR_PATH);
car=readtable(CAR_PATH,opts);
CUTOFF_FREQ = 20; %for lowpass 


%% 2. Data Import
time = data.("Time");
dt=0.01
[b, a] = butter(2, CUTOFF_FREQ / ((1/(time(2)-time(1)))/ 2), 'low');

d2r=1;

YawR=filtfilt(b, a,data.("YawR"))*(1/d2r);
PitchR=filtfilt(b, a,data.("PitchR"))*(1/d2r);
RollR=filtfilt(b, a, data.("RollR"))*(1/d2r);

ax = filtfilt(b, a, data.("Ax"));
%ax=data.("Ax");
ay =filtfilt(b, a, data.("Ay"));
%ay=data.("Ay");
az =filtfilt(b, a, data.("Az"));
%az=data.("Az");

IMU_Pitch=filtfilt(b,a,data.PitchPos);
IMU_Roll=filtfilt(b,a,-data.RollPos);


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
average_vel=FR_v+FL_v
average_vel=average_vel/2;

%% Pitch Roll and Yaw
%DO NOT TOUCH THIS IT IS TUNED FOR BILATERAL FILTER
alpha = 0.05; %DO NOT TOUCH THIS IT IS TUNED
% DO NOT TOUCH THIS IT IS TUNED 0.95 good

YawPos = zeros(length(time), 1);
PitchPos = zeros(length(time), 1);
RollPos = zeros(length(time), 1);

for i = 2:length(time)
    % Roll Reference (Average of front and rear roll angles)
    % Changed atand to atan
    RollLvdtF = atand(((FR_pos(i)*car.front_mr - FL_pos(i))*car.front_mr) / car.front_track);
    roll_lvdt_r = atand((RR_pos(i)*car.rear_mr - RL_pos(i)*car.rear_mr) / car.rear_track);
    roll_ref = (roll_lvdt_r);
    
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

% Identify anchor indices (near-zero lateral acceleration)
anchors = find(abs(ay) < 0.5);

% Ensure row/column consistency (force IMU_Roll to column vector)
IMU_Roll_raw = IMU_Roll(:);
IMU_Roll_corr = IMU_Roll_raw;

% Detrend segments between anchors using original uncorrected references
for i = 2:length(anchors)
    idx = anchors(i-1) : (anchors(i) - 1);
    
    if length(idx) > 1
        % Generate trendline matching vector shape
        trend = linspace(IMU_Roll_raw(anchors(i-1)), IMU_Roll_raw(anchors(i)), length(idx))';
        IMU_Roll_corr(idx) = IMU_Roll_raw(idx) - trend;
    else
        IMU_Roll_corr(idx) = 0;
    end
end

% Assign back to variable
IMU_Roll = IMU_Roll_corr;

% Linear fits against lateral acceleration
susfit  = polyfit(ay, RollPos, 1);
susgrad = susfit(1) * 9.81;

allfit  = polyfit(ay, IMU_Roll, 1);
allgrad = allfit(1) * 9.81;

% Plotting results
x = linspace(min(ay), max(ay), 100);
figure;
hold on;
plot(ay, RollPos, '.', 'DisplayName', 'Suspension Roll');
plot(x, polyval(susfit, x), 'r-', 'LineWidth', 1.5, 'DisplayName', sprintf('Suspension Fit (Grad: %.2f)', susgrad));

plot(ay, IMU_Roll, '.', 'DisplayName', 'Corrected IMU Roll');
plot(x, polyval(allfit, x), 'b-', 'LineWidth', 1.5, 'DisplayName', sprintf('Total Fit (Grad: %.2f)', allgrad));

grid on;
xlabel('Lateral Acceleration a_y [m/s^2]');
ylabel('Roll Angle [deg]');
legend('Location', 'best');
figure()
hold on;
plot(time,IMU_Roll)
plot(time,RollPos)

figure()
hold on;
plot(ay,IMU_Roll,'.','DisplayName','IMU Derived Roll','Color','white')
plot(x,polyval(allfit,x),'w-','DisplayName',sprintf('IMU Roll Gradient %g',allgrad))
plot(ay,RollPos,'.','DisplayName','Shockpot Derived Roll','Color','cyan')
plot(x,polyval(susfit,x),'c-','DisplayName',sprintf('Shockpot Roll Gradient %g',susgrad))
legend('show','Location', 'southeast')
fontsize(18,'points')