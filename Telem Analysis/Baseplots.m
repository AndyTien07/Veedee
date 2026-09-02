%% Gs
fig = figure('Color', 'k', 'Name', 'Gs breakdown');
s1 = subplot(4,1,1);
plot(time, longg, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
line([min(time) max(time)], [0 0], 'Color', 'k', 'LineStyle', '--'); % Zero line
ylabel('LongG (Gs)');
grid on; title('Longitudinal Acceleration');

s2 = subplot(4,1,2);
plot(time, latg, 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
line([min(time) max(time)], [0 0], 'Color', 'k', 'LineStyle', '--'); % Zero line
ylabel('LatG (Gs)');
grid on; title('Lateral Acceleration');

s3 = subplot(4,1,3);
plot(time, az/g, 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
line([min(time) max(time)], [0 0], 'Color', 'k', 'LineStyle', '--'); % Zero line
ylabel('Bump G (Gs)');
grid on; title('Vertical Acceleration');

s4 = subplot(4,1,4);
plot(time, gs, 'Color', [0.47 0.67 0.19], 'LineWidth', 1.2);
ylabel('Gs (Total)');
xlabel('Time (seconds)');
grid on; title('Combined G-Sum');

linkaxes([s1, s2, s3, s4], 'x');

h = zoom(fig);
set(h, 'Motion', 'horizontal', 'Enable', 'on'); 

%% LVDT Readout
fig2 = figure('Color', 'k', 'Name', 'Debug');

s1 = subplot(4,1,1);
plot(time, FL_pos*car.front_mr, 'LineWidth', 1.2);
ylabel('pos(m)'); grid on; title('FL');

s2 = subplot(4,1,2);
plot(time, FR_pos*car.front_mr, 'LineWidth', 1.2);
ylabel('pos(m)'); grid on; title('FR');

s3 = subplot(4,1,3);
plot(time, RR_pos*car.rear_mr, 'LineWidth', 1.2);
ylabel('pos(m)'); grid on; title('RR');

s4 = subplot(4,1,4);
plot(time, RL_pos*car.rear_mr, 'LineWidth', 1.2);
ylabel('pos(m)'); grid on; title('RL');
xlabel('Time (s)');

linkaxes([s1, s2, s3, s4], 'x');   
ylim([s1, s2, s3, s4],[-0.015,0.015])
h_zoom = zoom(fig2);
set(h_zoom, 'Motion', 'horizontal', 'Enable', 'on'); 


%% Pitch and roll

fig4 = figure('Color', 'k', 'Name', 'Roll');

s1 = subplot(2,1,1);
%derived from shock pot
plot(time, RollPos, 'LineWidth', 1.2,'DisplayName',"Pot derived");
hold on;
%derived from accel*gradient
plot(time, latg*car.roll_gradient, 'LineWidth', 1.2,'LineStyle', '--','DisplayName',"Accel derived");
ylabel('deg');
lgd1 = legend('TextColor', 'w', 'Color', 'none', 'EdgeColor', 'w');
grid on; title('Roll');

s2 = subplot(2,1,2);
%derived from shock pot
plot(time, PitchPos, 'LineWidth', 1.2,'DisplayName',"Pot derived");
hold on;
%derived from accel*gradient
plot(time, longg*car.pitch_gradient, 'LineWidth', 1.2,'LineStyle', '--','DisplayName',"accel derived");
ylabel('deg');
lgd2 = legend('TextColor', 'w', 'Color', 'none', 'EdgeColor', 'w');
grid on; title('Pitch');
linkaxes([s1, s2], 'x');   
h_zoom = zoom(fig4);
set(h_zoom, 'Motion', 'horizontal', 'Enable', 'on'); 

%% az vs pots
fig9=figure('Color','k','Name','az vs pots')
plot(time, goon,'DisplayName',"Goon");
hold on;
%plot(time,average_vel.*average_vel*2.7,'DisplayName','2.7 CL Downforce');
plot(time,average_vel*100,'DisplayName','mu (scaled)');
plot(time,gs*1000,'DisplayName','Combined gs (scaled)');
plot(time,YawR*10,'DisplayName','YawR (scaled)');
%plot(time,az*car.mass,'DisplayName',"Az Bump force");
legend('TextColor', 'w', 'Color', 'none', 'EdgeColor', 'w');
h_zoom = zoom(fig9);
set(h_zoom, 'Motion', 'horizontal', 'Enable', 'on'); 

%% stacked tire load
fig2 = figure('Color', 'k', 'Name', 'Front Left');

s1 = subplot(4,1,1);
plot(time, PotFlLoad, 'LineWidth', 1.2);
hold on;
plot(time, AccelFlLoad, 'LineWidth', 1.2);
ylabel('FL (N)'); grid on; title('Force');
ylim([0, 1800]);

s2 = subplot(4,1,2);
plot(time, PotFrLoad, 'LineWidth', 1.2);
hold on;
plot(time, AccelFrLoad, 'LineWidth', 1.2);
ylabel('FR (N)'); grid on; title('Force');
ylim([0, 1800]);

s3 = subplot(4,1,3);
plot(time, PotRlLoad, 'LineWidth', 1.2);
hold on;
plot(time, AccelRlLoad, 'LineWidth', 1.2);
ylabel('RL (N)'); grid on; title('Force');
ylim([0, 1800]);

s4 = subplot(4,1,4);
plot(time, PotRrLoad, 'LineWidth', 1.2);
hold on;
plot(time, AccelRrLoad, 'LineWidth', 1.2);
ylabel('RR (N)'); grid on; title('Force');
xlabel('Time (s)');
ylim([0, 1800]);

linkaxes([s1, s2, s3, s4], 'x');   
h_zoom = zoom(fig2);
set(h_zoom, 'Motion', 'horizontal', 'Enable', 'on'); 

%% Mu
fig69 = figure('Color', 'k', 'Name', 'mu');
% subplot(3,1,1);
% plot(time(Start_Time:End_Time),mu(Start_Time:End_Time));

subplot(1,1,1);
histogram(mu,'Normalization', 'percentage', 'BinLimits', [0, 4],'FaceColor','green');
grid on; title('mu');
xlabel('mu');
ylabel('Prercentage');
xlim([0,3]);
ylim([0,10])
% subplot(3,1,3);
% plot(time,total_load);
% xlabel('mu');
% ylabel('Probability Density');
% ylim([1500,4000])

set(h_zoom, 'Motion', 'horizontal', 'Enable', 'on'); 

%% Damper Histograms
% Create histograms for damper data
fig6 = figure('Color', 'k', 'Name', 'Yaw');
ax1=subplot(2,1,1);
histogram(FRVel, 'Normalization', 'percentage');
grid on; title('yaw');
xlabel('Ur Cousin');
ylabel('Percentage');
ax2=subplot(2,1,2);
understeer=abs(ay./average_vel)-abs(YawR*d2r);
understeer(understeer>2)=2;
histogram(understeer(Start_Time:End_Time), 'Normalization', 'percentage','BinWidth',0.025);
grid on; title('Optimal Yaw- Actual Yaw');
xlabel('Radians/s');
ylabel('Percentage');
ylim([0,10])
xlim([-2,2])

