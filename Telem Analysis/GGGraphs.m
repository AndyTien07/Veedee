%% Colored G-g
%volatility=gradient(total_load(Start_Time:End_Time), Start_Time:End_Time);
%volatility=movs    td(total_load(Start_Time:End_Time),50);
%volatility=sqrt(movmean(total_load(Start_Time:End_Time).^2,50));
local_max = movmax(goon,25);
local_min = movmin(goon,25);
volatility= local_max/4 - local_min/4;
lift=PotFlLoad<2 |PotFrLoad<2 |PotRlLoad<2 |PotRrLoad<2;
window=11;

YawR1=sgolayfilt(YawR, 2, window)*d2r;
Yaw_accel=gradient(YawR1, time);
Yaw_moment=Yaw_accel*40;

signals = {'volatility','average_vel','gs','goon','ay','ax','az','PitchPos','RollPos','YawR','Yaw_moment','lift'};

colors = lines(numel(signals));
for k = 1:numel(signals)
    figure('Name', signals{k});
    x = eval(signals{k});

    scatter(latg, ...
    longg, ...
    5, ...
    abs(x),"filled");

    colormap turbo;
    c=colorbar();
    c.Label.String=signals{k};
    xlabel('Lateral acceleration (gs)')
    ylabel('Longitudinal acceleration (gs)')
    fontsize(18,"points")
end


%% G-G+histogram 

fig3 = figure('Color', 'k', 'Name', 'G=Goon');

count = sum(gs(Start_Time:End_Time) < 1);         
total = numel(gs(Start_Time:End_Time));         
percentage = (count / total) * 100 %percentage of things under 1g
sgtitle(['G-G Breakdown: Endurance',round(percentage,2),"% under 1G"] );

subplot(2,2,1);
plot(latg(Start_Time:End_Time),longg(Start_Time:End_Time), '.');

rectangle('Position', [-1, -1, 2, 2], 'Curvature', [1, 1], 'EdgeColor', 'r');

grid on; title('g-g');

subplot(2,2,2);
histogram(latg(Start_Time:End_Time),'Normalization', 'percentage','BinWidth',0.1);
grid on; title('Lateral Gs');
xlabel('gs');
ylabel('Percentage (%)');

subplot(2,2,3);
histogram(longg(Start_Time:End_Time),'Normalization', 'percentage','BinWidth',0.1);
grid on; title('Longitudinal Gs');
xlabel('gs');
ylabel('Percentage (%)');

subplot(2,2,4);
histogram(gs(Start_Time:End_Time),'Normalization', 'percentage','BinWidth',0.1);
grid on; title('Combined Gs');
xlabel('gs');
ylabel('Percentage (%)');

h_zoom = zoom(fig3);
set(h_zoom, 'Motion', 'horizontal', 'Enable', 'on'); 

