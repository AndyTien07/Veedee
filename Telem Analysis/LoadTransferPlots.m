figure();
hold on
scatter(RollR2,abs((AccelFrLoad-AccelFlLoad)-(PotFrLoad-PotFlLoad)),5,abs(ay));
ylabel('Theoretical LT - Pot LT (N)');
xlabel('Roll Rate (deg/sec)');
colormap turbo;
c=colorbar();
c.Label.String='Lateral Acceleration (gs)';

% ==========================================
% ADDED: AY TREND ISOLINES
% ==========================================
X_data = RollR2;
Y_data = abs((AccelFrLoad-AccelFlLoad)-(PotFrLoad-PotFlLoad));
LATg_data = abs(latg);

% Define the specific ay steps you want to track (matching the turbo colorbar)
ay_levels = 1:0.2:2.5; 
tolerance = .05; % Window to capture enough points for each slice (+/- m/s^2)

for i = 1:length(ay_levels)
    target_ay = ay_levels(i);
    
    % Mask for the current ay band
    idx = abs(LATg_data - target_ay) < tolerance;
    
    if sum(idx) > 20 % Ensure there's enough data to build a trend
        x_slice = X_data(idx);
        y_slice = Y_data(idx);
        
        % Fit a 2nd-order\ polynomial to capture the "V" or "U" curve of the slice
        p = polyfit(x_slice, y_slice, 2);
        
        % Generate a smooth x-axis spanning the width of this specific data slice
        x_trend = linspace(min(x_slice), max(x_slice), 100);
        y_trend = polyval(p, x_trend);
        
        % Plot a thick white backing line so the trendline pops over the dense dots
        plot(x_trend, y_trend, 'w-', 'LineWidth', 3.5);
        
        % Plot the thin black dashed trendline
        plot(x_trend, y_trend, 'k--', 'LineWidth', 1.5);
        
        % Label the trendline at its furthest right point
        text(x_trend(end), y_trend(end), [num2str(target_ay) ' gs'], ...
            'Color', 'white', 'FontSize', 9, 'FontWeight', 'bold', ...
            'BackgroundColor', [0 0 0 0.5], 'Margin', 1);
    end
end
% ==========================================

hold off

wtf= abs(az/g)<0.1 & average_vel<20 & average_vel>5;
figure()
scatter(AccelTotalLoad(wtf)-PotTotalLoad(wtf),TotalLT(wtf),5,RollPos(wtf))
figure()
scatter(PotTotalLoad(wtf),TotalLT(wtf),5,max_load(wtf))
ylabel('Total Lateral Load Transfer (N)');
xlabel('Measured Load from Pots (N)');
c=colorbar();
c.Label.String='Max Load on a Tire (N)';

figure()
scatter(PotTotalLoad(wtf),LongLT(wtf),5,az(wtf)/g)