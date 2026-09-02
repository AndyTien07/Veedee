%% Visualize steady-state sweep results
% Run steady_state_sweep.m first so ayMap, yawMomentMap, convergedMap,
% Velocities, SteerAnglesDeg, BodySlipAnglesDeg, and velocityPeaks exist.

g = 9.80665;

requiredVariables = {'ayMap','yawMomentMap','convergedMap', ...
    'Velocities','SteerAnglesDeg','BodySlipAnglesDeg'};
for k = 1:numel(requiredVariables)
    assert(exist(requiredVariables{k},'var') == 1, ...
        'Run steady_state_sweep.m first. Missing variable: %s', ...
        requiredVariables{k});
end

[BetaGrid,SteerGrid] = meshgrid(BodySlipAnglesDeg,SteerAnglesDeg);


%% Combined yaw-moment diagram for every speed
allAyG = ayMap/g;
allYawMoment = yawMomentMap;
allValid = convergedMap & isfinite(allAyG) & isfinite(allYawMoment);

% Expand the velocity vector to the same size as the result maps.
velocityMap = repmat(reshape(Velocities,[],1,1), ...
    1,numel(SteerAnglesDeg),numel(BodySlipAnglesDeg));

figure('Name','Yaw Moment Diagram - All Speeds','Color','k');
scatter(allAyG(allValid),allYawMoment(allValid),24, ...
    velocityMap(allValid),'filled', ...
    'MarkerFaceAlpha',0.65,'MarkerEdgeAlpha',0.25);
hold on;
yline(0,'k-','LineWidth',1.8);
xline(0,'k:','LineWidth',1.0);
xlabel('Lateral acceleration [g]');
ylabel('Net yaw moment [N m]');
title('Yaw-Moment Diagram for All Speeds');
cb = colorbar;
cb.Label.String = 'Vehicle speed [m/s]';
colormap(turbo);
clim([min(Velocities),max(Velocities)]);
grid on;
box on;
fontsize(18,"points")

%% Maximum converged ay and maximum zero-moment ay versus speed
maxConvergedAyG = NaN(size(Velocities));
maxAyGMz=NaN(size(Velocities));
maxSteadyAyG = NaN(size(Velocities));
steadyMomentTolerance = 5; % [N m]

for velocityIndex = 1:numel(Velocities)
    localAyG = squeeze(ayMap(velocityIndex,:,:))/g;
    localMoment = squeeze(yawMomentMap(velocityIndex,:,:));
    localValid = squeeze(convergedMap(velocityIndex,:,:));

    allPoints = localValid & isfinite(localAyG) & isfinite(localMoment);
    if any(allPoints(:))
        maxConvergedAyG(velocityIndex) = max(localAyG(allPoints));
         maxAyGMz(velocityIndex) = localMoment(localAyG(allPoints)==maxConvergedAyG(velocityIndex));
    end

    steadyPoints = allPoints & abs(localMoment) <= steadyMomentTolerance;
    if any(steadyPoints(:))
        maxSteadyAyG(velocityIndex) = max(localAyG(steadyPoints));
    end
end

figure('Name','Lateral Acceleration versus Speed','Color','k');

ax(1)=subplot(2,1,1);
plot(Velocities,maxConvergedAyG,'o-','LineWidth',1.8, ...
    'DisplayName','Maximum converged');
hold on;
plot(Velocities,maxSteadyAyG,'s-','LineWidth',1.8, ...
    'DisplayName','Maximum Steady State');
xlabel('Velocity [m/s]');
ylabel('Lateral acceleration [g]');
ylim([1.5,3])
title('Maximum Lateral Acceleration versus Speed');
legend('Location','best');
ax(2)=subplot(2,1,2);
plot(Velocities,maxAyGMz);
title('Residual Yaw Moment at Peak AY vs speed');

% Finalize the plot with grid and labels
grid on;
xlabel('Velocity [m/s]');
ylabel('Residual Yaw Moment [N m]');
grid on;
linkaxes(ax,'x');
zoom xon;
fontsize(18,'points')

