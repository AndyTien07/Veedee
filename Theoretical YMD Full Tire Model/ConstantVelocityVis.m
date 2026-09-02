%% Visualize steady-state sweep results
% Run steady_state_sweep.m first so the result maps and sweep vectors exist.

g = 9.80665;

requiredVariables = {'ayMap','radiusMap','yawMomentMap','convergedMap', ...
    'Velocities','SteerAnglesDeg','BodySlipAnglesDeg'};
for k = 1:numel(requiredVariables)
    assert(exist(requiredVariables{k},'var') == 1, ...
        'Run steady_state_sweep.m first. Missing variable: %s', ...
        requiredVariables{k});
end

[BetaGrid,SteerGrid] = meshgrid(BodySlipAnglesDeg,SteerAnglesDeg);

%% Choose a velocity to inspect
plotVelocity =10; % [m/s]
[~,velocityIndex] = min(abs(Velocities-plotVelocity));
actualVelocity = Velocities(velocityIndex);

ayG = squeeze(ayMap(velocityIndex,:,:))/g;
yawMoment = squeeze(yawMomentMap(velocityIndex,:,:));
radius = abs(squeeze(radiusMap(velocityIndex,:,:)));
valid = squeeze(convergedMap(velocityIndex,:,:));

ayG(~valid) = NaN;
yawMoment(~valid) = NaN;
radius(~valid | ~isfinite(radius) | radius <= 0) = NaN;

hasZeroMomentContour = any(isfinite(yawMoment(:))) && ...
    min(yawMoment(:),[],'omitnan') <= 0 && ...
    max(yawMoment(:),[],'omitnan') >= 0;

%% Steering/body-slip maps
maximumRadiusColor = 11; % [m], clips color only
minimumRadiusColor=9; % [m], clips color only
radiusForPlot = min(radius,maximumRadiusColor);

figure('Name','Steady-State Steering and Body-Slip Maps','Color','k');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile;
contourf(BetaGrid,SteerGrid,ayG,20,'LineColor','none');
colorbar;
hold on;
if hasZeroMomentContour
    contour(BetaGrid,SteerGrid,yawMoment,[0 0], ...
        'k','LineWidth',2.5,'DisplayName','M_z = 0');
end
xlabel('Body slip \beta [deg]','FontSize',18);
ylabel('Handwheel angle [deg]','FontSize',18);
title(sprintf('Lateral Acceleration at %.1f m/s [g]',actualVelocity));
grid on;

nexttile;
contourf(BetaGrid,SteerGrid,yawMoment,20,'LineColor','none');
colorbar;
hold on;
if hasZeroMomentContour
    contour(BetaGrid,SteerGrid,yawMoment,[0 0],'k','LineWidth',2.5);
end
xlabel('Body slip \beta [deg]');
ylabel('Handwheel angle [deg]');
title(sprintf('Net Yaw Moment at %.1f m/s [N m]',actualVelocity));
colormap(turbo);
clim([-1000 1000]);
grid on;

nexttile;
contourf(BetaGrid,SteerGrid,radiusForPlot,20,'LineColor','none');
cb = colorbar;
cb.Label.String = sprintf('Turn radius [m], clipped at %.0f m', ...
    maximumRadiusColor);
hold on;
if hasZeroMomentContour
    contour(BetaGrid,SteerGrid,yawMoment,[0 0],'k','LineWidth',2.5);
end
xlabel('Body slip \beta [deg]');
ylabel('Handwheel angle [deg]');
title(sprintf('Turn Radius at %.1f m/s',actualVelocity));
colormap(turbo);
clim([minimumRadiusColor maximumRadiusColor]);
grid on;
fontsize(18,"points")