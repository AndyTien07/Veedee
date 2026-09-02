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
steadyMomentTolerance = 1; % [N m]

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

%% Required steer angle to reach steady state at each g increment
% Strategy: for each velocity and each steer angle, scan across body-slip
% angle and find where the net yaw moment crosses zero (the "trim" point,
% i.e. steady-state cornering at that steer angle). Interpolate to get
% the exact trim ay and store it. This produces a trim curve of ay vs.
% steer angle for each speed, which is then inverted (interpolated) to
% find the steer angle required to trim at fixed g increments.

gIncrement = 0.01;                 % [g] resolution for output table
nVelocity = numel(Velocities);
nSteer = numel(SteerAnglesDeg);

% Trim ay [g] achieved at each (velocity, steer angle) - first zero
% crossing of yaw moment scanning from beta = 0 upward.
trimAyG = NaN(nVelocity,nSteer);

for velocityIndex = 1:nVelocity
    for steerIndex = 1:nSteer
        localMoment = squeeze(yawMomentMap(velocityIndex,steerIndex,:));
        localAyG = squeeze(ayMap(velocityIndex,steerIndex,:))/g;
        localValid = squeeze(convergedMap(velocityIndex,steerIndex,:));

        valid = localValid & isfinite(localMoment) & isfinite(localAyG);
        validIdx = find(valid);
        if numel(validIdx) < 2
            continue
        end

        momentSeq = localMoment(validIdx);
        aySeq = localAyG(validIdx);

        % Look for the first sign change of yaw moment along increasing
        % body-slip angle (already ascending order in BodySlipAnglesDeg).
        signChangeIdx = find(diff(sign(momentSeq)) ~= 0,1,'first');
        if isempty(signChangeIdx)
            continue
        end

        m1 = momentSeq(signChangeIdx);
        m2 = momentSeq(signChangeIdx+1);
        a1 = aySeq(signChangeIdx);
        a2 = aySeq(signChangeIdx+1);

        if m2 == m1
            continue
        end
        interpFraction = -m1/(m2-m1);
        trimAyG(velocityIndex,steerIndex) = a1 + interpFraction*(a2-a1);
    end
end

% Determine common ay-increment axis spanning all trimmed data.
maxTrimAy = max(trimAyG(:),[],'omitnan');
if isempty(maxTrimAy) || ~isfinite(maxTrimAy)
    warning('No steady-state trim points found; skipping steer requirement table.');
else
    ayLevels = 0:gIncrement:floor(maxTrimAy/gIncrement)*gIncrement;
    nAyLevels = numel(ayLevels);

    % steerRequiredDeg(v,i) = handwheel steer angle needed to trim at
    % ayLevels(i) for Velocities(v). NaN where no trim exists at that
    % speed/ay combination.
    steerRequiredDeg = NaN(nVelocity,nAyLevels);

    for velocityIndex = 1:nVelocity
        curveAy = trimAyG(velocityIndex,:);
        curveSteer = SteerAnglesDeg;

        validCurve = isfinite(curveAy);
        if nnz(validCurve) < 2
            continue
        end

        % interp1 requires monotonic sample points - sort by ay.
        [sortedAy,sortOrder] = sort(curveAy(validCurve));
        sortedSteer = curveSteer(validCurve);
        sortedSteer = sortedSteer(sortOrder);

        % Remove duplicate ay values (keep first) to satisfy interp1.
        [sortedAy,uniqueIdx] = unique(sortedAy,'stable');
        sortedSteer = sortedSteer(uniqueIdx);

        if numel(sortedAy) < 2
            continue
        end

        inRange = ayLevels >= min(sortedAy) & ayLevels <= max(sortedAy);
        steerRequiredDeg(velocityIndex,inRange) = interp1( ...
            sortedAy,sortedSteer,ayLevels(inRange),'linear');
    end

    %% Plot: required steer angle vs lateral g, one line per speed
    figure('Name','Required Steer Angle for Steady State','Color','k');
    hold on;
    colors = turbo(nVelocity);
    for velocityIndex = 1:nVelocity
        plot(ayLevels,steerRequiredDeg(velocityIndex,:),'-o', ...
            'LineWidth',1.4,'MarkerSize',3, ...
            'Color',colors(velocityIndex,:), ...
            'DisplayName',sprintf('%.0f m/s',Velocities(velocityIndex)));
    end
    xlabel('Lateral acceleration [g]');
    ylabel('Required handwheel steer angle [deg]');
    title('Steer Angle Required for Steady-State Cornering');
    cb = colorbar;
    colormap(turbo);
    clim([min(Velocities),max(Velocities)]);
    cb.Label.String = 'Vehicle speed [m/s]';
    grid on;
    box on;
    fontsize(16,'points')

end
