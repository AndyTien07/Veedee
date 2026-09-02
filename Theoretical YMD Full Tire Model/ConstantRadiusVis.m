%% Constant-radius test visualization (SAE J266-style)
% Run steady_state_sweep.m first so ayMap, radiusMap, yawMomentMap,
% convergedMap, Velocities, SteerAnglesDeg, and BodySlipAnglesDeg exist.
%
% Strategy: for each velocity and each steer angle, scan across body-slip
% angle to find the steady-state trim point (net yaw moment = 0). At that
% trim point, record the trim lateral acceleration [g] and trim turn
% radius [m]. This gives, per velocity, a trim curve relating steer angle
% -> (ay, radius). For each target radius requested, interpolate along
% that trim curve to find the steer angle and lateral g required to hold
% that radius at each speed -- the standard constant-radius test output.

g = 9.80665;

requiredVariables = {'ayMap','radiusMap','yawMomentMap','convergedMap', ...
    'Velocities','SteerAnglesDeg','BodySlipAnglesDeg'};
for k = 1:numel(requiredVariables)
    assert(exist(requiredVariables{k},'var') == 1, ...
        'Run steady_state_sweep.m first. Missing variable: %s', ...
        requiredVariables{k});
end

targetRadii = [8.25,10,11,12,16,20];   % [m] radii to evaluate
nVelocity = numel(Velocities);
nSteer = numel(SteerAnglesDeg);
nRadii = numel(targetRadii);

%% Step 1: build trim curves (steer -> ay, radius) at each velocity
trimAyG = NaN(nVelocity,nSteer);
trimRadius = NaN(nVelocity,nSteer);

for velocityIndex = 1:nVelocity
    for steerIndex = 1:nSteer
        localMoment = squeeze(yawMomentMap(velocityIndex,steerIndex,:));
        localAyG = squeeze(ayMap(velocityIndex,steerIndex,:))/g;
        localRadius = squeeze(radiusMap(velocityIndex,steerIndex,:));
        localValid = squeeze(convergedMap(velocityIndex,steerIndex,:));

        valid = localValid & isfinite(localMoment) & isfinite(localAyG) ...
            & isfinite(localRadius);
        validIdx = find(valid);
        if numel(validIdx) < 2
            continue
        end

        momentSeq = localMoment(validIdx);
        aySeq = localAyG(validIdx);
        radiusSeq = localRadius(validIdx);

        % First sign change of yaw moment along increasing body-slip.
        signChangeIdx = find(diff(sign(momentSeq)) ~= 0,1,'first');
        if isempty(signChangeIdx)
            continue
        end

        m1 = momentSeq(signChangeIdx);
        m2 = momentSeq(signChangeIdx+1);
        if m2 == m1
            continue
        end
        interpFraction = -m1/(m2-m1);

        a1 = aySeq(signChangeIdx);
        a2 = aySeq(signChangeIdx+1);
        r1 = radiusSeq(signChangeIdx);
        r2 = radiusSeq(signChangeIdx+1);

        trimAyG(velocityIndex,steerIndex) = a1 + interpFraction*(a2-a1);
        trimRadius(velocityIndex,steerIndex) = r1 + interpFraction*(r2-r1);
    end
end

%% Step 2: for each target radius, interpolate steer angle and ay per speed
steerAtRadiusDeg = NaN(nVelocity,nRadii);
ayAtRadiusG = NaN(nVelocity,nRadii);

for velocityIndex = 1:nVelocity
    curveRadius = trimRadius(velocityIndex,:);
    curveSteer = SteerAnglesDeg;
    curveAyG = trimAyG(velocityIndex,:);

    validCurve = isfinite(curveRadius) & isfinite(curveAyG);
    if nnz(validCurve) < 2
        continue
    end

    % Radius decreases monotonically as steer increases (typically) --
    % sort by radius ascending for interp1.
    [sortedRadius,sortOrder] = sort(curveRadius(validCurve));
    sortedSteer = curveSteer(validCurve);
    sortedSteer = sortedSteer(sortOrder);
    sortedAyG = curveAyG(validCurve);
    sortedAyG = sortedAyG(sortOrder);

    [sortedRadius,uniqueIdx] = unique(sortedRadius,'stable');
    sortedSteer = sortedSteer(uniqueIdx);
    sortedAyG = sortedAyG(uniqueIdx);

    if numel(sortedRadius) < 2
        continue
    end

    inRange = targetRadii >= min(sortedRadius) & targetRadii <= max(sortedRadius);
    steerAtRadiusDeg(velocityIndex,inRange) = interp1( ...
        sortedRadius,sortedSteer,targetRadii(inRange),'linear');
    ayAtRadiusG(velocityIndex,inRange) = interp1( ...
        sortedRadius,sortedAyG,targetRadii(inRange),'linear');
end

%% Plot: steer angle vs velocity, one line per radius
figure('Name','Constant-Radius Test - Steer vs Speed','Color','k');
hold on;
colors = turbo(nRadii);
for radiusIndex = 1:nRadii
    plot(Velocities,steerAtRadiusDeg(:,radiusIndex),'-o', ...
        'LineWidth',1.6,'MarkerSize',4, ...
        'Color',colors(radiusIndex,:), ...
        'DisplayName',sprintf('R = %.0f m',targetRadii(radiusIndex)));
end
xlabel('Vehicle speed [m/s]');
ylabel('Required handwheel steer angle [deg]');
title('Constant-Radius Test: Steer Angle versus Speed');
legend('Location','best');
grid on;
box on;
fontsize(16,'points')

%% Plot: steer angle vs lateral g, one line per radius (classic SAE plot)
figure('Name','Constant-Radius Test - Steer vs Lateral G','Color','k');
hold on;
for radiusIndex = 1:nRadii
    plot(ayAtRadiusG(:,radiusIndex),steerAtRadiusDeg(:,radiusIndex),'-o', ...
        'LineWidth',1.6,'MarkerSize',4, ...
        'Color',colors(radiusIndex,:), ...
        'DisplayName',sprintf('R = %.0f m',targetRadii(radiusIndex)));
end
xlabel('Lateral acceleration [g]');
ylabel('Required handwheel steer angle [deg]');
title('Constant-Radius Test: Steer Angle versus Lateral Acceleration');
legend('Location','best');
grid on;
box on;
fontsize(16,'points')

%% Table form for quick lookup / export
steerAtRadiusTable = array2table(steerAtRadiusDeg, ...
    'VariableNames',compose('R_%.0fm',targetRadii), ...
    'RowNames',compose('%.0f_mps',Velocities));
disp('Required steer angle [deg] by radius and speed:');
disp(steerAtRadiusTable);
