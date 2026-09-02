%% Calculate understeer gradient using the constant-velocity method
% Run this script after the vehicle YMD sweep script.
%
% Required workspace variables:
%   Velocities, SteerAnglesDeg, BodySlipAnglesDeg
%   ayMap, yawMomentMap, convergedMap, car, g
%
% Method at each fixed velocity:
%   1. Find Mz = 0 equilibrium across body slip for each swept steer angle.
%   2. Calculate the exact kinematic steer for each equilibrium:
%
%          delta_kin = atan(L*ay/V^2)
%
%   3. Calculate excess road-wheel steer:
%
%          delta_excess = delta_actual - delta_kin
%
%   4. Fit:
%
%          delta_excess = K*(ay/g) + intercept
%
%      K is the constant-velocity understeer gradient in deg/g.
%      Positive K means understeer; negative K means oversteer.

%% Settings
angleType = "roadwheel";           % "roadwheel" or "handwheel"
steeringRatio = 1;                 % handwheel angle / road-wheel angle
minFitAyG = 0.05;
maxFitAyG = 0.40;                  % Keep this in the approximately linear range
maxAbsBetaDeg = Inf;
maxRootBetaGapDeg = 0.25;          % Do not bridge large gaps in beta data
minimumPointsPerVelocity = 4;
chooseRootClosestToZeroBeta = true;
makePlots = true;
writeCsvFiles = true;

%% Validate inputs

wheelbase = car.Wheelbase;

%% Find Mz = 0 equilibrium points
pointTemplate = struct( ...
    'velocity',NaN, ...
    'sweepSteerDeg',NaN, ...
    'roadWheelSteerDeg',NaN, ...
    'equilibriumBetaDeg',NaN, ...
    'ay',NaN, ...
    'ayG',NaN, ...
    'radius',NaN, ...
    'kinematicSteerDeg',NaN, ...
    'steerExcessDeg',NaN);
equilibriumPoints = repmat(pointTemplate,0,1);

for velocityIndex = 1:numel(Velocities)
    velocity = Velocities(velocityIndex);

    for steerIndex = 1:numel(SteerAnglesDeg)
        beta = BodySlipAnglesDeg(:);
        moment = squeeze(yawMomentMap(velocityIndex,steerIndex,:));
        ay = squeeze(ayMap(velocityIndex,steerIndex,:));
        valid = squeeze(convergedMap(velocityIndex,steerIndex,:));

        valid = valid & isfinite(beta) & isfinite(moment) & isfinite(ay);
        beta = beta(valid);
        moment = moment(valid);
        ay = ay(valid);

        if numel(beta) < 2
            continue
        end

        [beta,order] = sort(beta);
        moment = moment(order);
        ay = ay(order);
        roots = findMomentRoots(beta,moment,ay,maxRootBetaGapDeg);

        if isempty(roots)
            continue
        end

        if chooseRootClosestToZeroBeta
            [~,rootIndex] = min(abs([roots.betaDeg]));
            roots = roots(rootIndex);
        end

        for rootIndex = 1:numel(roots)
            ayEq = roots(rootIndex).ay;

            % This script analyzes the positive-turn branch. Use mirrored
            % data or modify this condition for negative-turn analysis.
            if ~isfinite(ayEq) || ayEq <= 0
                continue
            end

            if angleType == "handwheel"
                roadWheelSteerDeg = SteerAnglesDeg(steerIndex)/steeringRatio;
            else
                roadWheelSteerDeg = SteerAnglesDeg(steerIndex);
            end

            radiusEq = velocity^2/ayEq;
            kinematicSteerDeg = rad2deg(atan(wheelbase/radiusEq));

            point = pointTemplate;
            point.velocity = velocity;
            point.sweepSteerDeg = SteerAnglesDeg(steerIndex);
            point.roadWheelSteerDeg = roadWheelSteerDeg;
            point.equilibriumBetaDeg = roots(rootIndex).betaDeg;
            point.ay = ayEq;
            point.ayG = ayEq/g;
            point.radius = radiusEq;
            point.kinematicSteerDeg = kinematicSteerDeg;
            point.steerExcessDeg = roadWheelSteerDeg-kinematicSteerDeg;
            equilibriumPoints(end+1,1) = point; %#ok<SAGROW>
        end
    end
end

if isempty(equilibriumPoints)
    error(['No positive-lateral-acceleration Mz = 0 equilibria were found. ', ...
        'Check the yaw-moment sign convention and beta sweep range.']);
end

equilibriumTable = struct2table(equilibriumPoints);
equilibriumTable = sortrows(equilibriumTable,{'velocity','ayG'});

%% Select the approximately linear operating range
fitMask = isfinite(equilibriumTable.ayG) & ...
    isfinite(equilibriumTable.roadWheelSteerDeg) & ...
    isfinite(equilibriumTable.steerExcessDeg) & ...
    isfinite(equilibriumTable.equilibriumBetaDeg) & ...
    equilibriumTable.ayG >= minFitAyG & ...
    equilibriumTable.ayG <= maxFitAyG & ...
    abs(equilibriumTable.equilibriumBetaDeg) <= maxAbsBetaDeg;
fitPoints = equilibriumTable(fitMask,:);

if isempty(fitPoints)
    error(['No equilibrium points passed the fit filters. Adjust minFitAyG, ', ...
        'maxFitAyG, maxAbsBetaDeg, or the original sweep range.']);
end

%% Fit one understeer gradient at each constant velocity
fitTemplate = struct( ...
    'velocity',NaN, ...
    'pointCount',0, ...
    'minAyG',NaN, ...
    'maxAyG',NaN, ...
    'KdegPerG',NaN, ...
    'interceptDeg',NaN, ...
    'rSquared',NaN, ...
    'KZeroInterceptDegPerG',NaN, ...
    'rawSteerSlopeDegPerG',NaN, ...
    'geometricSlopeDegPerG',NaN, ...
    'KSmallAngleDegPerG',NaN, ...
    'handlingBalance',"");
velocityFits = repmat(fitTemplate,0,1);

for velocityIndex = 1:numel(Velocities)
    velocity = Velocities(velocityIndex);
    points = fitPoints(fitPoints.velocity == velocity,:);

    if height(points) < minimumPointsPerVelocity
        continue
    end

    % Sort and average duplicate ay/g values before fitting.
    [xSorted,order] = sort(points.ayG);
    excessSorted = points.steerExcessDeg(order);
    steerSorted = points.roadWheelSteerDeg(order);
    [x,~,group] = unique(xSorted);
    excessSteer = accumarray(group,excessSorted,[],@mean);
    actualSteer = accumarray(group,steerSorted,[],@mean);

    if numel(x) < minimumPointsPerVelocity || range(x) <= eps
        continue
    end

    % Primary fit: free intercept. This prevents a steering/toe offset from
    % biasing the slope used as the understeer gradient.
    freeFit = polyfit(x,excessSteer,1);
    KdegPerG = freeFit(1);
    interceptDeg = freeFit(2);
    predicted = polyval(freeFit,x);
    residual = excessSteer-predicted;
    ssResidual = sum(residual.^2);
    ssTotal = sum((excessSteer-mean(excessSteer)).^2);
    if ssTotal > eps
        rSquared = 1-ssResidual/ssTotal;
    else
        rSquared = NaN;
    end

    % Diagnostic fit forced through zero.
    KZeroInterceptDegPerG = (x.'*excessSteer)/(x.'*x);

    % Small-angle cross-check:
    % slope(delta versus ay/g) = K + rad2deg(g*L/V^2).
    rawFit = polyfit(x,actualSteer,1);
    rawSteerSlopeDegPerG = rawFit(1);
    geometricSlopeDegPerG = rad2deg(g*wheelbase/velocity^2);
    KSmallAngleDegPerG = rawSteerSlopeDegPerG-geometricSlopeDegPerG;

    if KdegPerG > 1e-6
        handlingBalance = "understeer";
    elseif KdegPerG < -1e-6
        handlingBalance = "oversteer";
    else
        handlingBalance = "neutral";
    end

    fit = fitTemplate;
    fit.velocity = velocity;
    fit.pointCount = numel(x);
    fit.minAyG = min(x);
    fit.maxAyG = max(x);
    fit.KdegPerG = KdegPerG;
    fit.interceptDeg = interceptDeg;
    fit.rSquared = rSquared;
    fit.KZeroInterceptDegPerG = KZeroInterceptDegPerG;
    fit.rawSteerSlopeDegPerG = rawSteerSlopeDegPerG;
    fit.geometricSlopeDegPerG = geometricSlopeDegPerG;
    fit.KSmallAngleDegPerG = KSmallAngleDegPerG;
    fit.handlingBalance = handlingBalance;
    velocityFits(end+1,1) = fit; %#ok<SAGROW>
end

if isempty(velocityFits)
    error(['No velocity had enough valid points for a fit. Reduce ', ...
        'minimumPointsPerVelocity or expand the original steer sweep.']);
end

perVelocityTable = struct2table(velocityFits);

%% Overall summary across the valid constant-velocity fits
validK = isfinite(perVelocityTable.KdegPerG);
weights = perVelocityTable.pointCount(validK);
Kvalues = perVelocityTable.KdegPerG(validK);
overallKdegPerG = sum(weights.*Kvalues)/sum(weights);
medianKdegPerG = median(Kvalues);
minimumKdegPerG = min(Kvalues);
maximumKdegPerG = max(Kvalues);

if overallKdegPerG > 1e-6
    overallBalance = "understeer";
elseif overallKdegPerG < -1e-6
    overallBalance = "oversteer";
else
    overallBalance = "neutral";
end

summaryTable = table( ...
    overallKdegPerG,medianKdegPerG,minimumKdegPerG,maximumKdegPerG, ...
    string(overallBalance),minFitAyG,maxFitAyG,height(fitPoints), ...
    'VariableNames',{ ...
    'WeightedMeanKdegPerG','MedianKdegPerG','MinimumKdegPerG', ...
    'MaximumKdegPerG','HandlingBalance','MinimumFitAyG', ...
    'MaximumFitAyG','TotalFitPointCount'});

fprintf('\nConstant-velocity understeer results\n');
fprintf('------------------------------------\n');
fprintf('Weighted mean K: %+.4f deg/g\n',overallKdegPerG);
fprintf('Median K:        %+.4f deg/g\n',medianKdegPerG);
fprintf('Range:           %+.4f to %+.4f deg/g\n', ...
    minimumKdegPerG,maximumKdegPerG);
fprintf('Classification:  %s\n\n',overallBalance);
disp(perVelocityTable(:,{'velocity','pointCount','KdegPerG', ...
    'interceptDeg','rSquared','handlingBalance'}));

%% Export tables
if writeCsvFiles
    writetable(equilibriumTable,'understeer_equilibrium_points.csv');
    writetable(fitPoints,'understeer_constant_velocity_fit_points.csv');
    writetable(perVelocityTable,'understeer_by_velocity.csv');
    writetable(summaryTable,'understeer_summary.csv');
end

%% Plots
if makePlots
    figure('Name','Constant-Velocity Understeer Fits');
    tiledlayout('flow');

    for fitIndex = 1:height(perVelocityTable)
        velocity = perVelocityTable.velocity(fitIndex);
        points = fitPoints(fitPoints.velocity == velocity,:);
        if isempty(points)
            continue
        end

        nexttile;
        scatter(points.ayG,points.steerExcessDeg,28,'filled');
        hold on;
        xLine = linspace(min(points.ayG),max(points.ayG),100).';
        yLine = perVelocityTable.KdegPerG(fitIndex)*xLine + ...
            perVelocityTable.interceptDeg(fitIndex);
        plot(xLine,yLine,'LineWidth',1.5);
        yline(0,'k:');
        grid on;
        xlabel('Lateral acceleration, a_y/g');
        ylabel('Excess road-wheel steer [deg]');
        title(sprintf('V = %.1f m/s, K = %+.3f deg/g', ...
            velocity,perVelocityTable.KdegPerG(fitIndex)));
    end

    figure('Name','Understeer Gradient versus Velocity');
    plot(perVelocityTable.velocity,perVelocityTable.KdegPerG, ...
        'o-','LineWidth',1.5,'MarkerFaceColor','auto');
    yline(0,'k:');
    grid on;
    xlabel('Velocity [m/s]');
    ylabel('Understeer gradient [deg/g]');
    title(sprintf('Weighted mean K = %+.3f deg/g (%s)', ...
        overallKdegPerG,overallBalance));
    fontsize(18,'points')
end

%% Local function: interpolate zero crossings of yaw moment
function roots = findMomentRoots(beta,moment,ay,maxBetaGapDeg)
rootTemplate = struct('betaDeg',NaN,'ay',NaN);
roots = repmat(rootTemplate,0,1);

for index = 1:numel(beta)-1
    beta1 = beta(index);
    beta2 = beta(index+1);
    moment1 = moment(index);
    moment2 = moment(index+1);
    ay1 = ay(index);
    ay2 = ay(index+1);

    if beta2-beta1 > maxBetaGapDeg
        continue
    end

    if moment1 == 0
        root = rootTemplate;
        root.betaDeg = beta1;
        root.ay = ay1;
        roots(end+1,1) = root; %#ok<AGROW>
    end

    if moment1*moment2 < 0
        fraction = -moment1/(moment2-moment1);
        root = rootTemplate;
        root.betaDeg = beta1+fraction*(beta2-beta1);
        root.ay = ay1+fraction*(ay2-ay1);
        roots(end+1,1) = root; %#ok<AGROW>
    end
end

if moment(end) == 0
    root = rootTemplate;
    root.betaDeg = beta(end);
    root.ay = ay(end);
    roots(end+1,1) = root;
end

% Remove duplicate roots caused by an exact zero shared by two intervals.
if numel(roots) > 1
    betaRoots = [roots.betaDeg].';
    [~,uniqueIndices] = unique(round(betaRoots,10),'stable');
    roots = roots(uniqueIndices);
end
end
