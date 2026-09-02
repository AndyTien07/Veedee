%% Curvature-compensated understeer-gradient analysis
% Removes geometric steering sample-by-sample:
%
%   delta_excess = delta_roadwheel - atan(L*yawRate/vehicleSpeed)
%
% Then fits, separately for each turn direction:
%
%   signed(delta_excess) = K*abs(ay/g) + intercept
%
% K is the understeer gradient in deg/g. Positive K indicates understeer.
% This method does not require a constant-radius segment.
%
% Outlier rejection is performed in two stages:
%   1. Samples within each lateral-g bin are rejected using a robust MAD test.
%   2. Bin means are iteratively rejected using residuals from the fitted line.

clearvars;
close all;
clc;

%% User settings
DATA_PATH = "T36 Comp Skidpad.csv";

wheelbaseM = 1.545;
steeringRatio = 1;          % steering-wheel angle / road-wheel angle
g = 9.80665;
rawYawRateUnits = "deg/s";    % "rad/s" or "deg/s"

% Change these only if a sensor has the opposite sign convention.
steerSign = -1;
yawRateSign = 1;
lateralAccelSign = 1;
longitudinalAccelSign = 1;

filterCutoffHz = 6;
startTimeS = 30;
endTimeS = 50;

% Quasi-steady sample filters.
minimumSpeedMS = 5;
minimumAbsLateralG = 0.6;
maximumAbsLateralG = 2.50;
maximumAbsLongitudinalG = 0.25;
maximumAbsYawAccelDegS2 = 50;
maximumAbsLateralJerkMS3 = 5;
minimumRadiusM = 8;
maximumRadiusM = 12;
minimumAbsYawRateDegS = 3;
requireAyYawRateSignAgreement = true;

% Bin settings. Binning prevents long sections from dominating the fit.
lateralGBinWidth = 0.1;
minimumSamplesPerBin = 1;
minimumBinsPerTurn = 5;
minimumFitGRange = 0.20;

% Robust outlier rejection.
enableOutlierRejection = true;
withinBinOutlierSigma = 3.5;   % MAD threshold for raw samples inside each bin
fitOutlierSigma = 3.5;         % MAD threshold for residuals of binned fit
maximumOutlierIterations = 10;
minimumRobustScaleDeg = 0.01;  % prevents rejection from a nearly zero MAD

makePlots = true;
writeCsvFiles = true;

%% Import and validate data
opts = detectImportOptions(DATA_PATH);
data = readtable(DATA_PATH,opts);

time = double(data.Time(:));

dt = median(diff(time));
sampleRateHz = 1/dt;
nyquistHz = sampleRateHz/2;
if filterCutoffHz <= 0 || filterCutoffHz >= nyquistHz
    error('filterCutoffHz must be between 0 and %.3f Hz.',nyquistHz);
end

[bFilter,aFilter] = butter(2,filterCutoffHz/nyquistHz,'low');

%% Filter measured channels
steeringWheelDeg = steerSign*filterSignal(double(data.Steer(:)),bFilter,aFilter);
steeringWheelDeg=steeringWheelDeg-steeringWheelDeg(100)
roadWheelSteerDeg = steeringWheelDeg/steeringRatio;

ax = longitudinalAccelSign*filterSignal(double(data.Ax(:)),bFilter,aFilter);
ay = lateralAccelSign*filterSignal(double(data.Ay(:)),bFilter,aFilter);

rawYawRate = yawRateSign*filterSignal(double(data.YawR(:)),bFilter,aFilter);

switch rawYawRateUnits
    case "rad/s"
        yawRateRadS = rawYawRate;
    case "deg/s"
        yawRateRadS = deg2rad(rawYawRate);
    otherwise
        error('rawYawRateUnits must be "rad/s" or "deg/s".');
end
yawRateDegS = rad2deg(yawRateRadS);

%% Wheel-speed-based vehicle speed
frontConversion = 1.24482467/(13.55*60);
rearConversion = 1.24482467/(13.55*60);

frontRightSpeedMS = -double(data.FRrpm(:))*frontConversion;
frontLeftSpeedMS  =  double(data.FLrpm(:))*frontConversion;
rearRightSpeedMS  = -double(data.RRrpm(:))*rearConversion;
rearLeftSpeedMS   =  double(data.RLrpm(:))*rearConversion;

% Average the non-driven/front wheels, matching the original analysis.
vehicleSpeedRawMS = abs((frontRightSpeedMS+frontLeftSpeedMS)/2);
vehicleSpeedMS = filterSignal(vehicleSpeedRawMS,bFilter,aFilter);

%% Derived channels
yawAccelDegS2 = gradient(yawRateDegS)./gradient(time);
lateralJerkMS3 = gradient(ay)./gradient(time);
lateralG = ay/g;
longitudinalG = ax/g;

% Signed path curvature and unsigned radius at the vehicle CG.
curvaturePerM = yawRateRadS./vehicleSpeedMS;
radiusM = 1./abs(curvaturePerM);

% Exact bicycle-model geometric steer for measured curvature.
geometricSteerDeg = rad2deg(atan(wheelbaseM*curvaturePerM));
excessSteerDeg = roadWheelSteerDeg-geometricSteerDeg;

% Normalize both turn directions to positive lateral-g coordinates.
turnSign = sign(yawRateRadS);
normalizedLateralG = turnSign.*lateralG;
normalizedExcessSteerDeg = turnSign.*excessSteerDeg;
normalizedRoadWheelSteerDeg = turnSign.*roadWheelSteerDeg;
normalizedGeometricSteerDeg = turnSign.*geometricSteerDeg;

%% Select quasi-steady cornering samples
timeMask = time >= startTimeS & time <= endTimeS;
finiteMask = isfinite(time) & isfinite(vehicleSpeedMS) & ...
    isfinite(roadWheelSteerDeg) & isfinite(yawRateRadS) & ...
    isfinite(yawAccelDegS2) & isfinite(lateralJerkMS3) & ...
    isfinite(lateralG) & isfinite(longitudinalG) & ...
    isfinite(radiusM) & isfinite(excessSteerDeg);

steadyMask = timeMask & finiteMask & ...
    vehicleSpeedMS >= minimumSpeedMS & ...
    abs(lateralG) >= minimumAbsLateralG & ...
    abs(lateralG) <= maximumAbsLateralG & ...
    abs(longitudinalG) <= maximumAbsLongitudinalG & ...
    abs(yawAccelDegS2) <= maximumAbsYawAccelDegS2 & ...
    abs(lateralJerkMS3) <= maximumAbsLateralJerkMS3 & ...
    radiusM >= minimumRadiusM & radiusM <= maximumRadiusM & ...
    abs(yawRateDegS) >= minimumAbsYawRateDegS;

if requireAyYawRateSignAgreement
    steadyMask = steadyMask & normalizedLateralG > 0;
end

leftMask = steadyMask & turnSign > 0;
rightMask = steadyMask & turnSign < 0;

if nnz(steadyMask) == 0
    error(['No samples passed the filters. Check sensor signs, yaw-rate units, ', ...
        'time limits, and quasi-steady thresholds.']);
end

%% Bin each turn direction
leftBins = makeLateralGBins(normalizedLateralG(leftMask), ...
    normalizedExcessSteerDeg(leftMask),lateralGBinWidth, ...
    minimumSamplesPerBin,enableOutlierRejection,withinBinOutlierSigma);
rightBins = makeLateralGBins(normalizedLateralG(rightMask), ...
    normalizedExcessSteerDeg(rightMask),lateralGBinWidth, ...
    minimumSamplesPerBin,enableOutlierRejection,withinBinOutlierSigma);

%% Fit separate gradients with iterative residual rejection
leftFit = fitTurn(leftBins,"positive-yaw",minimumBinsPerTurn, ...
    minimumFitGRange,enableOutlierRejection,fitOutlierSigma, ...
    maximumOutlierIterations,minimumRobustScaleDeg);
rightFit = fitTurn(rightBins,"negative-yaw",minimumBinsPerTurn, ...
    minimumFitGRange,enableOutlierRejection,fitOutlierSigma, ...
    maximumOutlierIterations,minimumRobustScaleDeg);

% Apply fit classifications back to exported bin tables.
if height(leftBins) > 0
    leftBins.FitInlier = leftFit.inlierMask;
    leftBins.FitResidualDeg = leftFit.residualAllDeg;
end
if height(rightBins) > 0
    rightBins.FitInlier = rightFit.inlierMask;
    rightBins.FitResidualDeg = rightFit.residualAllDeg;
end

%% Common-gradient fit with separate turn-direction intercepts
% This estimates one K while allowing steering zero, asymmetry, or banking to
% produce different positive-yaw and negative-yaw intercepts.
commonFit = fitCommonGradient(leftBins,rightBins,minimumBinsPerTurn, ...
    minimumFitGRange,enableOutlierRejection,fitOutlierSigma, ...
    maximumOutlierIterations,minimumRobustScaleDeg);

if commonFit.valid
    reportedKDegPerG = commonFit.KDegPerG;
elseif leftFit.valid && rightFit.valid
    reportedKDegPerG = mean([leftFit.KDegPerG,rightFit.KDegPerG]);
elseif leftFit.valid
    reportedKDegPerG = leftFit.KDegPerG;
elseif rightFit.valid
    reportedKDegPerG = rightFit.KDegPerG;
else
    error(['Neither turn direction contains enough independent lateral-g ', ...
        'bins for a reliable fit. Relax the filters or use more test data.']);
end

if reportedKDegPerG > 0
    handlingBalance = "understeer";
    characteristicVelocityMS = sqrt(g*wheelbaseM/deg2rad(reportedKDegPerG));
    characteristicVelocityKPH = 3.6*characteristicVelocityMS;
elseif reportedKDegPerG < 0
    handlingBalance = "oversteer";
    characteristicVelocityMS = NaN;
    characteristicVelocityKPH = NaN;
else
    handlingBalance = "neutral steer";
    characteristicVelocityMS = Inf;
    characteristicVelocityKPH = Inf;
end

%% Display results
fprintf('\nCurvature-compensated understeer analysis\n');
fprintf('-----------------------------------------\n');
fprintf('Selected raw samples:             %d\n',nnz(steadyMask));
fprintf('Positive-yaw samples/bins:        %d / %d\n',nnz(leftMask),height(leftBins));
fprintf('Negative-yaw samples/bins:        %d / %d\n',nnz(rightMask),height(rightBins));
fprintf('Outlier rejection enabled:        %s\n',string(enableOutlierRejection));

printTurnFit(leftFit);
printTurnFit(rightFit);

if commonFit.valid
    fprintf('\nCommon-gradient fit\n');
    fprintf('Common K:                        %9.5f deg/g\n',commonFit.KDegPerG);
    fprintf('Common K 95%% CI:                 %9.5f to %9.5f deg/g\n', ...
        commonFit.KCILowDegPerG,commonFit.KCIHighDegPerG);
    fprintf('Positive-yaw intercept:          %9.5f deg\n',commonFit.leftInterceptDeg);
    fprintf('Negative-yaw intercept:          %9.5f deg\n',commonFit.rightInterceptDeg);
    fprintf('R-squared / RMSE:                %9.5f / %.5f deg\n', ...
        commonFit.rSquared,commonFit.rmseDeg);
    fprintf('Used/rejected bins:              %d / %d\n', ...
        commonFit.usedBinCount,commonFit.rejectedBinCount);
else
    fprintf('\nCommon-gradient fit unavailable: %s\n',commonFit.message);
end

fprintf('\nReported K:                      %9.5f deg/g\n',reportedKDegPerG);
fprintf('Handling balance:                %s\n',handlingBalance);
if isfinite(characteristicVelocityMS)
    fprintf('Characteristic velocity:         %.3f m/s (%.3f km/h)\n', ...
        characteristicVelocityMS,characteristicVelocityKPH);
elseif isinf(characteristicVelocityMS)
    fprintf('Characteristic velocity:         Inf (neutral steer)\n');
else
    fprintf('Characteristic velocity:         N/A for oversteer\n');
end

%% Export results
selectedSamples = table(time(steadyMask),vehicleSpeedMS(steadyMask), ...
    lateralG(steadyMask),longitudinalG(steadyMask),yawRateDegS(steadyMask), ...
    radiusM(steadyMask),roadWheelSteerDeg(steadyMask), ...
    geometricSteerDeg(steadyMask),excessSteerDeg(steadyMask), ...
    normalizedLateralG(steadyMask),normalizedExcessSteerDeg(steadyMask), ...
    'VariableNames',{'TimeS','SpeedMS','LateralG','LongitudinalG', ...
    'YawRateDegS','RadiusM','RoadWheelSteerDeg','GeometricSteerDeg', ...
    'ExcessSteerDeg','NormalizedLateralG','NormalizedExcessSteerDeg'});

summaryTable = table(reportedKDegPerG,string(handlingBalance), ...
    characteristicVelocityMS,characteristicVelocityKPH,nnz(steadyMask), ...
    leftFit.KDegPerG,rightFit.KDegPerG,commonFit.KDegPerG, ...
    leftFit.usedBinCount,rightFit.usedBinCount,commonFit.usedBinCount, ...
    leftFit.rejectedBinCount,rightFit.rejectedBinCount,commonFit.rejectedBinCount, ...
    'VariableNames',{'ReportedKDegPerG','HandlingBalance', ...
    'CharacteristicVelocityMS','CharacteristicVelocityKPH','SelectedSampleCount', ...
    'PositiveYawKDegPerG','NegativeYawKDegPerG','CommonKDegPerG', ...
    'PositiveYawUsedBins','NegativeYawUsedBins','CommonUsedBins', ...
    'PositiveYawRejectedBins','NegativeYawRejectedBins','CommonRejectedBins'});

if writeCsvFiles
    writetable(selectedSamples,'curvature_compensated_selected_samples.csv');
    writetable(leftBins,'curvature_compensated_positive_yaw_bins.csv');
    writetable(rightBins,'curvature_compensated_negative_yaw_bins.csv');
    writetable(summaryTable,'curvature_compensated_summary.csv');
end

%% Plots
figure()
hold on
plot(time,steeringWheelDeg)
plot(time,yawRateDegS/5)
plot(time,ay)
if makePlots
    figure('Name','Curvature-Compensated Understeer Fit','Color','w');
    tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    nexttile;
    hold on;
    grid on;

    scatter(normalizedLateralG(leftMask),normalizedExcessSteerDeg(leftMask), ...
        8,[0.3 0.6 1.0],'filled','MarkerFaceAlpha',0.12, ...
        'DisplayName','Positive-yaw samples');
    scatter(normalizedLateralG(rightMask),normalizedExcessSteerDeg(rightMask), ...
        8,[1.0 0.5 0.3],'filled','MarkerFaceAlpha',0.12, ...
        'DisplayName','Negative-yaw samples');

    plotTurnBins(leftBins,[0 0.35 0.9],'Positive-yaw bins');
    plotTurnBins(rightBins,[0.9 0.25 0.05],'Negative-yaw bins');

    if leftFit.valid
        xLine = linspace(min(leftBins.MeanLateralG(leftFit.inlierMask)), ...
            max(leftBins.MeanLateralG(leftFit.inlierMask)),100);
        plot(xLine,leftFit.KDegPerG*xLine+leftFit.interceptDeg, ...
            '-','Color',[0 0.25 0.8],'LineWidth',2, ...
            'DisplayName',sprintf('Positive yaw: K = %.3f deg/g',leftFit.KDegPerG));
    end
    if rightFit.valid
        xLine = linspace(min(rightBins.MeanLateralG(rightFit.inlierMask)), ...
            max(rightBins.MeanLateralG(rightFit.inlierMask)),100);
        plot(xLine,rightFit.KDegPerG*xLine+rightFit.interceptDeg, ...
            '-','Color',[0.8 0.15 0],'LineWidth',2, ...
            'DisplayName',sprintf('Negative yaw: K = %.3f deg/g',rightFit.KDegPerG));
    end

    xlabel('Normalized lateral acceleration (g)');
    ylabel('Normalized excess steer (deg)');
    title('Curvature-compensated data and robust fits');
    legend('Location','best');

    nexttile;
    hold on;
    grid on;
    plot(time,roadWheelSteerDeg,'DisplayName','Road-wheel steer');
    plot(time,geometricSteerDeg,'DisplayName','Geometric steer');
    plot(time,excessSteerDeg,'DisplayName','Excess steer');
    scatter(time(steadyMask),excessSteerDeg(steadyMask),8,'filled', ...
        'DisplayName','Selected samples');
    xlabel('Time (s)');
    ylabel('Angle (deg)');
    title('Steering components and selected samples');
    legend('Location','best');
end

%% Local functions
function filtered = filterSignal(signal,b,a)
    signal = double(signal(:));
    finite = isfinite(signal);
    if nnz(finite) < 2
        filtered = NaN(size(signal));
        return
    end
    signal = fillmissing(signal,'linear','EndValues','nearest');
    filtered = filtfilt(b,a,signal);
end

function bins = makeLateralGBins(x,y,binWidth,minSamples,doReject,sigmaLimit)
    x = x(:);
    y = y(:);
    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    emptyBins = table([],[],[],[],[],[],[],[], ...
        'VariableNames',{'BinCenterG','MeanLateralG','MeanExcessSteerDeg', ...
        'StdExcessSteerDeg','RawSampleCount','UsedSampleCount', ...
        'WithinBinRejectedCount','WithinBinRobustScaleDeg'});
    if isempty(x)
        bins = emptyBins;
        return
    end

    binIndex = floor(x/binWidth);
    groups = unique(binIndex);
    rows = zeros(0,8);

    for groupIndex = 1:numel(groups)
        member = binIndex == groups(groupIndex);
        xb = x(member);
        yb = y(member);
        rawCount = numel(yb);
        if rawCount < minSamples
            continue
        end

        keep = true(rawCount,1);
        robustScale = robustMADScale(yb);
        if doReject && rawCount >= 4
            center = median(yb);
            threshold = sigmaLimit*max(robustScale,eps);
            keep = abs(yb-center) <= threshold;
        end

        if nnz(keep) < minSamples
            continue
        end

        rows(end+1,:) = [ ... %#ok<AGROW>
            (groups(groupIndex)+0.5)*binWidth, ...
            mean(xb(keep)), ...
            mean(yb(keep)), ...
            std(yb(keep),0), ...
            rawCount, ...
            nnz(keep), ...
            rawCount-nnz(keep), ...
            robustScale];
    end

    if isempty(rows)
        bins = emptyBins;
    else
        bins = array2table(rows,'VariableNames',emptyBins.Properties.VariableNames);
        bins = sortrows(bins,'MeanLateralG');
    end
end

function fit = fitTurn(bins,label,minBins,minGRange,doReject,sigmaLimit,maxIter,minScale)
    fit = emptyTurnFit(label,height(bins));
    if height(bins) < minBins
        fit.message = "Too few lateral-g bins.";
        return
    end

    x = bins.MeanLateralG;
    y = bins.MeanExcessSteerDeg;
    baseValid = isfinite(x) & isfinite(y);
    inlier = baseValid;

    for iteration = 1:maxIter
        if nnz(inlier) < minBins || range(x(inlier)) < minGRange
            break
        end

        coefficients = polyfit(x(inlier),y(inlier),1);
        residualAll = y-polyval(coefficients,x);
        if ~doReject
            break
        end

        centerResidual = median(residualAll(inlier));
        scale = max(robustMADScale(residualAll(inlier)),minScale);
        newInlier = baseValid & ...
            abs(residualAll-centerResidual) <= sigmaLimit*scale;

        if isequal(newInlier,inlier)
            break
        end
        if nnz(newInlier) < minBins || range(x(newInlier)) < minGRange
            break
        end
        inlier = newInlier;
    end

    if nnz(inlier) < minBins
        fit.message = "Too few bins remained after outlier rejection.";
        fit.inlierMask = inlier;
        return
    end
    if range(x(inlier)) < minGRange
        fit.message = "Insufficient lateral-g range after outlier rejection.";
        fit.inlierMask = inlier;
        return
    end

    [coefficients,S] = polyfit(x(inlier),y(inlier),1);
    predictedInlier = polyval(coefficients,x(inlier));
    residualInlier = y(inlier)-predictedInlier;
    residualAll = y-polyval(coefficients,x);

    n = nnz(inlier);
    dof = n-2;
    ssResidual = sum(residualInlier.^2);
    ssTotal = sum((y(inlier)-mean(y(inlier))).^2);
    if ssTotal > eps
        rSquared = 1-ssResidual/ssTotal;
    else
        rSquared = NaN;
    end
    rmse = sqrt(mean(residualInlier.^2));

    if dof > 0 && isfield(S,'R') && rcond(S.R) > eps
        covariance = (ssResidual/dof)*inv(S.R)*inv(S.R)'; %#ok<MINV>
        slopeSE = sqrt(max(covariance(1,1),0));
        tCritical = studentTCritical95(dof);
        ci = coefficients(1)+[-1 1]*tCritical*slopeSE;
    else
        slopeSE = NaN;
        ci = [NaN NaN];
    end

    fit.valid = true;
    fit.message = "OK";
    fit.KDegPerG = coefficients(1);
    fit.interceptDeg = coefficients(2);
    fit.KStandardErrorDegPerG = slopeSE;
    fit.KCILowDegPerG = ci(1);
    fit.KCIHighDegPerG = ci(2);
    fit.rSquared = rSquared;
    fit.rmseDeg = rmse;
    fit.usedBinCount = n;
    fit.rejectedBinCount = nnz(baseValid)-n;
    fit.inlierMask = inlier;
    fit.residualAllDeg = residualAll;
end

function fit = fitCommonGradient(leftBins,rightBins,minBins,minGRange, ...
        doReject,sigmaLimit,maxIter,minScale)
    nLeft = height(leftBins);
    nRight = height(rightBins);
    fit = struct('valid',false,'message',"",'KDegPerG',NaN, ...
        'leftInterceptDeg',NaN,'rightInterceptDeg',NaN, ...
        'KStandardErrorDegPerG',NaN,'KCILowDegPerG',NaN, ...
        'KCIHighDegPerG',NaN,'rSquared',NaN,'rmseDeg',NaN, ...
        'usedBinCount',0,'rejectedBinCount',0, ...
        'leftInlierMask',false(nLeft,1),'rightInlierMask',false(nRight,1));

    if nLeft < minBins || nRight < minBins
        fit.message = "Both turn directions need enough bins for the common fit.";
        return
    end

    x = [leftBins.MeanLateralG;rightBins.MeanLateralG];
    y = [leftBins.MeanExcessSteerDeg;rightBins.MeanExcessSteerDeg];
    isLeft = [true(nLeft,1);false(nRight,1)];
    isRight = ~isLeft;
    baseValid = isfinite(x) & isfinite(y);
    inlier = baseValid;

    for iteration = 1:maxIter
        if nnz(inlier & isLeft) < minBins || nnz(inlier & isRight) < minBins
            break
        end
        if range(x(inlier & isLeft)) < minGRange || ...
                range(x(inlier & isRight)) < minGRange
            break
        end

        design = [x double(isLeft) double(isRight)];
        coefficients = design(inlier,:)\y(inlier);
        residualAll = y-design*coefficients;
        if ~doReject
            break
        end

        % Center residuals independently by direction so one branch cannot
        % cause rejection in the other branch.
        centeredResidual = residualAll;
        centeredResidual(isLeft) = residualAll(isLeft)-median(residualAll(inlier & isLeft));
        centeredResidual(isRight) = residualAll(isRight)-median(residualAll(inlier & isRight));
        scale = max(robustMADScale(centeredResidual(inlier)),minScale);
        newInlier = baseValid & abs(centeredResidual) <= sigmaLimit*scale;

        if isequal(newInlier,inlier)
            break
        end
        if nnz(newInlier & isLeft) < minBins || nnz(newInlier & isRight) < minBins
            break
        end
        if range(x(newInlier & isLeft)) < minGRange || ...
                range(x(newInlier & isRight)) < minGRange
            break
        end
        inlier = newInlier;
    end

    if nnz(inlier & isLeft) < minBins || nnz(inlier & isRight) < minBins
        fit.message = "Too few bins remained in one turn direction.";
        return
    end
    if range(x(inlier & isLeft)) < minGRange || ...
            range(x(inlier & isRight)) < minGRange
        fit.message = "Insufficient lateral-g range in one turn direction.";
        return
    end

    design = [x double(isLeft) double(isRight)];
    X = design(inlier,:);
    coefficients = X\y(inlier);
    predicted = X*coefficients;
    residual = y(inlier)-predicted;
    residualAll = y-design*coefficients;

    n = nnz(inlier);
    parameterCount = 3;
    dof = n-parameterCount;
    ssResidual = sum(residual.^2);
    ssTotal = sum((y(inlier)-mean(y(inlier))).^2);
    if ssTotal > eps
        rSquared = 1-ssResidual/ssTotal;
    else
        rSquared = NaN;
    end

    if dof > 0 && rcond(X'*X) > eps
        covariance = (ssResidual/dof)*inv(X'*X); %#ok<MINV>
        slopeSE = sqrt(max(covariance(1,1),0));
        tCritical = studentTCritical95(dof);
        ci = coefficients(1)+[-1 1]*tCritical*slopeSE;
    else
        slopeSE = NaN;
        ci = [NaN NaN];
    end

    fit.valid = true;
    fit.message = "OK";
    fit.KDegPerG = coefficients(1);
    fit.leftInterceptDeg = coefficients(2);
    fit.rightInterceptDeg = coefficients(3);
    fit.KStandardErrorDegPerG = slopeSE;
    fit.KCILowDegPerG = ci(1);
    fit.KCIHighDegPerG = ci(2);
    fit.rSquared = rSquared;
    fit.rmseDeg = sqrt(mean(residual.^2));
    fit.usedBinCount = n;
    fit.rejectedBinCount = nnz(baseValid)-n;
    fit.leftInlierMask = inlier(1:nLeft);
    fit.rightInlierMask = inlier(nLeft+1:end);
    fit.residualAllDeg = residualAll;
end

function fit = emptyTurnFit(label,nBins)
    fit = struct('label',label,'valid',false,'message',"", ...
        'KDegPerG',NaN,'interceptDeg',NaN, ...
        'KStandardErrorDegPerG',NaN,'KCILowDegPerG',NaN, ...
        'KCIHighDegPerG',NaN,'rSquared',NaN,'rmseDeg',NaN, ...
        'usedBinCount',0,'rejectedBinCount',0, ...
        'inlierMask',false(nBins,1),'residualAllDeg',NaN(nBins,1));
end

function printTurnFit(fit)
    fprintf('\n%s fit\n',fit.label);
    if ~fit.valid
        fprintf('Unavailable:                      %s\n',fit.message);
        return
    end
    fprintf('K:                                %9.5f deg/g\n',fit.KDegPerG);
    fprintf('K 95%% CI:                         %9.5f to %9.5f deg/g\n', ...
        fit.KCILowDegPerG,fit.KCIHighDegPerG);
    fprintf('Intercept:                        %9.5f deg\n',fit.interceptDeg);
    fprintf('R-squared / RMSE:                %9.5f / %.5f deg\n', ...
        fit.rSquared,fit.rmseDeg);
    fprintf('Used/rejected bins:              %d / %d\n', ...
        fit.usedBinCount,fit.rejectedBinCount);
end

function scale = robustMADScale(values)
    values = values(isfinite(values));
    if isempty(values)
        scale = NaN;
        return
    end
    center = median(values);
    scale = 1.4826*median(abs(values-center));
    if ~isfinite(scale)
        scale = NaN;
    end
end

function critical = studentTCritical95(dof)
    % Uses tinv when available and a normal approximation otherwise.
    if exist('tinv','file') == 2
        critical = tinv(0.975,dof);
    else
        critical = 1.96;
    end
end

function plotTurnBins(bins,colorValue,labelText)
    if isempty(bins)
        return
    end
    if ismember('FitInlier',bins.Properties.VariableNames)
        inlier = bins.FitInlier;
    else
        inlier = true(height(bins),1);
    end

    scatter(bins.MeanLateralG(inlier),bins.MeanExcessSteerDeg(inlier), ...
        55,colorValue,'filled','DisplayName',labelText);
    if any(~inlier)
        scatter(bins.MeanLateralG(~inlier),bins.MeanExcessSteerDeg(~inlier), ...
            70,colorValue,'x','LineWidth',1.8, ...
            'DisplayName',[labelText ' rejected']);
    end
end
