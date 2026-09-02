%% Constant-radius understeer-gradient analysis
% Bundorf constant-radius method:
%
%   road-wheel steer [deg] = K [deg/g] * lateral acceleration [g] + intercept
%
% At constant radius, the slope K is directly the understeer gradient.
% This script analyzes one turn direction at a time so that opposite
% geometric-steer intercepts are never combined in the same fit.

clearvars;
close all;
clc;

%% User settings
DATA_PATH = "T36 Comp Skidpad.csv";

steeringRatio = 1;        % steering-wheel angle / road-wheel angle
wheelbaseM = 1.545;
g = 9.80665;

cutoffHz = 6;               % low-pass cutoff frequency
rawYawRateUnits = "deg/s";  % "rad/s" or "deg/s"

startTimeS = -inf;             % analysis window; use -Inf for entire file
endTimeS = inf;               % use Inf for entire file

turnDirection = 1;           % +1 for positive-yaw turn, -1 for negative-yaw turn
targetRadiusM = 8.2;
radiusToleranceM = 0.15;     % selects 9.5 to 10.0 m with settings above
minLateralG = 0.20;
maxLateralG = 2.00;
maxAbsLongitudinalG = 0.10;
maxAbsYawAccelDegS2 = 10;
minimumFitPoints = 10;

% Optional agreement check between R = V/yaw-rate and R = V^2/ay.
% Increase or set to Inf if accelerometer offsets make this too restrictive.
maxRadiusDisagreementM = 2.0;

writeCsvFiles = true;

%% Import data
opts = detectImportOptions(DATA_PATH);
data = readtable(DATA_PATH,opts);
time = double(data.Time(:));
dt = median(diff(time));
sampleRateHz = 1/dt;
nyquistHz = sampleRateHz/2;
if cutoffHz <= 0 || cutoffHz >= nyquistHz
    error('cutoffHz must be between 0 and the Nyquist frequency (%.3f Hz).',nyquistHz);
end

[bFilter,aFilter] = butter(2,cutoffHz/nyquistHz,'low');

%% Filter vehicle signals
% The CSV steering signal is assumed to be steering-wheel angle in degrees.
steeringWheelDeg = -filterTelemetry(double(data.Steer(:)),bFilter,aFilter);
roadWheelSteerDeg = steeringWheelDeg/steeringRatio;

ax = filterTelemetry(double(data.Ax(:)),bFilter,aFilter); % m/s^2
ay = filterTelemetry(double(data.Ay(:)),bFilter,aFilter); % m/s^2

yawRateRaw = filterTelemetry(double(data.YawR(:)),bFilter,aFilter);
switch rawYawRateUnits
    case "rad/s"
        yawRateRadS = yawRateRaw;
        yawRateDegS = rad2deg(yawRateRaw);
    case "deg/s"
        yawRateDegS = yawRateRaw;
        yawRateRadS = deg2rad(yawRateRaw);
    otherwise
        error('rawYawRateUnits must be "rad/s" or "deg/s".');
end

yawAccelDegS2 = gradient(yawRateDegS)./gradient(time);

%% Wheel-speed-based vehicle speed
% Conversion = tire circumference / (gear ratio * 60), giving m/s from rpm.
frontConversion = 1.24482467/(13.55*60);
rearConversion  = 1.24482467/(13.55*60);

frontRightSpeed = -double(data.FRrpm(:))*frontConversion;
frontLeftSpeed  =  double(data.FLrpm(:))*frontConversion;
rearRightSpeed  = -double(data.RRrpm(:))*rearConversion;
rearLeftSpeed   =  double(data.RLrpm(:))*rearConversion;

frontSpeed = (frontRightSpeed+frontLeftSpeed)/2;
rearSpeed = (rearRightSpeed+rearLeftSpeed)/2;
vehicleSpeedMS = filterTelemetry(abs(frontSpeed),bFilter,aFilter);

%% Derived channels
lateralG = ay/g;
longitudinalG = ax/g;

% Signed radii preserve turn direction. Yaw-rate radius is used for sample
% selection because it is independent of the lateral-acceleration fit axis.
radiusFromYawM = vehicleSpeedMS./yawRateRadS;
radiusFromAyM = vehicleSpeedMS.^2./ay;

% Convert the selected turn to positive coordinates. This allows the same
% fit logic to be used for either left or right turns.
selectedLateralG = turnDirection*lateralG;
selectedSteerDeg = turnDirection*roadWheelSteerDeg;
selectedRadiusYawM = turnDirection*radiusFromYawM;
selectedRadiusAyM = turnDirection*radiusFromAyM;
selectedYawRateDegS = turnDirection*yawRateDegS;

%% Select steady, constant-radius samples
radiusLowerM = targetRadiusM-radiusToleranceM;
radiusUpperM = targetRadiusM+radiusToleranceM;

timeMask = time >= startTimeS & time <= endTimeS;
finiteMask = isfinite(selectedLateralG) & isfinite(selectedSteerDeg) & ...
    isfinite(selectedRadiusYawM) & isfinite(selectedRadiusAyM) & ...
    isfinite(vehicleSpeedMS) & isfinite(yawAccelDegS2);

mask = timeMask & finiteMask & ...
    selectedRadiusYawM >= radiusLowerM & ...
    selectedRadiusYawM <= radiusUpperM & ...
    selectedLateralG >= minLateralG & ...
    selectedLateralG <= maxLateralG & ...
    abs(longitudinalG) <= maxAbsLongitudinalG & ...
    abs(yawAccelDegS2) <= maxAbsYawAccelDegS2 & ...
    abs(selectedRadiusYawM-selectedRadiusAyM) <= maxRadiusDisagreementM;

if nnz(mask) < minimumFitPoints
    error(['Only %d samples passed the filters; at least %d are required. ', ...
        'Check the time window, yaw-rate units, turnDirection, radius band, ', ...
        'and steady-state thresholds.'],nnz(mask),minimumFitPoints);
end

fitLateralG = selectedLateralG(mask);
fitSteerDeg = selectedSteerDeg(mask);

if range(fitLateralG) < 0.05
    warning(['The selected lateral-acceleration range is only %.3f g. ', ...
        'A narrow range makes the fitted gradient highly sensitive to noise.'], ...
        range(fitLateralG));
end

%% Constant-radius fit
% Use a free intercept. At zero lateral acceleration, the extrapolated
% steering is the geometric steering for the selected radius plus offsets.
fitCoefficients = polyfit(fitLateralG,fitSteerDeg,1);
understeerGradientDegPerG = fitCoefficients(1);
fitInterceptDeg = fitCoefficients(2);
fitPredictionDeg = polyval(fitCoefficients,fitLateralG);

residualDeg = fitSteerDeg-fitPredictionDeg;
ssResidual = sum(residualDeg.^2);
ssTotal = sum((fitSteerDeg-mean(fitSteerDeg)).^2);
if ssTotal > eps
    rSquared = 1-ssResidual/ssTotal;
else
    rSquared = NaN;
end
rmseDeg = sqrt(mean(residualDeg.^2));

meanRadiusM = mean(selectedRadiusYawM(mask));
meanSpeedMS = mean(vehicleSpeedMS(mask));
minSelectedG = min(fitLateralG);
maxSelectedG = max(fitLateralG);

% Bicycle-model geometric steer, shown only as an intercept cross-check.
expectedGeometricSteerDeg = rad2deg(atan(wheelbaseM/meanRadiusM));

if understeerGradientDegPerG > 0
    gradientRadPerG = deg2rad(understeerGradientDegPerG);
    characteristicVelocityMS = sqrt(g*wheelbaseM/gradientRadPerG);
    characteristicVelocityKPH = 3.6*characteristicVelocityMS;
    handlingBalance = "understeer";
elseif understeerGradientDegPerG < 0
    characteristicVelocityMS = NaN;
    characteristicVelocityKPH = NaN;
    handlingBalance = "oversteer";
else
    characteristicVelocityMS = Inf;
    characteristicVelocityKPH = Inf;
    handlingBalance = "neutral steer";
end

%% Display results
fprintf('\nConstant-radius understeer analysis\n');
fprintf('-----------------------------------\n');
fprintf('Selected samples:              %d\n',nnz(mask));
fprintf('Turn direction:                %+d\n',turnDirection);
fprintf('Mean radius:                   %.3f m\n',meanRadiusM);
fprintf('Mean speed:                    %.3f m/s\n',meanSpeedMS);
fprintf('Selected lateral-g range:      %.3f to %.3f g\n',minSelectedG,maxSelectedG);
fprintf('Understeer gradient:           %.5f deg/g\n',understeerGradientDegPerG);
fprintf('Fit intercept:                 %.5f deg\n',fitInterceptDeg);
fprintf('Expected geometric intercept:  %.5f deg\n',expectedGeometricSteerDeg);
fprintf('R-squared:                     %.5f\n',rSquared);
fprintf('Fit RMSE:                      %.5f deg\n',rmseDeg);
fprintf('Handling balance:              %s\n',handlingBalance);
if isfinite(characteristicVelocityMS)
    fprintf('Characteristic velocity:       %.3f m/s (%.3f km/h)\n', ...
        characteristicVelocityMS,characteristicVelocityKPH);
else
    fprintf('Characteristic velocity:       not defined for %s\n',handlingBalance);
end

%% Plot fit and diagnostics
figure()
histogram(radiusFromAyM,'BinLimits',[8,12],BinWidth=radiusToleranceM);
figure()
hold on
plot(time,steeringWheelDeg)
plot(time,ay)
plot(time,yawRateRaw)
figure('Name','Constant-Radius Understeer Fit');
scatter(fitLateralG,fitSteerDeg,10,'filled','MarkerFaceAlpha',1);
hold on;
plotX = linspace(min(fitLateralG),max(fitLateralG),200);
plot(plotX,polyval(fitCoefficients,plotX),'LineWidth',2);
grid on;
xlabel('Lateral acceleration (g)');
ylabel('Road-wheel steer angle (deg)');
title(sprintf(['Constant-radius fit: K = %.3f deg/g, R^2 = %.3f, ', ...
    'mean radius = %.2f m'],understeerGradientDegPerG,rSquared,meanRadiusM));
legend('Selected samples','Linear fit','Location','best');

figure('Name','Constant-Radius Diagnostics','Color','w');
tiledlayout(4,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(time,selectedRadiusYawM,'DisplayName','V / yaw rate');
hold on;
plot(time,selectedRadiusAyM,'DisplayName','V^2 / a_y');
yline(radiusLowerM,'--');
yline(radiusUpperM,'--');
scatter(time(mask),selectedRadiusYawM(mask),8,'filled','DisplayName','Selected');
ylabel('Radius (m)');
grid on;
legend('Location','best');

nexttile;
plot(time,selectedLateralG);
hold on;
scatter(time(mask),selectedLateralG(mask),8,'filled');
ylabel('Lateral accel. (g)');
grid on;

nexttile;
plot(time,selectedSteerDeg);
hold on;
scatter(time(mask),selectedSteerDeg(mask),8,'filled');
ylabel('Road-wheel steer (deg)');
grid on;

nexttile;
plot(time,selectedYawRateDegS,'DisplayName','Yaw rate');
hold on;
plot(time,yawAccelDegS2,'DisplayName','Yaw acceleration');
ylabel('deg/s or deg/s^2');
xlabel('Time (s)');
grid on;
legend('Location','best');

%% Optional CSV output
if writeCsvFiles
    selectedData = table( ...
        time(mask),vehicleSpeedMS(mask),selectedLateralG(mask), ...
        selectedSteerDeg(mask),selectedRadiusYawM(mask), ...
        selectedRadiusAyM(mask),selectedYawRateDegS(mask), ...
        yawAccelDegS2(mask),longitudinalG(mask), ...
        'VariableNames',{'timeS','speedMS','lateralG','roadWheelSteerDeg', ...
        'radiusFromYawM','radiusFromAyM','yawRateDegS', ...
        'yawAccelDegS2','longitudinalG'});
    writetable(selectedData,'constant_radius_selected_data.csv');

    summary = table(turnDirection,targetRadiusM,meanRadiusM,nnz(mask), ...
        minSelectedG,maxSelectedG,understeerGradientDegPerG, ...
        fitInterceptDeg,expectedGeometricSteerDeg,rSquared,rmseDeg, ...
        characteristicVelocityMS,characteristicVelocityKPH,handlingBalance, ...
        'VariableNames',{'turnDirection','targetRadiusM','meanRadiusM', ...
        'sampleCount','minLateralG','maxLateralG','KdegPerG', ...
        'fitInterceptDeg','expectedGeometricSteerDeg','rSquared','rmseDeg', ...
        'characteristicVelocityMS','characteristicVelocityKPH', ...
        'handlingBalance'});
    writetable(summary,'constant_radius_understeer_summary.csv');
end

%% Local function
function filtered = filterTelemetry(signal,bFilter,aFilter)
    signal = fillmissing(signal,'linear','EndValues','nearest');
    filtered = filtfilt(bFilter,aFilter,signal);
end
