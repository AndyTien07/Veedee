%% Load your car struct
CAR_PATH = fullfile('Car.csv');
warning('off', 'all');
opts = detectImportOptions(CAR_PATH);
car = readtable(CAR_PATH,opts);

car.b = car.Wheelbase .* car.WeightDist;
car.a = car.Wheelbase .* (1-car.WeightDist);
car.FrontMass = car.Weight .* car.WeightDist;
car.RearMass = car.Weight .* (1-car.WeightDist);
car = table2struct(car);
car.FrontToeRad = deg2rad(car.FrontToe);
car.RearToeRad = deg2rad(car.RearToe);

%% Pacejka inputs
    TIRfile = 'The One.TIR';
warning('off', 'all'); %sometimes one has to violate physics to model it
useMode = 111;
pressurePsi = 14.0;
pressurePa = pressurePsi * 6894.757293168;
fyOutputColumn = 2;
alphaSign = 1;       % Change to +1 only if required by your TIR convention
fySign = 1;           % Change to -1 only if required by your TIR convention
gamma = 0;
kappa = 0;
minimumTireSpeed = 0.75;

%% Sweep inputs
Velocities = 10:2:30;               % [m/s]
SteerAnglesDeg = 0:1:12;          % handwheel angle [deg]
BodySlipAnglesDeg = -2:1:3 ;      % vehicle body-slip angle [deg]

%% Iteration inputs
g = 9.80665;
mass = car.Weight/g;
forceToleranceN = 1;
maxIterations = 500;
relaxation = .5;
minimumFz = 5;
maxAbsAyG = 5;


%% Output storage
nVelocity = numel(Velocities);
nSteer = numel(SteerAnglesDeg);
nBeta = numel(BodySlipAnglesDeg);

ayMap = NaN(nVelocity,nSteer,nBeta);
radiusMap = NaN(nVelocity,nSteer,nBeta);
yawRateMap = NaN(nVelocity,nSteer,nBeta);
yawMomentMap = NaN(nVelocity,nSteer,nBeta);
forceErrorMap = NaN(nVelocity,nSteer,nBeta);
convergedMap = false(nVelocity,nSteer,nBeta);
iterationMap = zeros(nVelocity,nSteer,nBeta);

peakTemplate = struct( ...
    'velocity',[], ...
    'maxAy',NaN, ...
    'maxAyG',NaN, ...
    'configuration',[]);
velocityPeaks = repmat(peakTemplate,nVelocity,1);

resultsCell = cell(nVelocity,1);

%% Start a parallel pool 
p = gcp('nocreate');
if isempty(p)
    parpool; 
end

% ---- Parallel path (no live plotting) ----
parfor velocityIndex = 1:nVelocity
    Velocity = Velocities(velocityIndex);
    fprintf('Sweeping %.1f m/s\n',Velocity);

    [ayRow,radiusRow,yawRateRow,yawMomentRow,forceErrorRow, ...
        convergedRow,iterationRow,velResults,peak] = sweepOneVelocity( ...
        Velocity,SteerAnglesDeg,BodySlipAnglesDeg,car,g,mass, ...
        forceToleranceN,maxIterations,relaxation,minimumFz,maxAbsAyG, ...
        TIRfile,useMode,pressurePa,fyOutputColumn,alphaSign,fySign, ...
        gamma,kappa,minimumTireSpeed,false,[]);

    ayMap(velocityIndex,:,:) = ayRow;
    radiusMap(velocityIndex,:,:) = radiusRow;
    yawRateMap(velocityIndex,:,:) = yawRateRow;
    yawMomentMap(velocityIndex,:,:) = yawMomentRow;
    forceErrorMap(velocityIndex,:,:) = forceErrorRow;
    convergedMap(velocityIndex,:,:) = convergedRow;
    iterationMap(velocityIndex,:,:) = iterationRow;
    resultsCell{velocityIndex} = velResults;

    peak.velocity = Velocity;
    velocityPeaks(velocityIndex) = peak;
end

results = vertcat(resultsCell{:});
fprintf('Finished: %d configurations converged.\n',numel(results));

%% Visualization
%TITLE BLOCK this is very important if u dont wanna manually note down the setups
fields = fieldnames(car);
values = struct2cell(car);

fig = figure('Name','Car Parameters','Color','w');
uitable(fig,'Data',[fields values], ...
    'ColumnName',{'Parameter','Value'}, ...
    'Units','normalized','Position',[0 0 1 1]);
fontsize(18,"points")
%Graphs u wanna show
ymd();

%% ================= Local functions =================

function [ayRow,radiusRow,yawRateRow,yawMomentRow,forceErrorRow, ...
    convergedRow,iterationRow,velResults,peak] = sweepOneVelocity( ...
    Velocity,SteerAnglesDeg,BodySlipAnglesDeg,car,g,mass, ...
    forceToleranceN,maxIterations,relaxation,minimumFz,maxAbsAyG, ...
    TIRfile,useMode,pressurePa,fyOutputColumn,alphaSign,fySign, ...
    gamma,kappa,minimumTireSpeed,diagnosticOn,diagHandles)
% Runs the full steer x body-slip sweep for a single velocity.
% Returns sliced result matrices (nSteer x nBeta) plus a results struct
% array for all converged configurations and the peak-Ay configuration.

nSteer = numel(SteerAnglesDeg);
nBeta = numel(BodySlipAnglesDeg);

ayRow = NaN(nSteer,nBeta);
radiusRow = NaN(nSteer,nBeta);
yawRateRow = NaN(nSteer,nBeta);
yawMomentRow = NaN(nSteer,nBeta);
forceErrorRow = NaN(nSteer,nBeta);
convergedRow = false(nSteer,nBeta);
iterationRow = zeros(nSteer,nBeta);

resultTemplate = struct( ...
    'velocity',[], ...
    'handwheelAngleDeg',[], ...
    'bodySlipDeg',[], ...
    'lateralAcceleration',[], ...
    'lateralAccelerationG',[], ...
    'yawRate',[], ...
    'radius',[], ...
    'netYawMoment',[], ...
    'totalLateralForce',[], ...
    'forceBalanceError',[], ...
    'iterations',[], ...
    'vehicleState',[], ...
    'wheelState',[]);
velResults = repmat(resultTemplate,0,1);

frontDownforce = Velocity*Velocity*car.FrontalArea*car.Cl*1.204*car.COPX;
rearDownforce = Velocity*Velocity*car.FrontalArea*car.Cl*1.204*(1-car.COPX);
frontAxleLoad = car.FrontMass + frontDownforce;
rearAxleLoad = car.RearMass + rearDownforce;

for betaIndex = 1:nBeta
    BetaDeg = BodySlipAnglesDeg(betaIndex);
    Beta = deg2rad(BetaDeg)

    for steerIndex = 1:nSteer
        SteerDeg = SteerAnglesDeg(steerIndex);
        centerSteer = deg2rad(SteerDeg);
        [deltaFL,deltaFR] = ackermannSteer(centerSteer, ...
            car.Wheelbase,car.FrontTrack,car.Ackerman);

        wheelSteer = [deltaFL-car.FrontToeRad; deltaFR+car.FrontToeRad; ...
            -car.RearToeRad; car.RearToeRad];

        ay = car.MU*g;
        converged = false;
        forceError = Inf;
        totalFy = NaN;
        netYawMoment = NaN;
        wheelState = struct([]);

       

        for iteration = 1:maxIterations
            yawRate = ay/Velocity;

            [Fz,loadsValid] = findFz(frontAxleLoad,rearAxleLoad, ...
                mass,ay,car.CGZ,car.FrontTrack,car.RearTrack, ...
                car.TLLTD,minimumFz);
            
            if ~loadsValid
                break
            end

            [SA,Vx,phit] = findSlipAngles(Velocity,Beta,yawRate, ...
                wheelSteer,car.a,car.b,car.FrontTrack,car.RearTrack, ...
                alphaSign,minimumTireSpeed);

            [candidateFy,candidateMz,candidateWheelState] = findFy( ...
                Fz,SA,Vx,phit,wheelSteer,car.a,car.b, ...
                car.FrontTrack,car.RearTrack,TIRfile,useMode, ...
                pressurePa,fyOutputColumn,fySign,gamma,kappa);

            if ~isfinite(candidateFy) || ~isfinite(candidateMz)
                break
            end

            ayCandidate = candidateFy/mass;
            forceError = abs(candidateFy-mass*ay);

            if diagnosticOn && isgraphics(diagHandles.figure) && ...
                    (iteration == 1 || mod(iteration,diagHandles.updateInterval) == 0 || ...
                    forceError <= forceToleranceN)

                addpoints(diagHandles.convergenceLine,iteration,max(forceError,eps));
                addpoints(diagHandles.AyLine,iteration,candidateFy);
                diagHandles.loadBars.YData = Fz(:).';
                diagHandles.SaBars.YData = SA(:).';

                subtitleText = sprintf( ...
                    'a_y = %+.3f g | sum F_y = %+.1f N | residual = %.2f N', ...
                    ay/g,totalFy,forceError);
                diagHandles.convergenceAxes.Title.String = {sprintf( ...
                    'V = %.1f m/s, steer = %.1f deg, beta = %.1f deg', ...
                    Velocity,SteerDeg,BetaDeg),subtitleText};

                drawnow limitrate;
            end

            if abs(ayCandidate) > maxAbsAyG*g
                break
            end

            totalFy = candidateFy;
            netYawMoment = candidateMz;
            wheelState = candidateWheelState;

            if forceError <= forceToleranceN
                converged = true;
                break
            end

            ay = ay + relaxation*(ayCandidate-ay);
        end

        iterationRow(steerIndex,betaIndex) = iteration;
        forceErrorRow(steerIndex,betaIndex) = forceError;

        if ~converged
            continue
        end

        yawRate = ay/Velocity;
        if abs(yawRate) < 1e-12
            radius = Inf;
        else
            radius = Velocity/yawRate;
        end

        ayRow(steerIndex,betaIndex) = ay;
        radiusRow(steerIndex,betaIndex) = radius;
        yawRateRow(steerIndex,betaIndex) = yawRate;
        yawMomentRow(steerIndex,betaIndex) = netYawMoment;
        convergedRow(steerIndex,betaIndex) = true;

        vehicleState = struct( ...
            'bodySlipRad',Beta, ...
            'bodySlipDeg',BetaDeg, ...
            'yawRate',yawRate, ...
            'lateralAcceleration',ay, ...
            'lateralAccelerationG',ay/g, ...
            'radius',radius, ...
            'frontAxleLoad',frontAxleLoad, ...
            'rearAxleLoad',rearAxleLoad, ...
            'wheelSteerRad',wheelSteer, ...
            'wheelSteerDeg',rad2deg(wheelSteer));

        newResult = resultTemplate;
        newResult.velocity = Velocity;
        newResult.handwheelAngleDeg = SteerDeg;
        newResult.bodySlipDeg = BetaDeg;
        newResult.lateralAcceleration = ay;
        newResult.lateralAccelerationG = ay/g;
        newResult.yawRate = yawRate;
        newResult.radius = radius;
        newResult.netYawMoment = netYawMoment;
        newResult.totalLateralForce = totalFy;
        newResult.forceBalanceError = forceError;
        newResult.iterations = iteration;
        newResult.vehicleState = vehicleState;
        newResult.wheelState = wheelState;
        velResults(end+1,1) = newResult; %#ok<AGROW>
    end
end

peak = struct('velocity',[],'maxAy',NaN,'maxAyG',NaN,'configuration',[]);
if ~isempty(velResults)
    [maxAy,maxIndex] = max([velResults.lateralAcceleration]);
    peak.maxAy = maxAy;
    peak.maxAyG = maxAy/g;
    peak.configuration = velResults(maxIndex);
end
end

function [Fz,valid] = findFz(frontLoad,rearLoad,mass,ay,cGz, ...
    frontTrack,rearTrack,TLLTD,minimumFz)

frontTransfer = TLLTD*mass*ay*cGz/frontTrack;
rearTransfer = (1-TLLTD)*mass*ay*cGz/rearTrack;

Fz = [frontLoad/2-frontTransfer; ...
    frontLoad/2+frontTransfer; ...
    rearLoad/2-rearTransfer; ...
    rearLoad/2+rearTransfer];

valid = all(isfinite(Fz)) && all(Fz >= minimumFz);
end

function [SA,VxWheel,phit] = findSlipAngles(V,Beta,r,wheelSteer, ...
    a,b,frontTrack,rearTrack,alphaSign,minimumTireSpeed)

x = [a; a; -b; -b];
y = [frontTrack/2; -frontTrack/2; rearTrack/2; -rearTrack/2];

VxCG = V*cos(Beta);
VyCG = V*sin(Beta);
VxBody = VxCG-r.*y;
VyBody = VyCG+r.*x;

VxWheel = cos(wheelSteer).*VxBody + sin(wheelSteer).*VyBody;
VyWheel = -sin(wheelSteer).*VxBody + cos(wheelSteer).*VyBody;

SA = alphaSign*atan2(VyWheel,VxWheel);

safeVx = sign(VxWheel).*max(abs(VxWheel),minimumTireSpeed);
safeVx(safeVx == 0) = minimumTireSpeed;
phit = r./safeVx;
end

function [totalFy,netYawMoment,wheelState] = findFy(Fz,SA,Vx,phit, ...
    wheelSteer,a,b,frontTrack,rearTrack,TIRfile,useMode,pressurePa, ...
    fyOutputColumn,fySign,gamma,kappa)

nTires = 4;
inputs = [Fz, ...
    repmat(kappa,nTires,1), ...
    SA, ...
    repmat(gamma,nTires,1), ...
    phit, ...
    Vx, ...
    repmat(pressurePa,nTires,1)];

output = mfeval(TIRfile,inputs,useMode);
FyTire = fySign*output(:,fyOutputColumn);
MzTire = output(:,6);

FxBody = -FyTire.*sin(wheelSteer);
FyBody =  FyTire.*cos(wheelSteer);

x = [a; a; -b; -b];
y = [frontTrack/2; -frontTrack/2; rearTrack/2; -rearTrack/2];
wheelYawMoment = x.*FyBody-y.*FxBody;

totalFy = sum(FyBody);
netYawMoment = sum(wheelYawMoment)+sum(MzTire);

wheelNames = {'FL','FR','RL','RR'};
wheelState = repmat(struct( ...
    'name','', ...
    'Fz',0, ...
    'slipAngleRad',0, ...
    'slipAngleDeg',0, ...
    'Vx',0, ...
    'turnSlip',0, ...
    'steerRad',0, ...
    'steerDeg',0, ...
    'FyTire',0, ...
    'FxBody',0, ...
    'FyBody',0, ...
    'yawMoment',0),nTires,1);

for tireIndex = 1:nTires
    wheelState(tireIndex).name = wheelNames{tireIndex};
    wheelState(tireIndex).Fz = Fz(tireIndex);
    wheelState(tireIndex).slipAngleRad = SA(tireIndex);
    wheelState(tireIndex).slipAngleDeg = rad2deg(SA(tireIndex));
    wheelState(tireIndex).Vx = Vx(tireIndex);
    wheelState(tireIndex).turnSlip = phit(tireIndex);
    wheelState(tireIndex).steerRad = wheelSteer(tireIndex);
    wheelState(tireIndex).steerDeg = rad2deg(wheelSteer(tireIndex));
    wheelState(tireIndex).FyTire = FyTire(tireIndex);
    wheelState(tireIndex).FxBody = FxBody(tireIndex);
    wheelState(tireIndex).FyBody = FyBody(tireIndex);
    wheelState(tireIndex).yawMoment = wheelYawMoment(tireIndex);
end
end

function [deltaFL,deltaFR] = ackermannSteer(centerSteer,wheelbase, ...
    frontTrack,ackermannFraction)

if abs(centerSteer) < 1e-12
    deltaFL = 0;
    deltaFR = 0;
    return
end

steerSign = sign(centerSteer);
centerMagnitude = abs(centerSteer);
geometricRadius = wheelbase/tan(centerMagnitude);

innerSteer = atan(wheelbase/max(geometricRadius-frontTrack/2,eps));
outerSteer = atan(wheelbase/(geometricRadius+frontTrack/2));

innerSteer = centerMagnitude + ackermannFraction*(innerSteer-centerMagnitude);
outerSteer = centerMagnitude + ackermannFraction*(outerSteer-centerMagnitude);

if steerSign > 0
    deltaFL = innerSteer;
    deltaFR = outerSteer;
else
    deltaFL = -outerSteer;
    deltaFR = -innerSteer;
end
end
