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
ScaleFactor=0.67;
%% Front Models
FrontFyModelCSV='10 PSI FY.csv';
opts = detectImportOptions(FrontFyModelCSV);
FrontFyModelParams = readtable(FrontFyModelCSV, opts);

FrontMzModelCSV='10 PSI MZ.csv';
FrontMzTable = readtable(FrontMzModelCSV);
FrontMzModelParams = struct;
for name = string(FrontMzTable.Properties.VariableNames)
    x = FrontMzTable.(name);
    FrontMzModelParams.(name) = x(~ismissing(x));
end
%% Rear Models
RearFyModelCSV='10 PSI FY.csv';
opts = detectImportOptions(RearFyModelCSV);
RearFyModelParams = readtable(RearFyModelCSV, opts);

RearMzModelCSV='10 PSI MZ.csv';
RearMzTable = readtable(RearMzModelCSV);
RearMzModelParams = struct;
for name = string(RearMzTable.Properties.VariableNames)
    x = RearMzTable.(name);
    RearMzModelParams.(name) = x(~ismissing(x));
end

%% Sweep inputs
Velocities = 10:2:30;               % [m/s]
SteerAnglesDeg = 0:.1:20;          % handwheel angle [deg]
BodySlipAnglesDeg = -5:.1:5 ;      % vehicle body-slip angle [deg]

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
        FrontFyModelParams,RearFyModelParams,FrontMzModelParams,RearMzModelParams,ScaleFactor);

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
understeer();
% TODO: ymd() is a datetime function and errors with no arguments.
% Replace with your actual plotting call, e.g.:
% figure; contourf(SteerAnglesDeg, BodySlipAnglesDeg, squeeze(ayMap(1,:,:))'); ...
% ymd();

%% ================= Local functions =================

function [ayRow,radiusRow,yawRateRow,yawMomentRow,forceErrorRow, ...
    convergedRow,iterationRow,velResults,peak] = sweepOneVelocity( ...
    Velocity,SteerAnglesDeg,BodySlipAnglesDeg,car,g,mass, ...
    forceToleranceN,maxIterations,relaxation,minimumFz,maxAbsAyG, ...
    FrontFyModelParams,RearFyModelParams,FrontMzModelParams,RearMzModelParams,ScaleFactor)
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
    Beta = deg2rad(BetaDeg); % FIX: added semicolon to suppress console echo

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
                wheelSteer,car.a,car.b,car.FrontTrack,car.RearTrack);

            [candidateFy,candidateMz,candidateWheelState] = findFy( ...
                Fz,SA,Vx,phit,wheelSteer,car.a,car.b, ...
                car.FrontTrack,car.RearTrack, ...
                FrontFyModelParams,RearFyModelParams,FrontMzModelParams,RearMzModelParams,ScaleFactor);

            if ~isfinite(candidateFy) || ~isfinite(candidateMz)
                break
            end

            ayCandidate = candidateFy/mass;
            forceError = abs(candidateFy-mass*ay);

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
    a,b,frontTrack,rearTrack)

    x = [a; a; -b; -b];
    y = [frontTrack/2; -frontTrack/2; rearTrack/2; -rearTrack/2];
    
    VxCG = V*cos(Beta);
    VyCG = V*sin(Beta);
    VxBody = VxCG-r.*y;
    VyBody = VyCG+r.*x;
    
    VxWheel = cos(wheelSteer).*VxBody + sin(wheelSteer).*VyBody;
    VyWheel = -sin(wheelSteer).*VxBody + cos(wheelSteer).*VyBody;
    
    SA = -atan2(VyWheel,VxWheel);
    
    safeVx = sign(VxWheel).*max(abs(VxWheel),1);
    safeVx(safeVx == 0) = 1;
    phit = r./safeVx;
end

function [totalFy,netYawMoment,wheelState] = findFy(Fz,SA,Vx,phit, ...
    wheelSteer,a,b,frontTrack,rearTrack,FrontFyModel,RearFyModel,FrontMzModel,RearMzModel,ScaleFactor)
    nTires = 4;
    FyTire = zeros(nTires,1);
    MzTire = zeros(nTires,1);
    FyTire(1:2) = pacejka4FY(SA(1:2), Fz(1:2), FrontFyModel,ScaleFactor);
    FyTire(3:4) = pacejka4FY(SA(3:4), Fz(3:4), RearFyModel,ScaleFactor);
    MzTire(1:2) = pacejka4Mz(SA(1:2), Fz(1:2), FrontMzModel,ScaleFactor);
    MzTire(3:4) = pacejka4Mz(SA(3:4), Fz(3:4), RearMzModel,ScaleFactor);
    FxBody = -FyTire.*sin(wheelSteer);
    FyBody =  FyTire.*cos(wheelSteer);
    
    x = [a; a; -b; -b];
    y = [frontTrack/2; -frontTrack/2; rearTrack/2; -rearTrack/2];
    MzTire = MzTire+x.*FyBody-y.*FxBody;
    
    totalFy = sum(FyBody);
    netYawMoment = sum(MzTire);
    
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
        wheelState(tireIndex).yawMoment = MzTire(tireIndex);
    end
end
function fy=pacejka4FY(SA, Fz,model,ScaleFactor)
    fy=ScaleFactor*polyval(model.D, Fz) .* sin( polyval(model.C, Fz) .* atan( ...
        polyval(model.B, Fz).*SA - polyval(model.E, Fz).*( ...
        polyval(model.B, Fz).*SA - atan(polyval(model.B, Fz).*SA) ) ) );
end
function Mz = pacejka4Mz(SA,Fz,model,ScaleFactor)
    R0 = model.R0(1);
    Fz0 = model.Fz0(1);
    dfz = (Fz-Fz0)./Fz0;
    D = Fz.*R0.*polyval(model.D,dfz);
    C = polyval(model.C,dfz);
    E = polyval(model.E,dfz);
    K = Fz.*R0.*polyval(model.K,dfz);
    
    denominator = C.*D;
    tiny = 1e-12*max(1,max(abs(denominator)));
    bad = abs(denominator) < tiny;
    denominator(bad) = tiny.*sign(denominator(bad)+eps);
    B = K./denominator;
    
    z = B.*SA;
    q = z-E.*(z-atan(z));
    Mz = -ScaleFactor.*D.*sin(C.*atan(q)); %negative cuz of conventions
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
