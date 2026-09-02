clear;clc;
CAR_PATH= fullfile( 'Car.csv');
opts = detectImportOptions(CAR_PATH);
car=readtable(CAR_PATH,opts);
car.b=car.Wheelbase*car.WeightDist;
car.a=car.Wheelbase*(1-car.WeightDist);
car.FrontMass=car.Weight*car.WeightDist;
car.RearMass=car.Weight*(1-car.WeightDist);
car=table2struct(car);
%% STEADY_STATE_SWEEP
% Four-contact-patch steer/body-slip sweep using ordinary variables.
%
% Required before running:
%   car.Weight, car.Wheelbase, car.FrontTrack, car.RearTrack, car.a, car.b
%
% Every steering/body-slip configuration starts at ay = 0. The script then
% updates yaw rate, load transfer, tire slip angles, and Pacejka forces until
% the lateral-force balance changes by less than forceToleranceN.
%
% Net yaw moment is recorded, not forced to zero. The Mz = 0 contour in the
% resulting map identifies configurations that are also in yaw equilibrium.

assert(exist('car','var') == 1 && isstruct(car), ...
    'Create the car struct in the workspace before running this script.');

%% User variables
velocities = 8:1:10;                    % [m/s]
handwheelAnglesDeg = 20:4:60;          % [deg]
bodySlipAnglesDeg = 0:0.5:20;        % [deg]

steeringRatio = 4.0;                   % handwheel angle / road-wheel angle
steeringMapDeg = [];                   % [handwheel deg, road-wheel deg]
ackermannFraction = 0.30;              % 0 = parallel, 1 = full Ackermann

TIRfile = 'gooch.TIR';
useMode = 111;
pressurePsi = 12.0;
pressurePa = pressurePsi * 6894.757293168;
fyOutputColumn = 2;
mzOutputColumn = 0;                    % Set to mfeval Mz column, or 0 to omit
alphaSign = 1;                        % Reverse if tire SA convention differs
fySign = 1;                            % Reverse if tire Fy convention differs
mzSign = 1;
gamma = 0;                             % inclination angle [rad]
kappa = 0;                             % longitudinal slip
useTurnSlip = true;

weightUnits = 'N';                     % 'N', 'kg', 'lbf', or 'lbm'
g = 9.80665;
frontDownforceTable = [];              % [speed m/s, positive downforce N]
rearDownforceTable = [];               % [speed m/s, positive downforce N]
minimumFz = 5;                         % [N]
minimumTireSpeed = 0.5;                % [m/s]

forceToleranceN = 1.0;                 % lateral-force balance tolerance [N]
maxIterations = 10;
relaxation = 0.5;                      % 0 < relaxation <= 1
maxAbsAyG = 5.0;
stopAfterFailedSteer = true;           % stop when no beta converges at a steer
verbose = true;
outputFile = '';                       % e.g. 'steady_state_results.mat'

%% Vehicle constants
mass = convertWeightToMass(car.Weight, weightUnits, g);
staticWeight = mass * g;
wheelbase = car.Wheelbase;
a = car.a;
b = car.b;
frontTrack = car.FrontTrack;
rearTrack = car.RearTrack;

assert(mass > 0 && wheelbase > 0 && frontTrack > 0 && rearTrack > 0, ...
    'Vehicle mass and dimensions must be positive.');
assert(a > 0 && b > 0, 'car.a and car.b must be positive.');
assert(relaxation > 0 && relaxation <= 1, ...
    'relaxation must be greater than zero and no greater than one.');

if abs((a + b) - wheelbase) <= max(1e-6, 1e-3 * wheelbase)
    frontStaticFraction = b / wheelbase;
elseif isfield(car,'WeightDist')
    frontStaticFraction = car.WeightDist;
else
    frontStaticFraction = 0.5;
end

if isfield(car,'cGz')
    cgHeight = car.cGz;
elseif isfield(car,'CGHeight')
    cgHeight = car.CGHeight;
else
    cgHeight = 0;
end

if isfield(car,'TLLTD')
    TLLTD = car.TLLTD;
else
    TLLTD = 0.5;
end

assert(frontStaticFraction >= 0 && frontStaticFraction <= 1, ...
    'Front static weight fraction must be between zero and one.');
assert(TLLTD >= 0 && TLLTD <= 1, ...
    'TLLTD must be between zero and one.');

%% Output arrays
nVelocity = numel(velocities);
nSteer = numel(handwheelAnglesDeg);
nBeta = numel(bodySlipAnglesDeg);

ayGrid = NaN(nVelocity,nSteer,nBeta);
yawMomentGrid = NaN(nVelocity,nSteer,nBeta);
yawRateGrid = NaN(nVelocity,nSteer,nBeta);
radiusGrid = NaN(nVelocity,nSteer,nBeta);
forceErrorGrid = NaN(nVelocity,nSteer,nBeta);
iterationGrid = zeros(nVelocity,nSteer,nBeta);
convergedGrid = false(nVelocity,nSteer,nBeta);

resultTemplate = struct( ...
    'velocity',[], ...
    'handwheelAngleDeg',[], ...
    'centerRoadWheelAngleDeg',[], ...
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
results = repmat(resultTemplate,0,1);

peakTemplate = struct( ...
    'velocity',[], ...
    'maxAy',NaN, ...
    'maxAyG',NaN, ...
    'maxAyConfiguration',[], ...
    'maxAbsAy',NaN, ...
    'maxAbsAyG',NaN, ...
    'maxAbsAyConfiguration',[]);
velocityPeaks = repmat(peakTemplate,nVelocity,1);

%% Velocity, steering, and body-slip sweep
for velocityIndex = 1:nVelocity
    V = velocities(velocityIndex);
    assert(V > 0, 'All velocities must be positive.');

    frontDownforce = tableLookup(frontDownforceTable,V);
    rearDownforce = tableLookup(rearDownforceTable,V);
    frontAxleLoad = staticWeight * frontStaticFraction + frontDownforce;
    rearAxleLoad = staticWeight * (1 - frontStaticFraction) + rearDownforce;

    firstResultAtVelocity = numel(results) + 1;

    if verbose
        fprintf('Sweeping V = %.1f m/s\n',V);
    end

    for steeringIndex = 1:nSteer
        handwheelDeg = handwheelAnglesDeg(steeringIndex);
        centerSteer = roadWheelAngle(handwheelDeg,steeringRatio,steeringMapDeg);
        [deltaFL,deltaFR] = ackermannAngles(centerSteer,wheelbase, ...
            frontTrack,ackermannFraction);
        wheelSteer = [deltaFL;deltaFR;0;0];
        convergedAtThisSteer = false;

        for betaIndex = 1:nBeta
            beta = deg2rad(bodySlipAnglesDeg(betaIndex));
            ay = 0;
            converged = false;
            forceError = Inf;
            totalFy = NaN;
            netYawMoment = NaN;
            wheelState = struct([]);

            for iteration = 1:maxIterations
                yawRate = ay / V;

                [Fz,loadsValid] = calculateWheelLoads(frontAxleLoad, ...
                    rearAxleLoad,mass,ay,cgHeight,frontTrack,rearTrack, ...
                    TLLTD,minimumFz);
                if ~loadsValid
                    break
                end

                [candidateFy,candidateMz,candidateWheelState] = ...
                    calculateVehicleForces(V,beta,yawRate,wheelSteer,Fz, ...
                    a,b,frontTrack,rearTrack,TIRfile,useMode,pressurePa, ...
                    fyOutputColumn,mzOutputColumn,alphaSign,fySign,mzSign, ...
                    gamma,kappa,useTurnSlip,minimumTireSpeed);

                if ~isfinite(candidateFy) || ~isfinite(candidateMz)
                    break
                end

                % This is the difference between the tire force currently
                % produced and the force required by the current ay guess.
                forceError = abs(candidateFy - mass * ay);
                ayCandidate = candidateFy / mass;

                if abs(ayCandidate) > maxAbsAyG * g
                    break
                end

                totalFy = candidateFy;
                netYawMoment = candidateMz;
                wheelState = candidateWheelState;

                if forceError <= forceToleranceN
                    ay = ayCandidate;
                    converged = true;
                    break
                end

                ay = ay + relaxation * (ayCandidate - ay);
            end

            iterationGrid(velocityIndex,steeringIndex,betaIndex) = iteration;
            forceErrorGrid(velocityIndex,steeringIndex,betaIndex) = forceError;

            if ~converged
                continue
            end

            % Re-evaluate once at the accepted ay so every saved quantity
            % corresponds to the same final vehicle state.
            yawRate = ay / V;
            [Fz,loadsValid] = calculateWheelLoads(frontAxleLoad, ...
                rearAxleLoad,mass,ay,cgHeight,frontTrack,rearTrack, ...
                TLLTD,minimumFz);
            if ~loadsValid
                continue
            end

            [totalFy,netYawMoment,wheelState] = calculateVehicleForces( ...
                V,beta,yawRate,wheelSteer,Fz,a,b,frontTrack,rearTrack, ...
                TIRfile,useMode,pressurePa,fyOutputColumn,mzOutputColumn, ...
                alphaSign,fySign,mzSign,gamma,kappa,useTurnSlip, ...
                minimumTireSpeed);

            forceError = abs(totalFy - mass * ay);
            if ~isfinite(forceError) || forceError > forceToleranceN
                continue
            end

            if abs(yawRate) > 1e-10
                radius = V / yawRate;
            else
                radius = Inf;
            end

            convergedGrid(velocityIndex,steeringIndex,betaIndex) = true;
            ayGrid(velocityIndex,steeringIndex,betaIndex) = ay;
            yawMomentGrid(velocityIndex,steeringIndex,betaIndex) = netYawMoment;
            yawRateGrid(velocityIndex,steeringIndex,betaIndex) = yawRate;
            radiusGrid(velocityIndex,steeringIndex,betaIndex) = radius;
            forceErrorGrid(velocityIndex,steeringIndex,betaIndex) = forceError;
            convergedAtThisSteer = true;

            vehicleState = struct( ...
                'bodySlipRad',beta, ...
                'bodySlipDeg',rad2deg(beta), ...
                'yawRate',yawRate, ...
                'lateralAcceleration',ay, ...
                'lateralAccelerationG',ay/g, ...
                'radius',radius, ...
                'frontAxleLoad',frontAxleLoad, ...
                'rearAxleLoad',rearAxleLoad, ...
                'frontDownforce',frontDownforce, ...
                'rearDownforce',rearDownforce);

            newResult = resultTemplate;
            newResult.velocity = V;
            newResult.handwheelAngleDeg = handwheelDeg;
            newResult.centerRoadWheelAngleDeg = rad2deg(centerSteer);
            newResult.bodySlipDeg = rad2deg(beta);
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
            results(end+1,1) = newResult; %#ok<SAGROW>
        end

        if stopAfterFailedSteer && ~convergedAtThisSteer
            if verbose
                fprintf('  No converged beta at %.1f deg handwheel; stopping steer sweep.\n', ...
                    handwheelDeg);
            end
            break
        end
    end

    velocityPeaks(velocityIndex).velocity = V;
    lastResultAtVelocity = numel(results);

    if lastResultAtVelocity >= firstResultAtVelocity
        velocityResults = results(firstResultAtVelocity:lastResultAtVelocity);
        ayValues = [velocityResults.lateralAcceleration];

        [maxAy,maxIndex] = max(ayValues);
        [~,maxAbsIndex] = max(abs(ayValues));
        maxAbsAy = ayValues(maxAbsIndex);

        velocityPeaks(velocityIndex).maxAy = maxAy;
        velocityPeaks(velocityIndex).maxAyG = maxAy/g;
        velocityPeaks(velocityIndex).maxAyConfiguration = ...
            velocityResults(maxIndex);
        velocityPeaks(velocityIndex).maxAbsAy = maxAbsAy;
        velocityPeaks(velocityIndex).maxAbsAyG = maxAbsAy/g;
        velocityPeaks(velocityIndex).maxAbsAyConfiguration = ...
            velocityResults(maxAbsIndex);
    end
end

%% Gridded outputs
sweep = struct;
sweep.velocities = velocities;
sweep.handwheelAnglesDeg = handwheelAnglesDeg;
sweep.bodySlipAnglesDeg = bodySlipAnglesDeg;
sweep.ay = ayGrid;
sweep.yawMoment = yawMomentGrid;
sweep.yawRate = yawRateGrid;
sweep.radius = radiusGrid;
sweep.forceBalanceError = forceErrorGrid;
sweep.iterations = iterationGrid;
sweep.converged = convergedGrid;

if ~isempty(outputFile)
    save(outputFile,'results','velocityPeaks','sweep','car');
end

fprintf('Saved %d converged configurations in results.\n',numel(results));

%% Local functions
function mass = convertWeightToMass(value,units,g)
    switch lower(units)
        case 'n'
            mass = value / g;
        case 'kg'
            mass = value;
        case 'lbf'
            mass = value * 4.4482216152605 / g;
        case 'lbm'
            mass = value * 0.45359237;
        otherwise
            error('weightUnits must be N, kg, lbf, or lbm.');
    end
end

function force = tableLookup(table,V)
    if isempty(table)
        force = 0;
        return
    end
    assert(size(table,2) == 2, ...
        'Downforce tables must contain [speed, force] columns.');
    force = interp1(table(:,1),table(:,2),V,'linear','extrap');
end

function delta = roadWheelAngle(handwheelDeg,steeringRatio,steeringMapDeg)
    if isempty(steeringMapDeg)
        delta = deg2rad(handwheelDeg / steeringRatio);
    else
        assert(size(steeringMapDeg,2) == 2, ...
            'steeringMapDeg must contain [handwheel deg, road-wheel deg].');
        roadWheelDeg = interp1(steeringMapDeg(:,1),steeringMapDeg(:,2), ...
            handwheelDeg,'linear','extrap');
        delta = deg2rad(roadWheelDeg);
    end
end

function [deltaLeft,deltaRight] = ackermannAngles(centerDelta,L,track,fraction)
    if abs(centerDelta) < 1e-12
        deltaLeft = 0;
        deltaRight = 0;
        return
    end

    turnSign = sign(centerDelta);
    referenceRadius = L / tan(abs(centerDelta));
    inner = atan(L / max(referenceRadius - track/2,eps));
    outer = atan(L / (referenceRadius + track/2));

    inner = abs(centerDelta) + fraction * (inner - abs(centerDelta));
    outer = abs(centerDelta) + fraction * (outer - abs(centerDelta));

    if turnSign > 0
        deltaLeft = inner;
        deltaRight = outer;
    else
        deltaLeft = -outer;
        deltaRight = -inner;
    end
end

function [Fz,valid] = calculateWheelLoads(frontLoad,rearLoad,mass,ay, ...
        cgHeight,frontTrack,rearTrack,TLLTD,minimumFz)

    frontTransfer = TLLTD * mass * ay * cgHeight / frontTrack;
    rearTransfer = (1 - TLLTD) * mass * ay * cgHeight / rearTrack;

    % Positive ay is a left turn, so load transfers to the right wheels.
    Fz = [frontLoad/2 - frontTransfer; ...
          frontLoad/2 + frontTransfer; ...
          rearLoad/2 - rearTransfer; ...
          rearLoad/2 + rearTransfer];

    valid = all(isfinite(Fz)) && all(Fz >= minimumFz);
end

function [totalFy,netYawMoment,wheelState] = calculateVehicleForces( ...
        V,beta,yawRate,wheelSteer,Fz,a,b,frontTrack,rearTrack, ...
        TIRfile,useMode,pressurePa,fyColumn,mzColumn,alphaSign,fySign, ...
        mzSign,gamma,kappa,useTurnSlip,minimumTireSpeed)

    wheelNames = {'FL','FR','RL','RR'};
    x = [a;a;-b;-b];
    y = [frontTrack/2;-frontTrack/2;rearTrack/2;-rearTrack/2];

    cgVx = V * cos(beta);
    cgVy = V * sin(beta);
    totalFy = 0;
    netYawMoment = 0;

    stateTemplate = struct('name','','x',0,'y',0,'steerAngle',0, ...
        'Fz',0,'VxBody',0,'VyBody',0,'VxTire',0,'VyTire',0, ...
        'slipAngle',0,'turnSlip',0,'FyTire',0,'MzTire',0, ...
        'FxBody',0,'FyBody',0,'yawMoment',0);
    wheelState = repmat(stateTemplate,4,1);

    for wheel = 1:4
        delta = wheelSteer(wheel);
        vxBody = cgVx - yawRate * y(wheel);
        vyBody = cgVy + yawRate * x(wheel);

        vxTire = cos(delta) * vxBody + sin(delta) * vyBody;
        vyTire = -sin(delta) * vxBody + cos(delta) * vyBody;
        pacejkaVx = max(abs(vxTire),minimumTireSpeed);
        alpha = alphaSign * atan2(vyTire,pacejkaVx);

        if useTurnSlip
            phit = yawRate / pacejkaVx;
        else
            phit = 0;
        end

        inputs = [Fz(wheel),kappa,alpha,gamma,phit,pacejkaVx,pressurePa];
        output = mfeval(TIRfile,inputs,useMode);

        assert(size(output,2) >= fyColumn, ...
            'fyOutputColumn exceeds the number of mfeval output columns.');
        FyTire = fySign * output(1,fyColumn);

        if mzColumn > 0
            assert(size(output,2) >= mzColumn, ...
                'mzOutputColumn exceeds the number of mfeval output columns.');
            MzTire = mzSign * output(1,mzColumn);
        else
            MzTire = 0;
        end

        % Pure lateral tire force rotated into vehicle coordinates.
        FxBody = -sin(delta) * FyTire;
        FyBody = cos(delta) * FyTire;
        wheelMoment = x(wheel) * FyBody - y(wheel) * FxBody + MzTire;

        totalFy = totalFy + FyBody;
        netYawMoment = netYawMoment + wheelMoment;

        wheelState(wheel).name = wheelNames{wheel};
        wheelState(wheel).x = x(wheel);
        wheelState(wheel).y = y(wheel);
        wheelState(wheel).steerAngle = delta;
        wheelState(wheel).Fz = Fz(wheel);
        wheelState(wheel).VxBody = vxBody;
        wheelState(wheel).VyBody = vyBody;
        wheelState(wheel).VxTire = vxTire;
        wheelState(wheel).VyTire = vyTire;
        wheelState(wheel).slipAngle = alpha;
        wheelState(wheel).turnSlip = phit;
        wheelState(wheel).FyTire = FyTire;
        wheelState(wheel).MzTire = MzTire;
        wheelState(wheel).FxBody = FxBody;
        wheelState(wheel).FyBody = FyBody;
        wheelState(wheel).yawMoment = wheelMoment;
    end
end

