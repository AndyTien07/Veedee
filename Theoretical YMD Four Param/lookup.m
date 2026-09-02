%% Look up the closest converged car state to a requested (steer, V, ay[g])
% Run steady_state_sweep.m (or vehicle_sweep_parallel.m) first so that
% the "results" struct array exists in the workspace. "results" is the
% flat list of every converged configuration, each carrying full
% vehicleState and wheelState detail (Fz, slip angles, forces, etc.).

g = 9.80665;

assert(exist('results','var') == 1 && ~isempty(results), ...
    'Run steady_state_sweep.m (or vehicle_sweep_parallel.m) first so "results" exists.');

%% ---- USER INPUT: edit these three values ----
targetSteerDeg = 16.5;      % handwheel steer angle [deg]
targetVelocity = 10;     % vehicle speed [m/s]
targetAyG      = 1.63;      % lateral acceleration [g]
%% ----------------------------------------------

steerAll = [results.handwheelAngleDeg].';
velocityAll = [results.velocity].';
ayGAll = [results.lateralAccelerationG].';

% Normalize each axis by its own range so no single unit dominates the
% distance metric (degrees vs. m/s vs. g are not directly comparable).
steerRange = max(steerAll)-min(steerAll);
velocityRange = max(velocityAll)-min(velocityAll);
ayGRange = max(ayGAll)-min(ayGAll);
steerRange(steerRange == 0) = 1;
velocityRange(velocityRange == 0) = 1;
ayGRange(ayGRange == 0) = 1;

normSteerError = (steerAll-targetSteerDeg)/steerRange;
normVelocityError = (velocityAll-targetVelocity)/velocityRange;
normAyGError = (ayGAll-targetAyG)/ayGRange;

combinedDistance = sqrt(normSteerError.^2 + normVelocityError.^2 + normAyGError.^2);
[~,closestIndex] = min(combinedDistance);
closest = results(closestIndex);

fprintf('\n================ REQUESTED ================\n');
fprintf('Steer angle   : %.2f deg\n',targetSteerDeg);
fprintf('Velocity      : %.2f m/s\n',targetVelocity);
fprintf('Lateral accel : %.3f g\n',targetAyG);

fprintf('\n================ CLOSEST MATCH ================\n');
fprintf('Steer angle   : %.2f deg  (delta %.2f deg)\n', ...
    closest.handwheelAngleDeg,closest.handwheelAngleDeg-targetSteerDeg);
fprintf('Velocity      : %.2f m/s  (delta %.2f m/s)\n', ...
    closest.velocity,closest.velocity-targetVelocity);
fprintf('Lateral accel : %.3f g   (delta %.3f g)\n', ...
    closest.lateralAccelerationG,closest.lateralAccelerationG-targetAyG);
fprintf('Body-slip     : %.2f deg\n',closest.bodySlipDeg);

fprintf('\n---- Vehicle-level state ----\n');
fprintf('Yaw rate            : %.4f rad/s\n',closest.yawRate);
fprintf('Turn radius         : %.2f m\n',closest.radius);
fprintf('Net yaw moment      : %.2f N m\n',closest.netYawMoment);
fprintf('Total lateral force : %.2f N\n',closest.totalLateralForce);
fprintf('Force balance error : %.4f N\n',closest.forceBalanceError);
fprintf('Iterations to trim  : %d\n',closest.iterations);

fprintf('\n---- Wheel-level state ----\n');
wheelState = closest.wheelState;
wheelTable = table( ...
    {wheelState.name}.', ...
    [wheelState.Fz].', ...
    [wheelState.slipAngleDeg].', ...
    [wheelState.steerDeg].', ...
    [wheelState.Vx].', ...
    [wheelState.FyTire].', ...
    [wheelState.FxBody].', ...
    [wheelState.FyBody].', ...
    [wheelState.yawMoment].', ...
    'VariableNames',{'Wheel','Fz_N','SlipAngle_deg','Steer_deg', ...
    'Vx_mps','FyTire_N','FxBody_N','FyBody_N','YawMoment_Nm'});
disp(wheelTable);

fprintf('\nCombined normalized distance to requested point: %.4f\n\n', ...
    combinedDistance(closestIndex));
