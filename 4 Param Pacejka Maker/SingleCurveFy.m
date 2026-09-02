%% pacejka_6p2_to_4param.m
% Distills real MF6.2-generated (or physical) tire test data, loaded from a
% CSV file, down to the classic 4-parameter Magic Formula:
%
%   y(x) = D * sin( C * atan( B*x - E*(B*x - atan(B*x)) ) )
%
% Expects a CSV with columns (row 1 = header names, row 2 = units), e.g.
% TTC-style data:
%
%   Elapsed Time (ET)      [s]
%   Slip Angle (SA)        [rad]
%   Slip Ratio (SR)        [-]
%   Inclination Angle (IA) [rad]
%   Normal Load (FZ)       [N]
%   Speed (V)              [m/s]
%   Pressure (P)           [Pa]
%   Effective Rolling Radius (REIn) [m]
%   Turnslip (TS)          [1/m]
%   Longitudinal Force (FX)[N]
%   Lateral Force (FY)     [N]
%   Overturning Moment (MX)[Nm]
%   Rolling Resistance Moment (MY) [Nm]
%   Self-Aligning Moment (MZ)      [Nm]
%   Calculated Effective Rolling Radius (REOut) [m]
%   Loaded Radius (RL)     [m]
%   Pneumatic Trailing Arm (T)     [m]
%   Pneumatic Scrub Radius (TX)    [m]
%
% Workflow:
%   1) Load CSV, skip units row
%   2) Filter to a target Fz, IA (camber), and pure-slip condition
%   3) Extract initial B,C,D,E analytically from the filtered data
%   4) Refine with lsqcurvefit (nonlinear least squares)
%   5) Plot comparison + report goodness of fit
%
% Author: generated for atien07@uw.edu

%% ---------------- USER SETTINGS ----------------
csv_file    = '12 PSI.csv';   % <-- path to your CSV file

mode        = 'lateral';   % 'lateral'      -> fit Fy vs Slip Angle (SA)
                            % 'longitudinal' -> fit Fx vs Slip Ratio (SR)

Fz_target   = 200;        % [N] target normal load to distill at
Fz_tol      = 100;         % [N] tolerance window around Fz_target

IA_target   = 0;           % [rad] target inclination (camber) angle
IA_tol      = deg2rad(1);% [rad] tolerance window around IA_target

pureslip_tol = 0.01;       % tolerance for "other" slip channel ~ 0
                            % (SR tolerance for lateral fits, SA tolerance
                            %  [rad] for longitudinal fits)

%% ---------------- STEP 1: Load CSV --------------------------------------
% detectImportOptions handles the two-row header (names + units) cleanly
% by reading names from row 1 and treating row 2 (units) as data, which we
% then strip off manually.
opts = detectImportOptions(csv_file);
T = readtable(csv_file, opts);

% Standardize expected column names -> variables
% (Adjust these strings if your CSV header text differs slightly.)
varNames = T.Properties.VariableNames;
getCol = @(pat) T{:, find(contains(varNames, pat, 'IgnoreCase', true), 1)};
SA = T.SA;      % [rad]
SR = T.SR;      % [-]
IA = T.IA;     % [rad]
FZ = T.FZ;       % [N]  (matches "Normal Load (FZ)")
FX = T.FX;      % [N]
FY = T.FY;         % [N]

% FZ in TTC data is usually stored as negative (compression convention).
% Normalize to positive magnitude for filtering/fitting.
FZ = abs(FZ);
IA=abs(IA);
%% ---------------- STEP 2: Filter to target operating condition ---------
mask_Fz = abs(FZ - Fz_target) <= Fz_tol;
mask_IA = abs(IA - IA_target) <= IA_tol;

if strcmp(mode, 'lateral')
    mask_pure = abs(SR) <= pureslip_tol;
    x_data = SA(mask_Fz & mask_IA & mask_pure);
    y_data = FY(mask_Fz & mask_IA & mask_pure);
    xlabel_str = 'Slip Angle \alpha [deg]';
    ylabel_str = 'F_y [N]';
else
    mask_pure = abs(SA) <= pureslip_tol;
    x_data = SR(mask_Fz & mask_IA & mask_pure);
    y_data = FX(mask_Fz & mask_IA & mask_pure);
    xlabel_str = 'Slip Ratio \kappa [-]';
    ylabel_str = 'F_x [N]';
end

if numel(x_data) < 10
    error(['Not enough points matched the filter (Fz=%.0f\x00b1%.0fN, ' ...
        'IA=%.2f\x00b1%.2f deg). Found %d points. Widen tolerances or ' ...
        'check column names.'], Fz_target, Fz_tol, rad2deg(IA_target), ...
        rad2deg(IA_tol), numel(x_data));
end

% Sort by x for clean slope/asymptote extraction
[x_data, sortIdx] = sort(x_data);
y_data = y_data(sortIdx);

fprintf('Filtered dataset: %d points at Fz=%.0f\xb1%.0fN, IA=%.2f deg\n', ...
    numel(x_data), Fz_target, Fz_tol, rad2deg(IA_target));

%% ---------------- STEP 3: Analytical initial guess ----------------------
D0 = max(abs(y_data));

% Slope at origin (BCD) via linear fit of points nearest x=0
[~, i0] = min(abs(x_data));
win = 15; % points on each side, clipped to data bounds
idx = max(i0-win,1):min(i0+win,length(x_data));
p = polyfit(x_data(idx), y_data(idx), 1);
BCD0 = p(1);

% Asymptotic value from extremes of the sweep
n_tail = max(round(0.03*numel(x_data)), 5);
ya_pos = mean(y_data(end-n_tail+1:end));
ya_neg = mean(y_data(1:n_tail));
ya0 = (abs(ya_pos) + abs(ya_neg)) / 2;

C0 = (2/pi) * asin( min(max(ya0/D0, -1), 1) );
if C0 <= 0, C0 = 1.3; end % fallback if data doesn't show clear asymptote

B0 = BCD0 / (C0 * D0);

% Location of peak |y| to estimate E
[~, imax] = max(abs(y_data));
xm0 = abs(x_data(imax));
if xm0 < 1e-6, xm0 = max(abs(x_data))/2; end

denom = (B0*xm0 - atan(B0*xm0));
if abs(denom) > 1e-9
    E0 = (B0*xm0 - tan(pi/(2*C0))) / denom;
else
    E0 = 0;
end
E0 = min(max(E0, -5), 1); % clamp to sane range

p0 = [B0, C0, D0, E0];
fprintf('Initial guess:  B=%.4f  C=%.4f  D=%.4f  E=%.4f\n', p0);

%% ---------------- STEP 4: Nonlinear least-squares refinement -----------
pacejka4 = @(p, x) p(3) .* sin( p(2) .* atan( p(1).*x - p(4).*(p(1).*x - atan(p(1).*x)) ) );

fitopts = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 5000, 'MaxIterations', 2000);

lb = [1e-3, 0.1, 0, -10];
ub = [500,  3.0, 1.5*D0, 5];

[p_fit, resnorm] = lsqcurvefit(pacejka4, p0, x_data, y_data, lb, ub, fitopts);

B = p_fit(1); C = p_fit(2); D = p_fit(3); E = p_fit(4);

fprintf('\nFitted 4-parameter Pacejka (at Fz=%.0fN, IA=%.2f deg):\n', ...
    Fz_target, rad2deg(IA_target));
fprintf('  B = %.5f\n  C = %.5f\n  D = %.5f\n  E = %.5f\n', B, C, D, E);

y_fit = pacejka4(p_fit, x_data);
rmse = sqrt(mean((y_fit - y_data).^2));
r2 = 1 - sum((y_data-y_fit).^2) / sum((y_data-mean(y_data)).^2);
fprintf('\nRMSE = %.3f   R^2 = %.5f   (resnorm = %.3e)\n', rmse, r2, resnorm);

%% ---------------- STEP 5: Plot comparison -------------------------------
figure('Color','k');
if strcmp(mode,'lateral')
    xplot = rad2deg(x_data);
else
    xplot = x_data;
end
plot(xplot, y_data, 'b.', 'MarkerSize', 6, 'DisplayName', 'Measured / MF6.2 data'); hold on;

x_smooth = linspace(min(x_data), max(x_data), 400);
y_smooth = pacejka4(p_fit, x_smooth);
if strcmp(mode,'lateral')
    xplot_smooth = rad2deg(x_smooth);
else
    xplot_smooth = x_smooth;
end
plot(xplot_smooth, y_smooth, 'r-', 'LineWidth', 2, 'DisplayName', '4-Param Fit');

grid on; legend('Location','best');
xlabel(xlabel_str); ylabel(ylabel_str);
title(sprintf('Pacejka 6.2 \\rightarrow 4-Param Distillation @ F_z=%.0fN, IA=%.1f\\circ', ...
    Fz_target, rad2deg(IA_target)));
