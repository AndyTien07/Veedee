clear; close all; clc;

%% ---------------- USER SETTINGS ----------------
csv_file    = '10 PSI SAE.csv';   % <-- path to your CSV file

mode        = 'lateral';   % 'lateral'      -> fit Fy vs Slip Angle (SA)
                            % 'longitudinal' -> fit Fx vs Slip Ratio (SR)

IA_target   = deg2rad(2);           % [rad] target inclination (camber) angle
IA_tol      = deg2rad(.1);% [rad] tolerance window around IA_target

pureslip_tol = 0.01;       % tolerance for "other" slip channel ~ 0
                            % (SR tolerance for lateral fits, SA tolerance
                            %  [rad] for longitudinal fits)

Fz_bin_tol  = 50;          % [N] points within this range of each other
                            % are grouped into the same load level
Fz_min_pts  = 15;           % minimum points required to fit a load level

poly_order  = 2;            % polynomial order used to fit B(Fz), C(Fz),
                             % D(Fz), E(Fz) trends (2 = quadratic typical)

%% ---------------- STEP 1: Load CSV --------------------------------------
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
FY=-(FY)
IA=abs(IA);

if strcmp(mode, 'lateral')
    x_all = SA; y_all = FY; other_slip = SR;
    xlabel_str = 'Slip Angle \alpha [deg]';
    ylabel_str = 'F_y [N]';
else
    x_all = SR; y_all = FX; other_slip = SA;
    xlabel_str = 'Slip Ratio \kappa [-]';
    ylabel_str = 'F_x [N]';
end

mask_IA = abs(IA - IA_target) <= IA_tol;
mask_pure = abs(other_slip) <= pureslip_tol;
mask_base = mask_IA & mask_pure;

if nnz(mask_base) < Fz_min_pts
    error('Too few points after IA/pure-slip filtering. Widen IA_tol/pureslip_tol.');
end

FZ_f = FZ(mask_base); x_f = x_all(mask_base); y_f = y_all(mask_base);

%% ---------------- STEP 2: Auto-cluster distinct Fz load levels ----------
[FZ_sorted, sortIdx] = sort(FZ_f);
x_sorted = x_f(sortIdx); y_sorted = y_f(sortIdx);

load_levels = {}; % each cell: struct with Fz_nom, x, y
cur_start = 1;
for i = 2:length(FZ_sorted)+1
    if i > length(FZ_sorted) || (FZ_sorted(i) - FZ_sorted(cur_start)) > Fz_bin_tol
        idx = cur_start:(i-1);
        if numel(idx) >= Fz_min_pts
            lvl.Fz_nom = mean(FZ_sorted(idx));
            lvl.x = x_sorted(idx);
            lvl.y = y_sorted(idx);
            load_levels{end+1} = lvl; %#ok<SAGROW>
        end
        cur_start = i;
    end
end

nLevels = numel(load_levels);
if nLevels < 2
    error(['Only found %d distinct load level(s) with >= %d points. ' ...
        'Need >=2 levels to fit load sensitivity. Adjust Fz_bin_tol/Fz_min_pts.'], ...
        nLevels, Fz_min_pts);
end
fprintf('Found %d distinct load levels for fitting.\n', nLevels);

%% ---------------- STEP 3: Local 4-parameter fit per load level ----------
pacejka4 = @(p, x) p(3) .* sin( p(2) .* atan( p(1).*x - p(4).*(p(1).*x - atan(p(1).*x)) ) );
fitopts = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 5000, 'MaxIterations', 2000);

Fz_nom_arr = zeros(nLevels,1);
B_arr = zeros(nLevels,1); C_arr = zeros(nLevels,1);
D_arr = zeros(nLevels,1); E_arr = zeros(nLevels,1);
rmse_arr = zeros(nLevels,1); r2_arr = zeros(nLevels,1);

for k = 1:nLevels
    lvl = load_levels{k};
    [xk, sIdx] = sort(lvl.x); yk = lvl.y(sIdx);

    D0 = max(abs(yk));
    [~, i0] = min(abs(xk));
    win = 15;
    idx = max(i0-win,1):min(i0+win,length(xk));
    p = polyfit(xk(idx), yk(idx), 1);
    BCD0 = p(1);

    n_tail = max(round(0.03*numel(xk)), 5);
    ya_pos = mean(yk(end-n_tail+1:end));
    ya_neg = mean(yk(1:n_tail));
    ya0 = (abs(ya_pos) + abs(ya_neg)) / 2;

    C0 = (2/pi) * asin( min(max(ya0/D0, -1), 1) );
    if C0 <= 0, C0 = 1.3; end
    B0 = BCD0 / (C0 * D0);

    [~, imax] = max(abs(yk));
    xm0 = abs(xk(imax));
    if xm0 < 1e-6, xm0 = max(abs(xk))/2; end
    denom = (B0*xm0 - atan(B0*xm0));
    if abs(denom) > 1e-9
        E0 = (B0*xm0 - tan(pi/(2*C0))) / denom;
    else
        E0 = 0;
    end
    E0 = min(max(E0, -5), 1);

    p0 = [B0, C0, D0, E0];
    lb = [1e-3, 0.1, 0, -10];
    ub = [500,  3.0, 1.5*D0, 5];

    p_fit = lsqcurvefit(pacejka4, p0, xk, yk, lb, ub, fitopts);

    B_arr(k) = p_fit(1); C_arr(k) = p_fit(2);
    D_arr(k) = p_fit(3); E_arr(k) = p_fit(4);
    Fz_nom_arr(k) = lvl.Fz_nom;

    y_fit = pacejka4(p_fit, xk);
    rmse_arr(k) = sqrt(mean((y_fit - yk).^2));
    r2_arr(k) = 1 - sum((yk-y_fit).^2) / sum((yk-mean(yk)).^2);

    load_levels{k}.x_sorted = xk;
    load_levels{k}.y_sorted = yk;
    load_levels{k}.p_fit = p_fit;

    fprintf('Fz=%6.0fN | B=%7.4f C=%6.4f D=%8.2f E=%7.4f | RMSE=%7.2f R^2=%.4f\n', ...
        Fz_nom_arr(k), B_arr(k), C_arr(k), D_arr(k), E_arr(k), rmse_arr(k), r2_arr(k));
end

%% ---------------- STEP 4: Fit B,C,D,E as functions of Fz -----------------
% Polynomial regression captures load sensitivity while keeping the
% classic 4-param functional form usable at ANY Fz within the tested range.
B_poly = polyfit(Fz_nom_arr([1,2,3,4,5]), B_arr([1,2,3,4,5]), poly_order);
C_poly = polyfit(Fz_nom_arr, C_arr, poly_order);
D_poly = polyfit(Fz_nom_arr, D_arr, poly_order);
E_poly = polyfit(Fz_nom_arr([1,2,3,4,5]), E_arr([1,2,3,4,5]), poly_order);

fprintf('\n--- Load-sensitive parameter polynomials (order %d, in Fz [N]) ---\n', poly_order);
fprintf('B(Fz) coeffs: '); disp(B_poly);
fprintf('C(Fz) coeffs: '); disp(C_poly);
fprintf('D(Fz) coeffs: '); disp(D_poly);
fprintf('E(Fz) coeffs: '); disp(E_poly);

% Master load-sensitive model: y = f(x, Fz)
pacejka4_loadsens = @(x, Fz) ...
    polyval(D_poly, Fz) .* sin( polyval(C_poly, Fz) .* atan( ...
        polyval(B_poly, Fz).*x - polyval(E_poly, Fz).*( ...
        polyval(B_poly, Fz).*x - atan(polyval(B_poly, Fz).*x) ) ) );

%% ---------------- STEP 5: Plots ------------------------------------------
% (a) Per-load-level curve fits
figure(); hold on; grid on;
cmap = lines(nLevels);
for k = 1:nLevels
    xk = load_levels{k}.x_sorted; yk = load_levels{k}.y_sorted;
    if strcmp(mode,'lateral'), xp = rad2deg(xk); else, xp = xk; end
    plot(xp, yk, '.', 'Color', cmap(k,:), 'MarkerSize', 5, 'HandleVisibility','off');

    x_smooth = linspace(min(xk), max(xk), 300);
    y_smooth = pacejka4(load_levels{k}.p_fit, x_smooth);
    if strcmp(mode,'lateral'), xps = rad2deg(x_smooth); else, xps = x_smooth; end
    plot(xps, y_smooth, '-', 'Color', cmap(k,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Fz = %.0f N', Fz_nom_arr(k)));
end
xlabel(xlabel_str); ylabel(ylabel_str);
title('Local 4-Parameter Fits per Load Level'); legend('Location','best');

% (b) Parameter vs Fz trends
figure();
Fz_smooth = linspace(min(Fz_nom_arr), max(Fz_nom_arr), 200);
subplot(2,2,1); plot(Fz_nom_arr, B_arr, 'o', Fz_smooth, polyval(B_poly,Fz_smooth), 'r-');
xlabel('F_z [N]'); ylabel('B'); title('Stiffness Factor B(F_z)'); grid on;
subplot(2,2,2); plot(Fz_nom_arr, C_arr, 'o', Fz_smooth, polyval(C_poly,Fz_smooth), 'r-');
xlabel('F_z [N]'); ylabel('C'); title('Shape Factor C(F_z)'); grid on;
subplot(2,2,3); plot(Fz_nom_arr, D_arr, 'o', Fz_smooth, polyval(D_poly,Fz_smooth), 'r-');
xlabel('F_z [N]'); ylabel('D'); title('Peak Factor D(F_z)'); grid on;
subplot(2,2,4); plot(Fz_nom_arr, E_arr, 'o', Fz_smooth, polyval(E_poly,Fz_smooth), 'r-');
xlabel('F_z [N]'); ylabel('E'); title('Curvature Factor E(F_z)'); grid on;
sgtitle('Load-Sensitive Parameter Trends (dots = local fits, line = polynomial)');

% (c) Validate load-sensitive model against raw data across all levels

figure(); hold on; grid on;
for k = 1:nLevels
    xk = load_levels{k}.x_sorted; yk = load_levels{k}.y_sorted;
    [Fz_nom_arr(k)/1000,rad2deg(xk(yk==max(yk))),max(yk)/1000,max(yk)/Fz_nom_arr(k)]
    if strcmp(mode,'lateral'), xp = rad2deg(xk); else, xp = xk; end
    plot(xp, yk, '.', 'Color', cmap(k,:), 'MarkerSize', 5, 'HandleVisibility','off');

    x_smooth = linspace(min(xk), max(xk), 300);
    y_smooth = pacejka4_loadsens(x_smooth, Fz_nom_arr(k));
    if strcmp(mode,'lateral'), xps = rad2deg(x_smooth); else, xps = x_smooth; end
    plot(xps, y_smooth, '--', 'Color', cmap(k,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Fz=%.0fN (poly model)', Fz_nom_arr(k)));
end
xlabel(xlabel_str); ylabel(ylabel_str);
title('Load-Sensitive 4-Param Model vs Data (all levels)'); legend('Location','best');

fprintf('\nDone. Use pacejka4_loadsens(x, Fz) to evaluate the model at any Fz\n');
fprintf('within [%.0f, %.0f] N.\n', min(Fz_nom_arr), max(Fz_nom_arr));
