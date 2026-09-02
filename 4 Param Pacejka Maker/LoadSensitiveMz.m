%% pacejka_mz_4param_load_sensitive.m
% Distills MF6.2-generated or physical tire data into a load-sensitive,
% four-parameter Magic Formula for self-aligning moment:
%
%   Mz(alpha,Fz) = D(Fz)*sin(C(Fz)*atan(q))
%   q = B(Fz)*alpha - E(Fz)*(B(Fz)*alpha - atan(B(Fz)*alpha))
%
% Load sensitivity is imposed through normalized load:
%
%   dfz       = (Fz - Fz0)/Fz0
%   D(Fz)     = R0*Fz*(d0 + d1*dfz + ...)
%   C(Fz)     = c0                         [constant by default]
%   E(Fz)     = e0 + e1*dfz + ...
%   Kza(Fz)   = R0*Fz*(k0 + k1*dfz + ...)
%   B(Fz)     = Kza(Fz)/(C(Fz)*D(Fz))
%
% Consequently, the zero-slip slope is exactly
%
%   dMz/dalpha | alpha=0 = B*C*D = Kza.
%
% Local B,C,D,E fits are used only for outlier trimming and initialization.
% The reported coefficients come from one simultaneous global fit.

clear; close all; clc;

%% ---------------- USER SETTINGS ---------------------------------------
csv_file = '10 PSI SAE.csv';

IA_target = deg2rad(2);       % [rad]
IA_tol = deg2rad(0.5);        % [rad]
pureslip_tol = 0.01;          % [-]
SA_fit_range = deg2rad(12);   % [rad]; use Inf for the full sweep

Fz_bin_tol = 150;             % [N]
Fz_min_pts = 15;

% Reference dimensions. R0 only scales the dimensionless d and k
% coefficients; use the tire's unloaded/free radius.
R0 = 0.20574;                   % [m] <-- replace with the correct tire radius
Fz0_user = NaN;               % [N]; NaN selects median fitted load level

% Recommended minimum load-sensitive model:
D_order = 1;                  % D/(R0*Fz) = d0 + d1*dfz
C_order = 0;                  % C = c0; set to 1 only if clearly justified
E_order = 1;                  % E = e0 + e1*dfz
K_order = 1;                  % Kza/(R0*Fz) = k0 + k1*dfz

%% ---------------- OUTLIER REJECTION SETTINGS --------------------------
outlier_enable = true;
outlier_prefit_hampel = true;
hampel_window = 15;
hampel_nsigma = 3;

outlier_iters = 3;
outlier_resid_nsigma = 2.5;
outlier_min_keep_frac = 0.70;

%% ---------------- STEP 1: Load and filter CSV --------------------------
opts = detectImportOptions(csv_file);
T = readtable(csv_file, opts);

required = {'SA','SR','IA','FZ','MZ'};
missing = required(~ismember(required, T.Properties.VariableNames));
if ~isempty(missing)
    error('Missing required CSV variables: %s', strjoin(missing, ', '));
end

SA = T.SA;                    % [rad]
SR = T.SR;                    % [-]
IA = T.IA;                    % [rad]
FZ = abs(T.FZ);               % [N], positive compression magnitude
MZ = T.MZ;                    % [Nm]

finite_mask = isfinite(SA) & isfinite(SR) & isfinite(IA) & ...
              isfinite(FZ) & isfinite(MZ);
mask_IA = abs(IA - IA_target) <= IA_tol;
mask_pure = abs(SR) <= pureslip_tol;
mask_range = abs(SA) <= SA_fit_range;
mask_base = finite_mask & mask_IA & mask_pure & mask_range;

if nnz(mask_base) < Fz_min_pts
    error(['Too few points after IA, pure-slip, and SA-range filtering. ' ...
           'Widen IA_tol, pureslip_tol, or SA_fit_range.']);
end

FZ_f = FZ(mask_base);
x_f = SA(mask_base);
y_f = MZ(mask_base);

%% ---------------- STEP 1b: Gross-glitch Hampel screen -----------------
% This is only a coarse pre-screen. Residual trimming is performed again
% within each load level.
n_before_hampel = numel(x_f);
if outlier_enable && outlier_prefit_hampel
    [x_sort_h, ord_h] = sort(x_f);
    y_sort_h = y_f(ord_h);
    FZ_sort_h = FZ_f(ord_h);

    [~, hampel_outlier_mask] = hampel(y_sort_h, hampel_window, hampel_nsigma);

    x_f = x_sort_h(~hampel_outlier_mask);
    y_f = y_sort_h(~hampel_outlier_mask);
    FZ_f = FZ_sort_h(~hampel_outlier_mask);

    fprintf('Hampel pre-filter removed %d/%d points (%.1f%%).\n', ...
        nnz(hampel_outlier_mask), n_before_hampel, ...
        100*nnz(hampel_outlier_mask)/n_before_hampel);
end

%% ---------------- STEP 2: Cluster distinct Fz levels -------------------
[FZ_sorted, sortIdx] = sort(FZ_f);
x_sorted = x_f(sortIdx);
y_sorted = y_f(sortIdx);

load_levels = {};
cur_start = 1;
for i = 2:(length(FZ_sorted) + 1)
    if i > length(FZ_sorted) || ...
            (FZ_sorted(i) - FZ_sorted(cur_start)) > Fz_bin_tol
        idx = cur_start:(i - 1);
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
max_order = max([D_order, C_order, E_order, K_order]);
if nLevels < max_order + 1
    error(['Found %d usable load levels, but polynomial order %d requires ' ...
           'at least %d levels.'], nLevels, max_order, max_order + 1);
end
if nLevels < 2
    error('At least two distinct load levels are required.');
end
fprintf('Found %d distinct load levels for the Mz fit.\n', nLevels);

%% ---------------- STEP 3: Local fits for cleaning/seeding -------------
pacejka4 = @(p,x) p(3).*sin(p(2).*atan( ...
    p(1).*x - p(4).*(p(1).*x - atan(p(1).*x))));

fit_opts_local = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 5000, 'MaxIterations', 2000);

B_loc = zeros(nLevels,1);
C_loc = zeros(nLevels,1);
D_loc = zeros(nLevels,1);
E_loc = zeros(nLevels,1);
K_loc = zeros(nLevels,1);
Fz_loc = zeros(nLevels,1);

x_clean_all = [];
y_clean_all = [];
Fz_clean_all = [];
removed_x = [];
removed_y = [];
removed_c = [];

figure('Name','Local per-load Mz fits'); hold on;
colors = lines(nLevels);

for k = 1:nLevels
    x = load_levels{k}.x(:);
    y = load_levels{k}.y(:);
    Fz_nom = load_levels{k}.Fz_nom;
    n0 = numel(x);

    % Canonical local parameterization uses D>0. The sign of the initial
    % slope is carried by B, eliminating the equivalent (-B,-D) solution.
    D0 = max(max(abs(y)), 1e-6);
    n_slope = min(n0, max(5, ceil(0.20*n0)));
    [~, near_zero] = sort(abs(x), 'ascend');
    iz = near_zero(1:n_slope);
    slope0 = x(iz) \ y(iz);       % through-origin zero-slip slope estimate
    if ~isfinite(slope0)
        slope0 = 0;
    end

    C0 = 1.30;
    B0 = slope0/(C0*D0);
    if ~isfinite(B0) || abs(B0) < 1e-6
        B0 = sign(sum(x.*y) + eps)*10;
    end
    B0 = min(max(B0, -2000), 2000);
    E0 = 0;

    lb_local = [-2000, 0.50, max(1e-9, 0.01*D0), -5.0];
    ub_local = [ 2000, 2.50, max(1e-6, 10.0*D0),  2.0];
    p0 = min(max([B0, C0, D0, E0], lb_local), ub_local);

    keep = true(size(x));
    for iter = 1:outlier_iters
        xk = x(keep);
        yk = y(keep);
        try
            p_iter = lsqcurvefit(pacejka4, p0, xk, yk, ...
                lb_local, ub_local, fit_opts_local);
        catch
            break;
        end

        resid = y - pacejka4(p_iter, x);
        r_keep = resid(keep);
        sigma_r = 1.4826*mad(r_keep, 1);
        if ~isfinite(sigma_r) || sigma_r <= eps
            sigma_r = std(r_keep);
        end
        if ~isfinite(sigma_r) || sigma_r <= eps
            break;
        end

        new_keep = keep;
        if outlier_enable
            new_keep = keep & (abs(resid - median(r_keep)) <= ...
                outlier_resid_nsigma*sigma_r);
        end

        if nnz(new_keep) < max(4, ceil(outlier_min_keep_frac*n0))
            break;
        end
        if isequal(new_keep, keep)
            break;
        end
        keep = new_keep;
        p0 = p_iter;
    end

    xk = x(keep);
    yk = y(keep);
    p_final = lsqcurvefit(pacejka4, p0, xk, yk, ...
        lb_local, ub_local, fit_opts_local);

    B_loc(k) = p_final(1);
    C_loc(k) = p_final(2);
    D_loc(k) = p_final(3);
    E_loc(k) = p_final(4);
    K_loc(k) = prod(p_final(1:3));
    Fz_loc(k) = Fz_nom;

    x_clean_all = [x_clean_all; xk]; %#ok<AGROW>
    y_clean_all = [y_clean_all; yk]; %#ok<AGROW>
    Fz_clean_all = [Fz_clean_all; Fz_nom*ones(size(xk))]; %#ok<AGROW>

    removed_x = [removed_x; x(~keep)]; %#ok<AGROW>
    removed_y = [removed_y; y(~keep)]; %#ok<AGROW>
    removed_c = [removed_c; repmat(colors(k,:), nnz(~keep), 1)]; %#ok<AGROW>

    n_removed = n0 - nnz(keep);
    if n_removed > 0
        fprintf('  Fz=%6.0f N: removed %d/%d outliers (%.1f%%)\n', ...
            Fz_nom, n_removed, n0, 100*n_removed/n0);
    end

    x_plot = linspace(min(xk), max(xk), 200)';
    plot(rad2deg(xk), yk, 'o', 'Color', colors(k,:), ...
        'MarkerSize', 4, 'HandleVisibility','off');
    plot(rad2deg(x_plot), pacejka4(p_final, x_plot), '-', ...
        'Color', colors(k,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Fz=%.0f N local', Fz_nom));
end

if ~isempty(removed_x)
    scatter(rad2deg(removed_x), removed_y, 40, removed_c, 'x', ...
        'LineWidth', 1.5, 'HandleVisibility','off');
end
xlabel('Slip angle \alpha [deg]');
ylabel('M_z [Nm]');
title('Local per-load fits used for cleaning and initialization');
legend('show','Location','best'); grid on;

%% ---------------- STEP 4: Global load-sensitive fit --------------------
if isnan(Fz0_user)
    Fz0 = median(Fz_loc);
else
    Fz0 = Fz0_user;
end
if ~(isfinite(Fz0) && Fz0 > 0)
    error('Fz0 must be finite and positive.');
end
if ~(isfinite(R0) && R0 > 0)
    error('R0 must be finite and positive.');
end

dfz_loc = (Fz_loc - Fz0)/Fz0;

% Convert local fits to the proposed governing quantities.
d_norm_loc = D_loc./(R0*Fz_loc);
k_norm_loc = K_loc./(R0*Fz_loc);

pd0 = polyfit_ascending(dfz_loc, d_norm_loc, D_order);
pc0 = polyfit_ascending(dfz_loc, C_loc, C_order);
pe0 = polyfit_ascending(dfz_loc, E_loc, E_order);
pk0 = polyfit_ascending(dfz_loc, k_norm_loc, K_order);

nD = D_order + 1;
nC = C_order + 1;
nE = E_order + 1;
nK = K_order + 1;
p0_global = [pd0, pc0, pe0, pk0];

% Coefficient bounds remove the worst non-identifiable local solutions.
% Pointwise trend checks after fitting still verify D>0 and C>0 over the
% measured load range.
d_scale = max([median(abs(d_norm_loc)), abs(pd0(1)), 1e-6]);
k_scale = max([median(abs(k_norm_loc)), abs(pk0(1)), 1e-6]);

lb_d = -50*d_scale*ones(1,nD);
ub_d =  50*d_scale*ones(1,nD);
lb_d(1) = max(1e-10, 0.01*d_scale);
ub_d(1) = 50*d_scale;

lb_c = -3*ones(1,nC);
ub_c =  3*ones(1,nC);
lb_c(1) = 0.30;
ub_c(1) = 3.00;

lb_e = -10*ones(1,nE);
ub_e =  10*ones(1,nE);
lb_e(1) = -5.0;
ub_e(1) =  2.0;

lb_k = -50*k_scale*ones(1,nK);
ub_k =  50*k_scale*ones(1,nK);

lb = [lb_d, lb_c, lb_e, lb_k];
ub = [ub_d, ub_c, ub_e, ub_k];
p0_global = min(max(p0_global, lb), ub);

X_global = [x_clean_all, Fz_clean_all];
model_global = @(p,X) eval_load_sensitive_mz( ...
    p, X, nD, nC, nE, nK, R0, Fz0);

fit_opts_global = optimoptions('lsqcurvefit', 'Display', 'off', ...
    'MaxFunctionEvaluations', 30000, 'MaxIterations', 8000, ...
    'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);

p_global = lsqcurvefit(model_global, p0_global, X_global, y_clean_all, ...
    lb, ub, fit_opts_global);

id = 0;
pDnorm = p_global(id + (1:nD)); id = id + nD;
pC = p_global(id + (1:nC));     id = id + nC;
pE = p_global(id + (1:nE));     id = id + nE;
pKnorm = p_global(id + (1:nK));

D_of_Fz = @(Fz) R0.*Fz.*polyval_ascending(pDnorm, (Fz-Fz0)./Fz0);
C_of_Fz = @(Fz) polyval_ascending(pC, (Fz-Fz0)./Fz0);
E_of_Fz = @(Fz) polyval_ascending(pE, (Fz-Fz0)./Fz0);
K_of_Fz = @(Fz) R0.*Fz.*polyval_ascending(pKnorm, (Fz-Fz0)./Fz0);
B_of_Fz = @(Fz) K_of_Fz(Fz)./(C_of_Fz(Fz).*D_of_Fz(Fz));
pacejka4_loadsens = @(alpha,Fz) eval_load_sensitive_mz( ...
    p_global, [alpha(:), Fz(:)], nD, nC, nE, nK, R0, Fz0);

fprintf('\nGlobal load-sensitive Mz fit complete.\n');
fprintf('R0  = %.6g m\n', R0);
fprintf('Fz0 = %.6g N\n', Fz0);
fprintf('d coefficients [d0 d1 ...]: '); disp(pDnorm);
fprintf('c coefficients [c0 c1 ...]: '); disp(pC);
fprintf('e coefficients [e0 e1 ...]: '); disp(pE);
fprintf('k coefficients [k0 k1 ...]: '); disp(pKnorm);

%% ---------------- STEP 5: Trend plots and physical checks -------------
Fz_plot = linspace(min(Fz_loc), max(Fz_loc), 300)';
D_plot = D_of_Fz(Fz_plot);
C_plot = C_of_Fz(Fz_plot);
E_plot = E_of_Fz(Fz_plot);
K_plot = K_of_Fz(Fz_plot);
B_plot = B_of_Fz(Fz_plot);

if any(D_plot <= 0)
    warning(['D(Fz) became non-positive inside the measured load range. ' ...
             'Reduce D_order or tighten the D coefficient bounds.']);
end
if any(C_plot <= 0)
    warning(['C(Fz) became non-positive inside the measured load range. ' ...
             'Use C_order=0 or tighten the C coefficient bounds.']);
end
if any(~isfinite(B_plot)) || max(abs(B_plot)) > 5000
    warning(['Derived B(Fz) is singular or extremely large. Check D/C ' ...
             'trends, sign conventions, and the fitted SA range.']);
end

figure('Name','Mz load-sensitive parameter trends');
subplot(2,3,1);
plot(Fz_loc, B_loc, 'o', 'DisplayName','local seed'); hold on;
plot(Fz_plot, B_plot, '-', 'LineWidth',1.6, 'DisplayName','global');
xlabel('F_z [N]'); ylabel('B'); title('Derived B=K_{z\alpha}/(CD)');
grid on; legend show;

subplot(2,3,2);
plot(Fz_loc, C_loc, 'o', 'DisplayName','local seed'); hold on;
plot(Fz_plot, C_plot, '-', 'LineWidth',1.6, 'DisplayName','global');
xlabel('F_z [N]'); ylabel('C'); title('C(F_z)'); grid on; legend show;

subplot(2,3,3);
plot(Fz_loc, D_loc, 'o', 'DisplayName','local seed'); hold on;
plot(Fz_plot, D_plot, '-', 'LineWidth',1.6, 'DisplayName','global');
xlabel('F_z [N]'); ylabel('D [Nm]'); title('D(F_z)'); grid on; legend show;

subplot(2,3,4);
plot(Fz_loc, E_loc, 'o', 'DisplayName','local seed'); hold on;
plot(Fz_plot, E_plot, '-', 'LineWidth',1.6, 'DisplayName','global');
xlabel('F_z [N]'); ylabel('E'); title('E(F_z)'); grid on; legend show;

subplot(2,3,5);
plot(Fz_loc, K_loc, 'o', 'DisplayName','local BCD'); hold on;
plot(Fz_plot, K_plot, '-', 'LineWidth',1.6, 'DisplayName','global');
xlabel('F_z [N]'); ylabel('K_{z\alpha} [Nm/rad]');
title('Zero-slip stiffness'); grid on; legend show;

subplot(2,3,6);
plot(Fz_loc, d_norm_loc, 'o', 'DisplayName','local D/(R_0F_z)'); hold on;
plot(Fz_plot, D_plot./(R0*Fz_plot), '-', 'LineWidth',1.6, ...
    'DisplayName','global');
xlabel('F_z [N]'); ylabel('D/(R_0F_z)');
title('Normalized peak moment'); grid on; legend show;

%% ---------------- STEP 6: Validation ----------------------------------
y_pred_all = pacejka4_loadsens(x_clean_all, Fz_clean_all);
resid_all = y_clean_all - y_pred_all;
rmse = sqrt(mean(resid_all.^2));
mae = mean(abs(resid_all));
ss_res = sum(resid_all.^2);
ss_tot = sum((y_clean_all - mean(y_clean_all)).^2);
r2 = 1 - ss_res/max(ss_tot, eps);

fprintf('\nGlobal model validation on cleaned data:\n');
fprintf('  RMSE = %.4f Nm\n', rmse);
fprintf('  MAE  = %.4f Nm\n', mae);
fprintf('  R^2  = %.6f\n', r2);

figure('Name','Mz global-model validation'); hold on;
for k = 1:nLevels
    Fz_nom = Fz_loc(k);
    level_mask = abs(Fz_clean_all - Fz_nom) < 1e-9;
    x_lvl = x_clean_all(level_mask);
    y_lvl = y_clean_all(level_mask);
    x_plot = linspace(min(x_lvl), max(x_lvl), 250)';
    Fz_plot_level = Fz_nom*ones(size(x_plot));

    plot(rad2deg(x_lvl), y_lvl, 'o', 'Color', colors(k,:), ...
        'MarkerSize',4, 'HandleVisibility','off');
    plot(rad2deg(x_plot), pacejka4_loadsens(x_plot, Fz_plot_level), '-', ...
        'Color', colors(k,:), 'LineWidth',1.8, ...
        'DisplayName', sprintf('Fz=%.0f N global', Fz_nom));
end
xlabel('Slip angle \alpha [deg]');
ylabel('M_z [Nm]');
title(sprintf('Load-sensitive M_z model | RMSE=%.3f Nm, R^2=%.4f', ...
    rmse, r2));
legend('show','Location','best'); grid on;

figure('Name','Mz residual validation');
subplot(1,2,1);
scatter(y_pred_all, resid_all, 14, Fz_clean_all, 'filled');
yline(0,'k--'); grid on; colorbar;
xlabel('Predicted M_z [Nm]'); ylabel('Residual [Nm]');
title('Residual vs prediction');

subplot(1,2,2);
scatter(rad2deg(x_clean_all), resid_all, 14, Fz_clean_all, 'filled');
yline(0,'k--'); grid on; colorbar;
xlabel('Slip angle \alpha [deg]'); ylabel('Residual [Nm]');
title('Residual vs slip angle');

%% ---------------- Local helper functions ------------------------------
function y = eval_load_sensitive_mz(p, X, nD, nC, nE, nK, R0, Fz0)
    alpha = X(:,1);
    Fz = X(:,2);
    dfz = (Fz - Fz0)./Fz0;

    id = 0;
    pd = p(id + (1:nD)); id = id + nD;
    pc = p(id + (1:nC)); id = id + nC;
    pe = p(id + (1:nE)); id = id + nE;
    pk = p(id + (1:nK));

    D = R0.*Fz.*polyval_ascending(pd, dfz);
    C = polyval_ascending(pc, dfz);
    E = polyval_ascending(pe, dfz);
    K = R0.*Fz.*polyval_ascending(pk, dfz);

    denominator = C.*D;
    tiny = 1e-12*max(1, max(abs(denominator)));
    bad = abs(denominator) < tiny;
    denominator(bad) = tiny.*sign(denominator(bad) + eps);
    B = K./denominator;

    z = B.*alpha;
    q = z - E.*(z - atan(z));
    y = D.*sin(C.*atan(q));
end

function y = polyval_ascending(a, x)
    y = zeros(size(x));
    for j = 1:numel(a)
        y = y + a(j).*x.^(j-1);
    end
end

function a = polyfit_ascending(x, y, order)
    p_descending = polyfit(x(:), y(:), order);
    a = fliplr(p_descending);
end
