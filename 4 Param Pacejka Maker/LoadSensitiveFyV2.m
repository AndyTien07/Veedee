%% pacejka_fy_4param_load_sensitive.m
% Fits a load-sensitive four-parameter Magic Formula for lateral force:
%
%   Fy(alpha,Fz) = D(Fz)*sin(C(Fz)*atan(q))
%   q            = B(Fz)*alpha - E(Fz)*(B(Fz)*alpha-atan(B(Fz)*alpha))
%
% using the load laws
%
%   dfz    = (Fz-Fz0)/Fz0
%   D(Fz)  = Fz*(d0+d1*dfz+...)
%   C(Fz)  = c0+c1*dfz+...
%   E(Fz)  = e0+e1*dfz+...
%   Ky(Fz) = Fz*(k0+k1*dfz+...)
%   B(Fz)  = Ky(Fz)/(C(Fz)*D(Fz))
%
% Therefore the zero-slip slope is exactly
%
%   dFy/dalpha | alpha=0 = B*C*D = Ky.
%
% Local B,C,D,E fits are used only for outlier removal and initialization.
% The reported coefficients are obtained from one simultaneous global fit.
% Friction scaling is applied during evaluation, not fitting:
%
%   Fy_scaled = muy*Fy_reference.
%
% The saved model is compatible with eval_load_sensitive_fy.m.

clear; close all; clc;

%% ---------------- USER SETTINGS ---------------------------------------
csv_file = '12 PSI SAE.csv';

IA_target = deg2rad(4);       % [rad]
IA_tol = deg2rad(0.5);        % [rad]
pureslip_tol = 0.01;          % [-]
SA_fit_range = deg2rad(12);   % [rad]; use Inf for the full sweep

Fz_bin_tol = 150;             % [N]
Fz_min_pts = 15;
Fz0_user = NaN;               % [N]; NaN selects median fitted load level

% Recommended minimum load-sensitive model.
D_order = 1;                  % D/Fz = d0+d1*dfz
C_order = 0;                  % C = c0
E_order = 1;                  % E = e0+e1*dfz
K_order = 1;                  % Ky/Fz = k0+k1*dfz

output_mat_file = 'pacejka_fy_4param_model.mat';

%% ---------------- OUTLIER SETTINGS ------------------------------------
outlier_enable = true;
outlier_prefit_hampel = true;
hampel_window = 15;
hampel_nsigma = 3;

outlier_iters = 3;
outlier_resid_nsigma = 2.5;
outlier_min_keep_frac = 0.70;

%% ---------------- STEP 1: LOAD AND FILTER CSV -------------------------
opts = detectImportOptions(csv_file);
T = readtable(csv_file, opts);

required = {'SA','SR','IA','FZ','FY'};
missing = required(~ismember(required, T.Properties.VariableNames));
if ~isempty(missing)
    error('Missing required CSV variables: %s', strjoin(missing, ', '));
end

SA = T.SA;                    % [rad]
SR = T.SR;                    % [-]
IA = T.IA;                    % [rad]
FZ = abs(T.FZ);               % [N], positive compression magnitude
FY = T.FY;                    % [N]

finite_mask = isfinite(SA) & isfinite(SR) & isfinite(IA) & ...
              isfinite(FZ) & isfinite(FY);
mask_IA = abs(IA-IA_target) <= IA_tol;
mask_pure = abs(SR) <= pureslip_tol;
mask_range = abs(SA) <= SA_fit_range;
mask_base = finite_mask & mask_IA & mask_pure & mask_range;

if nnz(mask_base) < Fz_min_pts
    error(['Too few points after IA, pure-slip, and SA-range filtering. ' ...
           'Widen IA_tol, pureslip_tol, or SA_fit_range.']);
end

FZ_f = FZ(mask_base);
x_f = SA(mask_base);
y_f = FY(mask_base);

%% ---------------- STEP 2: CLUSTER DISTINCT LOAD LEVELS ----------------
[FZ_sorted, sort_idx] = sort(FZ_f);
x_sorted = x_f(sort_idx);
y_sorted = y_f(sort_idx);

load_levels = {};
cur_start = 1;
for i = 2:(numel(FZ_sorted)+1)
    end_of_level = i > numel(FZ_sorted);
    if ~end_of_level
        end_of_level = (FZ_sorted(i)-FZ_sorted(cur_start)) > Fz_bin_tol;
    end

    if end_of_level
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
max_order = max([D_order,C_order,E_order,K_order]);
if nLevels < max_order+1
    error(['Found %d usable load levels, but polynomial order %d requires ' ...
           'at least %d levels.'], nLevels, max_order, max_order+1);
end
if nLevels < 2
    error('At least two distinct load levels are required.');
end
fprintf('Found %d distinct load levels for the Fy fit.\n', nLevels);

%% ---------------- STEP 3: LOCAL FITS FOR CLEANING AND SEEDING ---------
pacejka4 = @(p,x) p(3).*sin(p(2).*atan( ...
    p(1).*x-p(4).*(p(1).*x-atan(p(1).*x))));

fit_opts_local = optimoptions('lsqcurvefit', ...
    'Display','off', ...
    'MaxFunctionEvaluations',5000, ...
    'MaxIterations',2000);

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

figure('Name','Local per-load Fy fits'); hold on;
colors = lines(nLevels);

for j = 1:nLevels
    x = load_levels{j}.x(:);
    y = load_levels{j}.y(:);
    Fz_nom = load_levels{j}.Fz_nom;

    % Hampel screening is performed within each load level so that force
    % differences between load levels are not mistaken for outliers.
    if outlier_enable && outlier_prefit_hampel && numel(y) >= hampel_window
        [x, order] = sort(x);
        y = y(order);
        [~, is_hampel_outlier] = hampel(y,hampel_window,hampel_nsigma);
        x_hampel_removed = x(is_hampel_outlier);
        y_hampel_removed = y(is_hampel_outlier);
        x = x(~is_hampel_outlier);
        y = y(~is_hampel_outlier);
    else
        x_hampel_removed = [];
        y_hampel_removed = [];
    end

    n0 = numel(x);
    if n0 < Fz_min_pts
        error('Too few points remain at Fz = %.0f N after Hampel filtering.',Fz_nom);
    end

    % Canonical local parameterization: D>0, while B carries the sign of
    % the initial Fy slope. This avoids the equivalent (-B,-D) solution.
    D0 = max(max(abs(y)),1e-6);
    n_slope = min(n0,max(5,ceil(0.20*n0)));
    [~,near_zero] = sort(abs(x),'ascend');
    iz = near_zero(1:n_slope);
    slope0 = x(iz)\y(iz);     % through-origin estimate near alpha=0
    if ~isfinite(slope0)
        slope0 = 0;
    end

    C0 = 1.30;
    B0 = slope0/(C0*D0);
    if ~isfinite(B0) || abs(B0) < 1e-6
        B0 = sign(sum(x.*y)+eps)*10;
    end
    B0 = min(max(B0,-2000),2000);
    E0 = 0;

    lb_local = [-2000,0.50,max(1e-9,0.01*D0),-5.0];
    ub_local = [ 2000,2.50,max(1e-6,10.0*D0), 2.0];
    p0 = min(max([B0,C0,D0,E0],lb_local),ub_local);

    keep = true(size(x));
    for iter = 1:outlier_iters
        try
            p_iter = lsqcurvefit(pacejka4,p0,x(keep),y(keep), ...
                lb_local,ub_local,fit_opts_local);
        catch
            break;
        end

        resid = y-pacejka4(p_iter,x);
        r_keep = resid(keep);
        sigma_r = 1.4826*mad(r_keep,1);
        if ~isfinite(sigma_r) || sigma_r <= eps
            sigma_r = std(r_keep);
        end
        if ~isfinite(sigma_r) || sigma_r <= eps
            break;
        end

        if outlier_enable
            center_r = median(r_keep);
            new_keep = keep & ...
                (abs(resid-center_r) <= outlier_resid_nsigma*sigma_r);
        else
            new_keep = keep;
        end

        min_keep = max(4,ceil(outlier_min_keep_frac*n0));
        if nnz(new_keep) < min_keep || isequal(new_keep,keep)
            break;
        end

        keep = new_keep;
        p0 = p_iter;
    end

    xk = x(keep);
    yk = y(keep);
    p_final = lsqcurvefit(pacejka4,p0,xk,yk, ...
        lb_local,ub_local,fit_opts_local);

    B_loc(j) = p_final(1);
    C_loc(j) = p_final(2);
    D_loc(j) = p_final(3);
    E_loc(j) = p_final(4);
    K_loc(j) = prod(p_final(1:3));
    Fz_loc(j) = Fz_nom;

    x_clean_all = [x_clean_all;xk]; %#ok<AGROW>
    y_clean_all = [y_clean_all;yk]; %#ok<AGROW>
    Fz_clean_all = [Fz_clean_all;Fz_nom*ones(size(xk))]; %#ok<AGROW>

    x_removed = [x_hampel_removed;x(~keep)];
    y_removed = [y_hampel_removed;y(~keep)];
    removed_x = [removed_x;x_removed]; %#ok<AGROW>
    removed_y = [removed_y;y_removed]; %#ok<AGROW>
    removed_c = [removed_c;repmat(colors(j,:),numel(x_removed),1)]; %#ok<AGROW>

    fprintf('Fz=%7.1f N: kept %d points, removed %d points.\n', ...
        Fz_nom,numel(xk),numel(x_removed));

    x_plot = linspace(min(xk),max(xk),250).';
    plot(rad2deg(xk),yk,'o','Color',colors(j,:),'MarkerSize',4, ...
        'HandleVisibility','off');
    plot(rad2deg(x_plot),pacejka4(p_final,x_plot),'-', ...
        'Color',colors(j,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('Fz=%.0f N local',Fz_nom));
end

if ~isempty(removed_x)
    scatter(rad2deg(removed_x),removed_y,40,removed_c,'x', ...
        'LineWidth',1.5,'HandleVisibility','off');
end
xlabel('Slip angle \alpha [deg]');
ylabel('F_y [N]');
title('Local per-load F_y fits; x marks rejected data');
grid on; legend('show','Location','best');

%% ---------------- STEP 4: BUILD GLOBAL INITIAL GUESS ------------------
if isnan(Fz0_user)
    Fz0 = median(Fz_loc);
else
    Fz0 = Fz0_user;
end
if ~isfinite(Fz0) || Fz0 <= 0
    error('Fz0 must be positive and finite.');
end

dfz_loc = (Fz_loc-Fz0)./Fz0;

% Normalize the dimensional local amplitude and stiffness before fitting
% their load-law coefficients.
d_norm_loc = D_loc./Fz_loc;
k_norm_loc = K_loc./Fz_loc;

pd0 = fit_poly_ascending(dfz_loc,d_norm_loc,D_order);
pc0 = fit_poly_ascending(dfz_loc,C_loc,C_order);
pe0 = fit_poly_ascending(dfz_loc,E_loc,E_order);
pk0 = fit_poly_ascending(dfz_loc,k_norm_loc,K_order);

nD = numel(pd0);
nC = numel(pc0);
nE = numel(pe0);
nK = numel(pk0);
p0_global = [pd0,pc0,pe0,pk0];

% Coefficient bounds mainly establish a canonical, numerically stable
% solution. Pointwise load-law behavior is checked after fitting.
lb = -Inf(size(p0_global));
ub = Inf(size(p0_global));

iD0 = 1;
iC0 = nD+1;
iE0 = nD+nC+1;

% D/Fz must be positive at the reference load; C and E use common
% four-parameter Magic Formula ranges.
lb(iD0) = max(1e-8,0.01*max(abs(pd0(1)),1e-6));
ub(iD0) = max(10*max(abs(pd0(1)),1e-3),lb(iD0)*10);
lb(iC0) = 0.50;
ub(iC0) = 2.50;
lb(iE0) = -5.0;
ub(iE0) = 2.0;

p0_global = min(max(p0_global,lb),ub);
X_global = [x_clean_all,Fz_clean_all];

global_model = @(p,X) eval_fy_global(p,X,nD,nC,nE,nK,Fz0);
fit_opts_global = optimoptions('lsqcurvefit', ...
    'Display','iter', ...
    'MaxFunctionEvaluations',30000, ...
    'MaxIterations',5000);

p_global = lsqcurvefit(global_model,p0_global,X_global,y_clean_all, ...
    lb,ub,fit_opts_global);

id = 0;
pd = p_global(id+(1:nD)); id = id+nD;
pc = p_global(id+(1:nC)); id = id+nC;
pe = p_global(id+(1:nE)); id = id+nE;
pk = p_global(id+(1:nK));

fprintf('\nGlobal load-sensitive Fy fit complete.\n');
fprintf('Reference load Fz0 = %.6g N\n',Fz0);
fprintf('D/Fz coefficients [ascending dfz powers]:\n'); disp(pd);
fprintf('C coefficients [ascending dfz powers]:\n'); disp(pc);
fprintf('E coefficients [ascending dfz powers]:\n'); disp(pe);
fprintf('Ky/Fz coefficients [ascending dfz powers]:\n'); disp(pk);

%% ---------------- STEP 5: PARAMETER TRENDS ----------------------------
Fz_plot = linspace(0.9*min(Fz_loc),1.1*max(Fz_loc),300).';
dfz_plot = (Fz_plot-Fz0)./Fz0;

D_plot = Fz_plot.*polyval_ascending(pd,dfz_plot);
C_plot = polyval_ascending(pc,dfz_plot);
E_plot = polyval_ascending(pe,dfz_plot);
K_plot = Fz_plot.*polyval_ascending(pk,dfz_plot);
B_plot = safe_divide(K_plot,C_plot.*D_plot);

figure('Name','Fy parameter trends versus Fz');
subplot(2,2,1);
plot(Fz_loc,B_loc,'o','DisplayName','local'); hold on;
plot(Fz_plot,B_plot,'-','LineWidth',1.5,'DisplayName','global');
xlabel('F_z [N]'); ylabel('B'); title('B(F_z)=K_y/(CD)'); grid on; legend show;

subplot(2,2,2);
plot(Fz_loc,C_loc,'o','DisplayName','local'); hold on;
plot(Fz_plot,C_plot,'-','LineWidth',1.5,'DisplayName','global');
xlabel('F_z [N]'); ylabel('C'); title('Shape factor C(F_z)'); grid on; legend show;

subplot(2,2,3);
plot(Fz_loc,D_loc,'o','DisplayName','local'); hold on;
plot(Fz_plot,D_plot,'-','LineWidth',1.5,'DisplayName','global');
xlabel('F_z [N]'); ylabel('D [N]'); title('Amplitude D(F_z)'); grid on; legend show;

subplot(2,2,4);
plot(Fz_loc,E_loc,'o','DisplayName','local'); hold on;
plot(Fz_plot,E_plot,'-','LineWidth',1.5,'DisplayName','global');
xlabel('F_z [N]'); ylabel('E'); title('Curvature E(F_z)'); grid on; legend show;

figure('Name','Fy cornering-stiffness trend');
plot(Fz_loc,K_loc,'o','DisplayName','local BCD'); hold on;
plot(Fz_plot,K_plot,'-','LineWidth',1.5,'DisplayName','global K_y');
xlabel('F_z [N]'); ylabel('K_y = dF_y/d\alpha|_0 [N/rad]');
title('Zero-slip lateral stiffness'); grid on; legend show;

%% ---------------- STEP 6: VALIDATION ----------------------------------
y_pred_all = global_model(p_global,X_global);
resid_all = y_clean_all-y_pred_all;
rmse = sqrt(mean(resid_all.^2));
mae = mean(abs(resid_all));
ss_res = sum(resid_all.^2);
ss_tot = sum((y_clean_all-mean(y_clean_all)).^2);
r2 = 1-ss_res/max(ss_tot,eps);

fprintf('\nGlobal-model validation on cleaned data:\n');
fprintf('RMSE = %.6g N\n',rmse);
fprintf('MAE  = %.6g N\n',mae);
fprintf('R^2  = %.6f\n',r2);

figure('Name','Global load-sensitive Fy validation'); hold on;
for j = 1:nLevels
    Fz_nom = Fz_loc(j);
    idx = abs(Fz_clean_all-Fz_nom) < max(1e-9,1e-10*Fz_nom);
    xj = x_clean_all(idx);
    yj = y_clean_all(idx);
    x_plot = linspace(min(xj),max(xj),250).';
    X_plot = [x_plot,Fz_nom*ones(size(x_plot))];

    plot(rad2deg(xj),yj,'o','Color',colors(j,:),'MarkerSize',4, ...
        'HandleVisibility','off');
    plot(rad2deg(x_plot),global_model(p_global,X_plot),'-', ...
        'Color',colors(j,:),'LineWidth',1.8, ...
        'DisplayName',sprintf('Fz=%.0f N global',Fz_nom));
end
xlabel('Slip angle \alpha [deg]'); ylabel('F_y [N]');
title(sprintf('Global load-sensitive F_y model: RMSE=%.2f N, R^2=%.4f', ...
    rmse,r2));
grid on; legend('show','Location','best');

figure('Name','Fy residual diagnostics');
subplot(1,2,1);
scatter(y_pred_all,resid_all,10,Fz_clean_all,'filled');
yline(0,'k--'); xlabel('Predicted F_y [N]'); ylabel('Residual [N]');
title('Residual versus prediction'); grid on; colorbar;

subplot(1,2,2);
histogram(resid_all,40);
xlabel('Residual [N]'); ylabel('Count'); title('Residual distribution'); grid on;

%% ---------------- STEP 7: SAVE MODEL ----------------------------------
model.type = 'load-sensitive four-parameter Pacejka Fy';
model.p = p_global;
model.D.param = pd;
model.C.param = pc;
model.E.param = pe;
model.K.param = pk;
model.nD = nD;
model.nC = nC;
model.nE = nE;
model.nK = nK;
model.D_order = D_order;
model.C_order = C_order;
model.E_order = E_order;
model.K_order = K_order;
model.Fz0 = Fz0;
model.Fz_range = [min(Fz_loc),max(Fz_loc)];
model.alpha_range = [min(x_clean_all),max(x_clean_all)];
model.coefficient_order = 'ascending powers of dfz';
model.fit_rmse = rmse;
model.fit_mae = mae;
model.fit_r2 = r2;
model.reference_muy = 1;

save(output_mat_file,'model');
fprintf('Saved fitted model to %s\n',output_mat_file);

% Example evaluation after loading the saved model:
%
%   load('pacejka_fy_4param_model.mat','model');
%   X = [alpha(:),Fz(:)];
%   muy = 0.80;
%   Fy = eval_load_sensitive_fy(model.p,X,model.nD,model.nC, ...
%       model.nE,model.nK,model.Fz0,muy);

%% ---------------- LOCAL FUNCTIONS -------------------------------------
function y = eval_fy_global(p,X,nD,nC,nE,nK,Fz0)
alpha = X(:,1);
Fz = X(:,2);
dfz = (Fz-Fz0)./Fz0;

id = 0;
pd = p(id+(1:nD)); id = id+nD;
pc = p(id+(1:nC)); id = id+nC;
pe = p(id+(1:nE)); id = id+nE;
pk = p(id+(1:nK));

D = Fz.*polyval_ascending(pd,dfz);
C = polyval_ascending(pc,dfz);
E = polyval_ascending(pe,dfz);
K = Fz.*polyval_ascending(pk,dfz);
B = safe_divide(K,C.*D);

z = B.*alpha;
q = z-E.*(z-atan(z));
y = D.*sin(C.*atan(q));
end

function p_ascending = fit_poly_ascending(x,y,order)
p_descending = polyfit(x(:),y(:),order);
p_ascending = fliplr(p_descending);
end

function y = polyval_ascending(p,x)
y = zeros(size(x));
for i = numel(p):-1:1
    y = y.*x+p(i);
end
end

function q = safe_divide(numerator,denominator)
tiny = 1e-12*max(1,max(abs(denominator(:))));
bad = abs(denominator) < tiny;
denominator(bad) = tiny.*sign(denominator(bad)+eps);
q = numerator./denominator;
end
