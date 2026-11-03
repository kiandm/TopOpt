% MMA Optimization for Topology Optimization with Filtering
% By Zahur Ullah 21/5/2025 and edited to add Hashin stress-constraint by Kian Das 18/8/2026
% Here it is optimising both density and theta for a fibre-reinforced composite in 2D
% With Heaviside
clearvars -except opts; clc;
tStart = tic;
addpath('Functions')                 % forward slash: works on Windows, Linux and Octave
%% Options: every run option lives in Functions/setOpts.m (defaults + validation).
% Run with the defaults:   main_heaviside
% Run with other options:  opts = struct('agg','smoothmax','p',8);  main_heaviside
if ~exist('opts','var'), opts = struct(); end
opts = setOpts(opts);
penal = opts.penal;  rmin_phys = opts.rmin_phys;  maxiter = opts.maxiter;  theta_init = opts.theta_init;
if strcmp(opts.projection, 'heaviside')
    beta = opts.beta0;  beta_max = opts.beta_max;      % projection on: beta continuation from beta0 to beta_max
else
    beta = 0;           beta_max = 0;                  % projection off: heavisideProjection(x,0,eta) returns x, so xphy = filtered density
end
eta = opts.eta;
fprintf('Options: agg=%s  p=%g  rmin_phys=%g  c_max=%g  projection=%s  beta_max=%g  maxiter=%d   (opts stays in the workspace: "clear opts" returns to the defaults)\n', ...
        opts.agg, opts.p, opts.rmin_phys, opts.c_max, opts.projection, beta_max, opts.maxiter);
c_max = opts.c_max;                  % compliance limit F'*U
crow  = opts.crow;                   % compliance-row scale for MMA ([] = set automatically at iteration 1)
if ~isempty(opts.seed), rng(opts.seed); end
%Material properties composites (from Guowei Ma)
matprop.E1=39e3;                                 % Young's modulus in fiber direction
matprop.E2=8.4e3;                                % Young's modulus perpendicular to fiber direction
matprop.nu12=0.26;                               % Major Poisson's ratio
matprop.nu21=matprop.nu12*matprop.E2/matprop.E1; % Minor Poisson's ratio
matprop.G12=4.2e3;                               % Shear modulus
% Strength allowables for Hashin
strength.Xt=1062;                                % Fibre direction tension (MPa)
strength.Xc=610;                                 % Fibre direction compression (MPa)
strength.Yt=31;                                  % Transverse direction tension (MPa)
strength.Yc=118;                                 % Transverse direction compression (MPa)
strength.S=72;                                   % Shear term (MPa)
%%
if strcmp(opts.problem, 'Lbrac60'), probArgs = {rmin_phys, opts.ndiv}; else, probArgs = {rmin_phys}; end
[coords, conn, edofMat, numnode, numele, freedofs, F, H] = feval(['problem_setup_' opts.problem], probArgs{:});

Hs = full(sum(H,2));      % full vector (sum of a sparse matrix is sparse)
Ht = H';                  % transpose: used for every ADJOINT (sensitivity) filter below
fscale = opts.fscale_coef*numele;     % Objective scale (gradient per element)
U = zeros(2*numnode,1);
gs=gauss_domain(coords,numele,conn,2);
gs1=gauss_domain(coords,numele,conn,1);    %1 gauss point at the centre
% Pre-computing shape function derivatives (ONLY VALID FOR STRUCTURED MESH)
dphix_ref = zeros(4,4); dphiy_ref = zeros(4,4); % (node, gauss-point)
for gp_idx = 1:4
    [~, dphix_ref(:,gp_idx), dphiy_ref(:,gp_idx)] = SF_FE(gs(:,gp_idx), coords, conn);
end
% Calculate the volume of each element
ve=zeros(numele,1);
gcount=0;
for i=1:numele
    for gg=1:2*2 %for 2x2 gauss points
        gcount=gcount+1;
        weight=gs(6,gcount); jac=gs(7,gcount);
        ve(i)=ve(i) + jac*weight;
    end
end
% Initialise design variables and combine
x = 1 * ones(numele,1);                % Density variables
theta = theta_init * ones(numele,1);         % Fiber direction variables
if opts.theta_pert > 0, theta = theta + opts.theta_pert*(rand(numele,1) - 0.5); end
xval = [x; theta];                           % Combine design variables
% Bounds for densities and fiber directions
xmin_x = 1e-4 * ones(numele,1);              % Lower bound for densities
xmax_x = 1 * ones(numele,1);                 % Upper bound for densities
xmin_theta = (0) * ones(numele,1);           % Lower bound for fiber directions
xmax_theta =  (pi) * ones(numele,1);         % Upper bound for fiber directions
%%
% INITIALIZE MMA OPTIMIZER
%Reference from: https://www.top3d.app/tutorials/3d-topology-optimization-using-method-of-moving-asymptotes-top3dmma
% m     = 2;                          % The number of general constraints.
nHS   = 4;  if strcmp(opts.agg, 'smoothmax'), nHS = 1; end   % number of failure constraints
m     = 1 + nHS;                    % general constraints: compliance + failure constraint(s)
n     = numel(xval);                % The number of design variables x_j.
xmin  = [xmin_x; xmin_theta];       % Column vector with the lower bounds for the variables x_j.
xmax  = [xmax_x; xmax_theta];       % Column vector with the upper bounds for the variables x_j.
xold1 = xval;                       % xval, one iteration ago (provided that iter>1).
xold2 = xval;                       % xval, two iterations ago (provided that iter>2).
low   = xmin;  %ones(n,1);          % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = xmax;  %ones(n,1);          % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                          % The constants a_0 in the term a_0*z.
a     = zeros(m,1);                 % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 10000*ones(m,1);            % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);                 % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
xphy=xval;                          % Filter design variable 
% (xphy is used only in FE_analysis, objective_function, and final plotting but is not not used in the MMA. In the MMA unfiltered x i used)
%%
% Precompute element centroids
xrow = coords(1,:); yrow = coords(2,:);
x_cen = mean(xrow(conn), 1)';   y_cen = mean(yrow(conn), 1)';
barLength = 1;              % Total bar length for fibre angle plotting
halfL = barLength / 2;      % plot from middle of element
% Heaviside projection
x_tilde = (H*xval(1:numele))./Hs;
[x_proj, ~] = heavisideProjection(x_tilde, beta, eta);
xphy(1:numele) = x_proj;
% p1 = cos(xval(numele+1:end)); p2 = sin(xval(numele+1:end));
% xphy(numele+1:end) = atan2((H*p2)./Hs, (H*p1)./Hs);
xphy(numele+1:end) = filterTheta(xval(numele+1:end), H, Hs);
%% Optimisation loop
iterationHistory = zeros(maxiter, 10);   % [iter, c, V, change, g(1:4), M, seconds]; g padded with NaN for smoothmax
itb = inf;  exitflag = 'maxiter';
change = 1; iter = 0; M = 100; converged = false;
while ~converged && iter < maxiter
    iter = iter + 1;  tIter = tic;
    % Heaviside projection
    x_tilde = (H*xval(1:numele))./Hs;
    [x_proj,dxphy] = heavisideProjection(x_tilde,beta,eta);
    xphy(1:numele) = x_proj;
    xphy(numele+1:end) = filterTheta(xval(numele+1:end), H, Hs);
    % FE Analysis
    [U, K, KE0, dK] = FE_analysis(xphy, penal, numnode, numele, gs, edofMat, coords, conn, freedofs, F, matprop, dphix_ref, dphiy_ref); % ADDED DPHI
    % Hashin constraint
    [g_hs, dgh_dx_raw, dgh_dtheta, FailIdx, vonMises] = Hashin(U, dK, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs, dphix_ref, dphiy_ref, opts);
    % Objective function and sensitivities
    [c, dc_dx_raw, dc_theta] = objective_function(U, xphy, penal, numele, gs, edofMat, coords, conn, matprop, dphix_ref, dphiy_ref); % ADDED DPHI 
    % Volume constraint and sensitivities
    [v, dv_dx_raw, dv_theta] = volume_constraint(xphy, 1, numele, ve); 
%%
        % filtering of sensitivities: each adjoint filter is applied with Ht = H' (NOT H)
    th      = xval(numele+1:end);
    q1      = (H*cos(2*th))./Hs;   q2 = (H*sin(2*th))./Hs;   R2 = max(q1.^2 + q2.^2, 1e-6);
    dth_dq1 = -0.5*q2./R2;         dth_dq2 = 0.5*q1./R2;     % d(theta_tilde)/d(q1,q2),  theta_tilde = 0.5*atan2(q2,q1)
    sc      = -2*sin(2*th);        cc2     = 2*cos(2*th);    % d(q1,q2)/d(theta)
    % dv_dtheta is zero anyway since volume doesn't depend on fibre direction
    dc_theta   = sc.*(Ht*((dc_theta  .*dth_dq1)./Hs)) + cc2.*(Ht*((dc_theta  .*dth_dq2)./Hs));
    dgh_dtheta = sc.*(Ht*((dgh_dtheta.*dth_dq1)./Hs)) + cc2.*(Ht*((dgh_dtheta.*dth_dq2)./Hs));
    % sensitivities in x: chain rule through the Heaviside projection, then the adjoint of the density filter
    dc_dx  = Ht*((dc_dx_raw .*dxphy)./Hs);
    dv_dx  = Ht*((dv_dx_raw .*dxphy)./Hs);
    dgh_dx = Ht*((dgh_dx_raw.*dxphy)./Hs);
    if isempty(crow)                 % automatic compliance-row scale, fixed at iteration 1 (same for every aggregation)
        crow = norm(fscale*[dv_dx; dv_theta]) / norm([dc_dx; dc_theta]/c_max);
        fprintf('crow = %.3f (automatic)\n', crow);
    end
    % Combine sensitivities
    g_c = c/c_max - 1;                             % Compliance constraints c <= c_max
    df0dx = fscale*[dv_dx; dv_theta];                     % Combined objective function sensitivities (dv_theta = 0)
    dfdx = [ crow*dc_dx(:).'/c_max,      crow*dc_theta(:).'/c_max  ;
            dgh_dx.',     dgh_dtheta.' ];    % Combined constraint sensitivities 
 %%
    % Initial values for MMA
    % f0val = c;             % Initial objective function value
    % fval = [v; g_hs];      % Initial volume constraint value 
    f0val = fscale*(v+1);          % objective: volume fraction (scaled)
    fval  = [crow*g_c; g_hs];      % constraints: compliance row, then the failure constraint row(s)
    % MMA update
    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low1, upp1] = mmasub(m, n, iter, xval, xmin,...
        xmax, xold1, xold2, f0val, df0dx, fval, dfdx, low, upp, a0, a, c_MMA, d);
    low=low1;      upp=upp1;
    xold2 = xold1; xold1 = xval; % Update old values
    xval = xmma;                 % current values of the design variables
    % %filter theta with Cartesian components
    % p1 = cos(xval(numele+1:end)); p2 = sin(xval(numele+1:end));
    % xphy(numele+1:end) = atan2((H*p2)./Hs, (H*p1)./Hs);
    % Print results
    change_x = max(abs(xval(1:numele) - xold1(1:numele)));
    change_t = max(abs(xval(numele+1:end) - xold1(numele+1:end))) / pi;
    change = max(change_x, change_t);
    g4 = nan(4,1);  g4(1:numel(g_hs)) = g_hs;                    % pad to 4 entries: same columns for both aggregations
    M  = 100 * sum(4*xphy(1:numele).*(1-xphy(1:numele))) / numele;   % greyness [%] of the design analysed in this iteration
    if nHS == 4, gstr = sprintf('g_ft = %f, g_fc = %f, g_mt = %f, g_mc = %f', g_hs); else, gstr = sprintf('g_hs = %f', g_hs); end
    fprintf('It %d: C = %f, obj(V) = %f, M = %f, %s, Change = %f, Change in x = %f, Change in theta = %f\n', iter, c, v+1, M, gstr, change, change_x, change_t);
    iterationHistory(iter, :) = [iter, c, v+1, change, g4.', M, toc(tIter)];
    % Plot design (x and theta)
    if opts.plot && mod(iter, 5) == 0
        figure(9); clf;
        patch('Faces',conn','Vertices',coords','FaceVertexCData',xphy(1:numele),...
              'FaceColor','flat','EdgeColor','none'); 
        axis equal tight off; colormap(flipud(gray)); colorbar;
        hold on;
        ind = find(xphy(1:numele) > 0.49); % Only show fibers where there is material
        theta_curr = xphy(numele+1:end);
        x_plot = [x_cen(ind) - halfL*cos(theta_curr(ind)), ...
                  x_cen(ind) + halfL*cos(theta_curr(ind)), ...
                  nan(length(ind),1)]';
        y_plot = [y_cen(ind) - halfL*sin(theta_curr(ind)), ...
                  y_cen(ind) + halfL*sin(theta_curr(ind)), ...
                  nan(length(ind),1)]';
        line(x_plot(:), y_plot(:), 'Color', [1 0 0], 'LineWidth', 0.5); % Red fibers
        title(sprintf('%s \n Iter: %d | Obj(V): %.2f | c: %.1f | max(g): %.2f', opts.agg, iter, v+1, c, max(g_hs)));
        drawnow;
    end
    % Beta continuation block
    if mod(iter, opts.beta_step) == 0 && beta < beta_max
        beta = min(beta*opts.beta_factor, beta_max);
        fprintf('   >>> Beta updated to: %d\n',beta)
    end
    % Stopping rule: volume plateau after beta has reached beta_max, with the constraints satisfied and, if opts.tol_M is
    % finite, the greyness M of the design analysed in this iteration at or below it. M is always logged (history column 9, results.M).
    if beta >= beta_max && ~isfinite(itb), itb = max(iter, opts.min_iter); end   % without projection beta >= beta_max from iteration 1
    win        = opts.win;
    stationary = (iter - itb >= win) && ...
                 (max(iterationHistory(iter-win+1:iter,3)) - min(iterationHistory(iter-win+1:iter,3))) / iterationHistory(iter,3) < opts.tol_plateau;
    feasible   = max([g_c; g_hs]) <= opts.tol_g;
    greyOK     = M <= opts.tol_M;
    converged  = (beta >= beta_max) && stationary && feasible && greyOK;
end
if converged, exitflag = 'converged'; end
%%
% Measure of non-discreteness
%x = xphy(1:numele);
% M = 100 * sum(4 * x .* (1 - x))/numele;
% disp(M) % percentage of average greyness (i.e. design is M2% grey )
fprintf('Exit: %s after %d iterations | V = %.4f | c/c_max = %.4f | max(g) = %.4f | M = %.2f %%\n', exitflag, iter, v+1, c/c_max, max([g_c; g_hs]), M);
if opts.plot
% Plot orientation
theta_rad = xphy(numele+1:end);
theta_deg = mod(rad2deg(theta_rad), 180); % Extract physical angles and convert from radians to degrees [0, 180]
x_dens = xphy(1:numele);
theta_plot = theta_deg;
theta_plot(x_dens <= 0.5) = NaN; % Hide void elements
figure(10); clf;
patch('Faces', conn', ...
      'Vertices', coords', ...
      'FaceVertexCData', theta_plot, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'none');% Plot elements colored by fiber angle
axis equal tight off;
colormap(hsv);             % 'hsv' or 'jet' work well for periodic angles
cb = colorbar;
cb.Label.String = 'Fiber Angle (degrees)';
clim([0 180]);             % Fixed scale from 0° to 180°
set(gcf, 'Color', 'w');
title('Fiber Orientation Field'); % Format colormap, limits, and colorbar
hold on;
barLength = 1;              % Set appropriate length relative to element size
halfL = barLength / 2;
ind = find(x_dens > 0.5);   % Solid elements index
x_lines = [x_cen(ind) - halfL*cos(theta_rad(ind)), ...
           x_cen(ind) + halfL*cos(theta_rad(ind)), ...
           nan(length(ind),1)]';
y_lines = [y_cen(ind) - halfL*sin(theta_rad(ind)), ... % Fixed: y_cen instead of x_cen
           y_cen(ind) + halfL*sin(theta_rad(ind)), ...
           nan(length(ind),1)]';
line(x_lines(:), y_lines(:), 'Color', [0 0 0 0.5], 'LineWidth', 0.8); % Overlay fiber direction vector lines
% Failure-index plot: one panel per mode (separated) or a single panel (smoothmax)
figure(11); clf;
modeNames = {'Fibre Tension', 'Fibre Compression', 'Matrix Tension', 'Matrix Compression'};
if nHS == 1, modeNames = {'Smooth-max Hashin index'}; end
mask = xphy(1:numele) < 0.3;
for k = 1:nHS
    subplot(2,2,k);
    field_plot2 = FailIdx(:,k);
    field_plot2(mask) = NaN;
    patch('Faces', conn', 'Vertices', coords', 'FaceVertexCData', field_plot2, 'FaceColor', 'flat', 'EdgeColor', 'none');
    axis equal off;
    colorbar;
    clim([0 1.2]);   % 1 = failure limit
    title(sprintf('%s (iter %d)', modeNames{k}, iter));
end
set(gcf, 'Color', 'white');
drawnow;
% plot iteration convergence history
figure(12); clf;
yyaxis left
plot(iterationHistory(1:iter, 1), iterationHistory(1:iter, 3), '-o');
xlabel('Iteration');
ylabel('Objective Function (Volume Fraction)');
yyaxis right
gNames = {'g_{ft}', 'g_{fc}', 'g_{mt}', 'g_{mc}'};  if nHS == 1, gNames = {'g_{hs}'}; end
for k = 1:nHS
    plot(iterationHistory(1:iter, 1), iterationHistory(1:iter, 4+k) ...
        , '-o', 'DisplayName', gNames{k}); hold on;
end
ylabel('Failure constraints, g');
legend('Location', 'best');
ax = gca;
ax.YAxis(2).Color = 'r';
title('Convergence History');
grid on;
end   % opts.plot

% Save the final state with everything needed to repeat the run (see saveOptResults.m and checkTsaiWu.m)
if opts.save
    saveOptResults(xphy, U, numele, numnode, gs, edofMat, coords, conn, matprop, strength, ...
        penal, dphix_ref, dphiy_ref, freedofs, F, FailIdx, g_hs, vonMises, iter, M, ...
        opts, exitflag, iterationHistory(1:iter,:), c, v+1, toc(tStart));
end

fprintf('Elapsed time: %.1f s\n', toc(tStart));