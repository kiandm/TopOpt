% MMA Optimization for Topology Optimization with Filtering
% By Zahur Ullah 21/5/2025 and edited to add Hashin stress-constraint by Kian Das 18/8/2026
% Here it is optimising both density and theta for a fibre-reinforced composite in 2D
% With Heaviside
tic
clear; clc; 
% close all;
warning off
addpath('Functions\')
%% Parameters
volfrac = 0.50; penal = 3.0; rmin_phys = 5; 
maxiter = 3000; theta_init = pi/2;
beta = 1; beta_max = 32; eta = 0.5;
c_max = 200;
crow = 25;
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
[coords, conn, edofMat, numnode, numele, freedofs, F, H]= problem_setup_Lbrac60(rmin_phys);
Hs = full(sum(H,2));      % full vector (sum of a sparse matrix is sparse)
Ht = H';                  % transpose: used for every ADJOINT (sensitivity) filter below
fscale = 0.02*numele;     % Objective scale (gradient per element)
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
xval = [x; theta];                           % Combine design variables
% Bounds for densities and fiber directions
xmin_x = 1e-4 * ones(numele,1);              % Lower bound for densities
xmax_x = 1 * ones(numele,1);                 % Upper bound for densities
xmin_theta = (0) * ones(numele,1);           % Lower bound for fiber directions
xmax_theta =  (pi) * ones(numele,1);         % Upper bound for fiber directions
%%
% INITIALIZE MMA OPTIMIZER
%Reference from: https://www.top3d.app/tutorials/3d-topology-optimization-using-method-of-moving-asymptotes-top3dmma
m     = 2;                          % The number of general constraints.
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
iterationHistory = zeros(maxiter, 5);
itb = inf;
change = 1; iter = 0; M = 100;
converged = (beta >= beta_max) && (change <= 1e-3) && (M <= 5);
while ~converged && iter < maxiter 
% iterationHistory = zeros(maxiter, 5);
% change = 1; iter = 0;
% while change > 1e-3 && iter < maxiter
    iter = iter + 1;
    % Heaviside projection
    x_tilde = (H*xval(1:numele))./Hs;
    [x_proj,dxphy] = heavisideProjection(x_tilde,beta,eta);
    xphy(1:numele) = x_proj;
    xphy(numele+1:end) = filterTheta(xval(numele+1:end), H, Hs);
    % FE Analysis
    [U, K, KE0, dK] = FE_analysis(xphy, penal, numnode, numele, gs, edofMat, coords, conn, freedofs, F, matprop, dphix_ref, dphiy_ref); % ADDED DPHI
    % Hashin constraint
    [g_hs, dgh_dx_raw, dgh_dtheta, TW, ~, vonMises] = Hashin(U, dK, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs, dphix_ref, dphiy_ref); % ADDED DPHI
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
    % Combine sensitivities
    g_c = c-c_max - 1;                             % Compliance constraints c <= c_max
    df0dx = fscale*[dv_dx; dv_theta];                     % Combined objective function sensitivities (dv_theta = 0)
    dfdx = [ crow*dc_dx(:).'/c_max,      crow*dc_theta(:).'/c_max  ;
            dgh_dx(:).',     dgh_dtheta(:).' ];    % Combined constraint sensitivities 
 %%
    % Initial values for MMA
    % f0val = c;             % Initial objective function value
    % fval = [v; g_hs];      % Initial volume constraint value 
    f0val = fscale*(v+1);          % objective: volume fraction (scaled)
    fval  = [crow*g_c; g_hs];      % constraints: compliance row, then the collapsed Hashin row
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
    fprintf('It %d: Obj = %f, V = %f, g_hs = %f, Change = %f, Change in x = %f, Change in theta = %f\n', iter, c, v+1, g_hs, change, change_x, change_t);
    iterationHistory(iter, :) = [iter, c, v+1, change, g_hs];
    % Plot design (x and theta)
    if mod(iter, 5) == 0 || iter == 0
        figure(9); clf;
        patch('Faces',conn','Vertices',coords','FaceVertexCData',xphy(1:numele),...
              'FaceColor','flat','EdgeColor','none'); 
        axis equal tight off; colormap(flipud(gray)); colorbar;
        hold on;
        ind = find(xphy(1:numele) > 0.2); % Only show fibers where there is material
        theta_curr = xphy(numele+1:end);
        x_plot = [x_cen(ind) - halfL*cos(theta_curr(ind)), ...
                  x_cen(ind) + halfL*cos(theta_curr(ind)), ...
                  nan(length(ind),1)]';
        y_plot = [y_cen(ind) - halfL*sin(theta_curr(ind)), ...
                  y_cen(ind) + halfL*sin(theta_curr(ind)), ...
                  nan(length(ind),1)]';
        line(x_plot(:), y_plot(:), 'Color', [1 0 0], 'LineWidth', 0.5); % Red fibers
        title(sprintf('Iter: %d | V: %.3f | Obj: %.2f | Stress: %.2f', iter, v+1, c, g_hs));
        drawnow;
    end
    % Beta continuation block
    if mod(iter, 25) == 0 && beta < beta_max
        beta = min(beta*1.5, beta_max);
        fprintf('   >>> Beta updated to: %d\n',beta)
    end
    M = 100 * sum(4*xphy(1:numele).*(1-xphy(1:numele))) / numele;
    % converged = (beta >= beta_max) && (change <= 1e-3) && (M <= 5);
        if beta >= beta_max && ~isfinite(itb), itb = iter; end
    win        = 20;
    stationary = (iter - itb >= win) && ...
                 (max(iterationHistory(iter-win+1:iter,3)) - min(iterationHistory(iter-win+1:iter,3))) / iterationHistory(iter,3) < 1e-3;   % column 3 = volume fraction
    feasible   = max([g_c; g_hs]) <= 1e-2;
    converged  = (beta >= beta_max) && (M <= 5) && stationary && feasible;
end
warning on
%%
% Measure of non-discreteness
x = xphy(1:numele);
M = 100 * sum(4 * x .* (1 - x))/numele;
disp(M) % percentage of average greyness (i.e. design is M2% grey )
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
c = colorbar;
c.Label.String = 'Fiber Angle (degrees)';
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
% Hashin failure plot
figure(11); clf;
mask = xphy(1:numele) < 0.3;
field_plot2 = TW;
field_plot2(mask) = NaN;
patch('Faces', conn', ...
      'Vertices', coords', ...
      'FaceVertexCData', field_plot2, ...
      'FaceColor', 'flat', ...
      'EdgeColor', 'none');
set(gcf, 'Color', 'white')
axis equal off;
colorbar;
clim([0 1.2]);   % 1 = failure limit
title(sprintf('Hashin Index (iteration %d)', iter));
drawnow;
% plot iteration convergence history
figure(12); clf;
yyaxis left
plot(iterationHistory(1:iter, 1), iterationHistory(1:iter, 3), '-o');
xlabel('Iteration');
ylabel('Objective Function (Volume Fraction)');
yyaxis right
plot(iterationHistory(1:iter, 1), iterationHistory(1:iter, 5), '-o', 'Color', 'r');
ylabel('Hashin Index, g_{hs}');
ax = gca;
ax.YAxis(2).Color = 'r';
title('Convergence History');
grid on;

% Save converged results for later checking against other failure criteria
% (see checkTsaiWu.m)
saveOptResults(xphy, U, numele, numnode, gs, edofMat, coords, conn, matprop, strength, ...
    penal, dphix_ref, dphiy_ref, freedofs, F, TW, g_hs, vonMises, iter, M);

toc