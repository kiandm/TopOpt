% MMA Optimization for Topology Optimization with Filtering
% By Zahur Ullah 21/5/2025 and edited to add Hashin stress-constraint by Kian Das 18/8/2026
% Here it is optimising both density and theta for a fibre-reinforced composite in 2D
% With Heaviside
tic
clear; clc; 
close all;
warning off
%% Parameters
volfrac = 0.4; penal = 3.0; rmin = 3; maxiter = 500; theta_init = pi/2;
beta = 1; beta_max = 32; eta = 0.5;
%Material properties composites (from Guowei Ma)
matprop.E1=39e3;                                 % Young's modulus in fiber direction
matprop.E2=8.4e3;                                % Young's modulus perpendicular to fiber direction
matprop.nu12=0.26;                               % Major Poisson's ratio
matprop.nu21=matprop.nu12*matprop.E2/matprop.E1; % Minor Poisson's ratio
matprop.G12=4.2e3;                               % Shear modulus
% Strength allowables for Hashin
strength.Xt=1062;                                
strength.Xc=610; 
strength.Yt=31; 
strength.Yc=118; 
strength.S=72; 
%%
% [coords, conn, edofMat, numnode, numele, freedofs, F, H]= problem_setup_mbb(rmin);
[coords, conn, edofMat, numnode, numele, freedofs, F, H]= problem_setup_Lbrac60(rmin);
Hs = sum(H,2);
U = zeros(2*numnode,1);
gs=gauss_domain(coords,numele,conn,2);
gs1=gauss_domain(coords,numele,conn,1);    %1 gauss point at the centre
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
x = volfrac * ones(numele,1);                % Density variables
theta = theta_init * ones(numele,1);           % Fiber direction variables
xval = [x; theta];                         % Combine design variables
% Bounds for densities and fiber directions
xmin_x = 1e-4 * ones(numele,1); % Lower bound for densities
xmax_x = 1 * ones(numele,1); % Upper bound for densities
xmin_theta = (0) * ones(numele,1); % Lower bound for fiber directions
xmax_theta =  (pi) * ones(numele,1); % Upper bound for fiber directions
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
x_cen=zeros(numele,1); y_cen=x_cen; 
for i=1:numele
    x_cen(i)=mean(coords(1, conn(:,i)));
    y_cen(i)=mean(coords(2, conn(:,i)));
end 
barLength = 1;              % TOTAL bar length for fibre angle plotting
halfL = barLength / 2;      % plot from middle of element
% Heaviside projection
x_tilde = (H*xval(1:numele))./Hs;
[x_proj, ~] = heavisideProjection(x_tilde, beta, eta);
xphy(1:numele) = x_proj;
%xphy(numele+1:end) = xval(numele+1:end);
p1 = cos(xval(numele+1:end)); p2 = sin(xval(numele+1:end));
xphy(numele+1:end) = atan2((H*p2)./Hs, (H*p1)./Hs);
%% Optimisation loop
iterationHistory = zeros(maxiter, 5);
change = 1; iter = 0;
while change > 1e-3 && iter < maxiter
    iter = iter + 1;
    % Heaviside projection
    x_tilde = (H*xval(1:numele))./Hs;
    [x_proj,dxphy] = heavisideProjection(x_tilde,beta,eta);
    xphy(1:numele) = x_proj;
    % FE Analysis
    [U, K, KE0, dK] = FE_analysis(xphy, penal, numnode, numele, gs, edofMat, coords, conn, freedofs, F, matprop);
    % Tsai-Wu constraint
    %[g_tw, dgtw_dx_raw, dgtw_dtheta, TW, ~, vonMises] = TsaiWu(U, K, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs);
    % Hashin constraint
    [g_hs, dgh_dx_raw, dgh_dtheta, TW, ~, vonMises] = Hashin(U, dK, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs);
    % Objective function and sensitivities
    [c, dc_dx_raw, dc_theta] = objective_function(U, xphy, penal, numele, gs, edofMat, coords, conn, matprop);   
    % Volume constraint and sensitivities
    [v, dv_dx_raw, dv_theta] = volume_constraint(xphy, volfrac, numele, ve); 
%%
    % filtering of sensitivites 
    % sensitivities in theta
    p1_tilde = (H * cos(xval(numele+1:end))) ./ Hs;
    p2_tilde = (H * sin(xval(numele+1:end))) ./ Hs;
    R2 = max(p1_tilde.^2 + p2_tilde.^2, 1e-6);
    dtheta_dp1 = -p2_tilde ./ R2;
    dtheta_dp2 =  p1_tilde ./ R2;
    dc_dp1 = dc_theta .* dtheta_dp1;
    dc_dp2 = dc_theta .* dtheta_dp2;
    dc_theta = -sin(xval(numele+1:end)) .* (H * (dc_dp1 ./ Hs)) ...
              + cos(xval(numele+1:end)) .* (H * (dc_dp2 ./ Hs));
    % dv_dtheta is zero anyway since volume doesn't depend on fibre
    % direction
    dgh_dp1 = dgh_dtheta .* dtheta_dp1;
    dgh_dp2 = dgh_dtheta .* dtheta_dp2;
    dgh_dtheta = -sin(xval(numele+1:end)) .* (H * (dgh_dp1 ./ Hs)) ...
                 + cos(xval(numele+1:end)) .* (H * (dgh_dp2 ./ Hs));
    % sensitivites in x
    dc_dx_chain   = dc_dx_raw   .* dxphy; % Chain rule
    dv_dx_chain   = dv_dx_raw   .* dxphy; % Chain rule
    dgh_dx_chain = dgh_dx_raw .* dxphy; % Chain rule
    dc_dx   = H * (dc_dx_chain   ./ Hs);  % Filter
    dv_dx   = H * (dv_dx_chain   ./ Hs);  % Filter
    dgh_dx = H * (dgh_dx_chain ./ Hs);  % Filter
    % Combine sensitivities
    df0dx = [dc_dx; dc_theta];                     % Combined objective function sensitivities
    dfdx = [ dv_dx(:).',      dv_theta(:).'  ;
            dgh_dx(:).',     dgh_dtheta(:).' ];  % Combined constraint sensitivities 
 %%
    % Initial values for MMA
    f0val = c;             % Initial objective function value
    fval = [v; g_hs];      % Initial volume constraint value    
    % MMA update
    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low1, upp1] = mmasub(m, n, iter, xval, xmin,...
        xmax, xold1, xold2, f0val, df0dx, fval, dfdx, low, upp, a0, a, c_MMA, d);
    low=low1;      upp=upp1;
    xold2 = xold1; xold1 = xval; % Update old values
    xval = xmma;                 % current values of the design variables
    % %filter theta with Cartesian components
    p1 = cos(xval(numele+1:end)); p2 = sin(xval(numele+1:end));
    xphy(numele+1:end) = atan2((H*p2)./Hs, (H*p1)./Hs);
    % xphy(numele+1:end) = (H*xval(numele+1:end))./Hs;
    % Print results
    change_x = max(abs(xval(1:numele) - xold1(1:numele)));
    change_t = max(abs(xval(numele+1:end) - xold1(numele+1:end))) / pi;
    change = max(change_x, change_t);
    fprintf('It %d: Obj = %f, V = %f, g_hs = %f, Change = %f, Change in x = %f, Change in theta = %f\n', iter, c, v, g_hs, change, change_x, change_t);
    iterationHistory(iter, :) = [iter, c, v, change, g_hs];
    % Plot design (x and theta)
    if mod(iter, 5) == 0 || iter == 0
        figure(1); clf;
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
        title(sprintf('Iter: %d | Obj: %.2f | Stress: %.2f', iter, c, g_hs));
        drawnow;
    end
    % Beta continuation block
    if mod(iter, 50) == 0 && beta < beta_max
        beta = min(beta*2, beta_max);
        fprintf('   >>> Beta updated to: %d\n',beta)
    end
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
figure(2); clf;
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
figure(3); clf;
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
figure(4); clf;
yyaxis left
plot(iterationHistory(1:iter, 1), iterationHistory(1:iter, 2), '-o');
xlabel('Iteration');
ylabel('Objective Function (Compliance)');
yyaxis right
plot(iterationHistory(1:iter, 1), iterationHistory(1:iter, 5), '-o', 'Color', 'r');
ylabel('Hashin Index, g_{hs}');
ax = gca;
ax.YAxis(2).Color = 'r';
title('Convergence History');
grid on;

toc