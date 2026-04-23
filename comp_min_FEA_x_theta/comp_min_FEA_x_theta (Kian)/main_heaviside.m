% MMA Optimization for Topology Optimization with Filtering
% By Zahur Ullah 21/5/2025 and edited to add stress-constraint by Kian Das 20/4/2026
% Here it is optimising both density and theta for a fibre-reinforced composite in 2D
clear; clc; close all; warning off
%% Parameters
volfrac = 0.4;                                   % Volume fraction
penal = 3.0;                                     % Penalization factor
rmin = 3.0;                                      % Filter radius
maxiter = 250;
beta     = 2;                                    % Initial sharpness
beta_max = 32;                                   % Maximum sharpness
eta      = 0.5;                                  % Threshold (0.5 = standard, increase to erode, decrease to dilate)
%Material properties composites (from Guowei Ma)
matprop.E1=39e3;                                 % Young's modulus in fiber direction
matprop.E2=8.4e3;                                % Young's modulus perpendicular to fiber direction
matprop.nu12=0.26;                               % Major Poisson's ratio
matprop.nu21=matprop.nu12*matprop.E2/matprop.E1; % Minor Poisson's ratio
matprop.G12=4.2e3;                               % Shear modulus
% Strength allowables for Tsai-Wu
strength.Xt=1062;                                %
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
% Initialize design variables and combine
x = volfrac * ones(numele,1);              % Density variables
theta = (pi/4) * ones(numele,1);           % Fiber direction variables
xval = [x; theta];                         % Combine design variables
% Bounds for densities and fiber directions
xmin_x = 1e-4 * ones(numele,1); % Lower bound for densities
xmax_x = 1 * ones(numele,1); % Upper bound for densities
xmin_theta = -(pi/2) * ones(numele,1); % Lower bound for fiber directions
xmax_theta =  (pi/2) * ones(numele,1); % Upper bound for fiber directions
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
xphy(numele+1:end) = xval(numele+1:end);
%% Optimisation loop
iterationHistory = zeros(maxiter, 5);
change = 1; iter = 0;
while change > 1e-3 && iter < maxiter
    iter = iter + 1;

    % Heaviside projection
    x_tilde = (H*xval(1:numele))./Hs;
    [x_proj,dxphy] = heavisideProjection(x_tilde,beta,eta);
    xphy(1:numele) = x_proj;

    % FE Analysis (Extracting K and KE0 too now for TW)
    [U, K, KE0] = FE_analysis(xphy, penal, numnode, numele, gs, edofMat, coords, conn, freedofs, F, matprop);
    
    % Tsai-Wu constraint (NEW)
    [g_tw, dgtw_dx_raw, dgtw_dtheta,TW,sigmax, TW_gp, vonMises] = TsaiWu(U, K, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs);

    % Objective function and sensitivities
    [c, dc_dx_raw, dc_theta] = objective_function(U, xphy, penal, numele, gs, edofMat, coords, conn, matprop);
    
    % Volume constraint and sensitivities
    [v, dv_dx_raw, dv_theta] = volume_constraint(xphy, volfrac, numele, ve);
    
    % filtering of sensitivites
    % dc_dx = H*(dc_dx./Hs);
    % dv_dx = H*(dv_dx./Hs);
    %dc_theta = H*(dc_theta./Hs);
    %dv_theta = H*(dv_theta./Hs);
    % Filter gradients (NEW)
    % dgtw_dx     = H*(dgtw_dx./Hs);
    %dgtw_dtheta = H*(dgtw_dtheta./Hs);
    % New sensitivities for Heaviside
    % dc_dx   = H*((dc_dx   .* dxphy)./Hs);
    % dv_dx   = H*((dv_dx   .* dxphy)./Hs);
    % dgtw_dx = H*((dgtw_dx .* dxphy)./Hs);
    
    dc_dx_chain   = dc_dx_raw   .* dxphy;
    dv_dx_chain   = dv_dx_raw   .* dxphy;
    dgtw_dx_chain = dgtw_dx_raw .* dxphy;
    dc_dx   = H * (dc_dx_chain   ./ Hs);
    dv_dx   = H * (dv_dx_chain   ./ Hs);
    dgtw_dx = H * (dgtw_dx_chain ./ Hs);
    % Combine sensitivities
    df0dx = [dc_dx; dc_theta];                     % Objective function sensitivities
    dfdx = [ dv_dx(:).',      dv_theta(:).'  ;
            dgtw_dx(:).',     dgtw_dtheta(:).' ];  % Combined constraint sensitivities 

    % Initial values for MMA
    f0val = c;             % Initial objective function value
    fval = [v; g_tw];      % Initial volume constraint value
    
    % fprintf('m = %d, n = %d\n', m, n);
    % disp(size(fval))
    % disp(size(dfdx))
    % disp(size(a))
    % disp(size(c_MMA))
    % disp(size(d))

    % MMA update
    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low1, upp1] = mmasub(m, n, iter, xval, xmin,...
        xmax, xold1, xold2, f0val, df0dx, fval, dfdx, low, upp, a0, a, c_MMA, d);
    low=low1;
    upp=upp1;
    % Update old values
    xold2 = xold1; 
    xold1 = xval;
    %current values of the design variables
    xval = xmma;
    
    % %filter design variables both density and 
    % xphy=xval; 
    %xphy(1:numele) =     (H*xval(1:numele))./Hs; % old density filter 
    %xphy(numele+1:end) = (H*xval(numele+1:end))./Hs; % theta
    xphy(numele+1:end) = xval(numele+1:end);
    % xval=xphy; 

    % Print results
    change_x = max(abs(xval(1:numele) - xold1(1:numele)));
    change_t = max(abs(xval(numele+1:end) - xold1(numele+1:end))) / pi;
    change = max(change_x, change_t);
    fprintf('It %d: Obj = %f, V = %f, g_tw = %f, Change = %f\n', iter, c, v, g_tw, change);
    iterationHistory(iter, :) = [iter, c, v, change, g_tw];

    % Plot design (x and theta)
    % colormap(jet);
    % patch('Faces',conn','Vertices',coords','FaceVertexCData',xphy(1:numele),...
    %       'FaceColor','flat','EdgeColor',[0,0,0]); axis equal off; colorbar
    % clim([0 1]);  drawnow;
    if mod(iter, 5) == 0 || iter == 1
        figure(1); clf;
        patch('Faces',conn','Vertices',coords','FaceVertexCData',xphy(1:numele),...
              'FaceColor','flat','EdgeColor','none'); 
        axis equal tight off; colormap(flipud(gray)); colorbar;
        hold on;
        ind = find(xphy(1:numele) > volfrac); % Only show fibers where there is material
        theta_curr = xphy(numele+1:end);
        x_plot = [x_cen(ind) - halfL*cos(theta_curr(ind)), ...
                  x_cen(ind) + halfL*cos(theta_curr(ind)), ...
                  nan(length(ind),1)]';
        y_plot = [y_cen(ind) - halfL*sin(theta_curr(ind)), ...
                  y_cen(ind) + halfL*sin(theta_curr(ind)), ...
                  nan(length(ind),1)]';
        line(x_plot(:), y_plot(:), 'Color', [1 0 0], 'LineWidth', 0.5); % Red fibers
        title(sprintf('Iter: %d | Obj: %.2f | Stress: %.2f', iter, c, g_tw));
        drawnow;
    end
    % Beta continuation block
    if mod(iter, 50) == 0 && beta < beta_max
        beta = min(beta*2, beta_max);
        fprintf('   >>> Beta updated to: %d\n',beta)
        % Reset MMA internal state to prevent spikes
        %xold1 = xval;
        %xold2 = xval;
        %low   = max(0, xval - move);
        %upp   = min(1, xval + move);
    end
end
warning on
%%
% Measure of non-discreteness (NEW)
nGrey = sum(x > 0.05 & x < 0.95);
M1 = nGrey / numele;
M2 = sum(4 * x .* (1 - x))/n;
disp(M1) % percentage of elements with density between 0.05<x<0.95
disp(M2) % percentage of average greyness (i.e. design is M2% grey )

% plot angles
% Update design variables
x = xphy(1:numele);
theta = xphy(numele+1:end);

figure(2)
patch('Faces',conn','Vertices',coords','FaceVertexCData',x,...
      'FaceColor','flat','EdgeColor','none'); colorbar; %axis equal off;
axis equal; hold on 
len=0.05; 
u_len=len*cos(theta); 
v_len=len*sin(theta); 

ind = find(x > 0.5);        % only plot in solid regions

for k = 1:length(ind)
    e = ind(k);
    theta_e = theta(e);

    x1 = x_cen(e) - halfL*cos(theta_e);
    x2 = x_cen(e) + halfL*cos(theta_e);
    y1 = y_cen(e) - halfL*sin(theta_e);
    y2 = y_cen(e) + halfL*sin(theta_e);

    line([x1 x2], [y1 y2], ...
         'Color','k', ...
         'LineWidth',1.2);
end

% von Mises and Tsai-Wu stress plots
figure(3);
mask = xphy(1:numele) < 0.3;
field_plot1 = vonMises;
field_plot1(mask) = NaN;
patch('Faces', conn', 'Vertices', coords', ...
      'FaceVertexCData', field_plot1, ...
      'FaceColor','flat','EdgeColor','none');
axis equal off; colorbar;
set(gcf, 'Color', 'white')
title('von Mises stress');

figure(4); clf;
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
title(sprintf('Tsai–Wu Index (iteration %d)', iter));
drawnow;

% plot iteration convergence history
figure(5);
plot(iterationHistory(1:iter, 1), iterationHistory(1:iter, 2), '-o');
xlabel('Iteration');
ylabel('Objective Function (Compliance)');
title('Convergence History');
grid on;



