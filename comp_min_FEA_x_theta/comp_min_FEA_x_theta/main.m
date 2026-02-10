% MMA Optimization for Topology Optimization with Filtering
% By Zahur Ullah 21/5/2025
% Here it is optimising both density and theta for compoiste in 2D

clear; clc; close all; warning off

% Parameters
volfrac = 0.3; % Volume fraction
penal = 3.0; % Penalization factor
rmin = 1.5; % Filter radius

% % Material properties isotropic
% matprop.E1=1; % Young's modulus in fiber direction
% matprop.E2=1; % Young's modulus perpendicular to fiber direction
% matprop.nu12=0.3; % Major Poisson's ratio
% matprop.nu21=matprop.nu12*matprop.E2/matprop.E1; % Minor Poisson's ratio
% matprop.G12=matprop.E1/(2*(1+matprop.nu12)); % Shear modulus

% Material properties composites
matprop.E1=10; % Young's modulus in fiber direction
matprop.E2=1; % Young's modulus perpendicular to fiber direction
matprop.nu12=0.3; % Major Poisson's ratio
matprop.nu21=matprop.nu12*matprop.E2/matprop.E1; % Minor Poisson's ratio
matprop.G12=matprop.E1/(2*(1+matprop.nu12)); % Shear modulus

% %Material properties composites (from Guowei Ma)
% matprop.E1=39e3; % Young's modulus in fiber direction
% matprop.E2=8.4e3; % Young's modulus perpendicular to fiber direction
% matprop.nu12=0.26; % Major Poisson's ratio
% matprop.nu21=matprop.nu12*matprop.E2/matprop.E1; % Minor Poisson's ratio
% matprop.G12=4.2e3; % Shear modulus

% [coords, conn, edofMat, numnode, numele, freedofs, F, H]= problem_setup_mbb(rmin);
[coords, conn, edofMat, numnode, numele, freedofs, F, H]= problem_setup_Lbrac60(rmin);
Hs = sum(H,2);

U = zeros(2*numnode,1);
gs=gauss_domain(coords,numele,conn,2);
gs1=gauss_domain(coords,numele,conn,1);  %1 gauss pint at the centre

% volume of each element
ve=zeros(numele,1);
gcount=0;
for i=1:numele
    for gg=1:2*2 %for 2x2 gauss points
        gcount=gcount+1;
        weight=gs(6,gcount); jac=gs(7,gcount);
        ve(i)=ve(i) + jac*weight;
    end
end


% Initialize design variables
x = volfrac * ones(numele,1); % Density variables
theta = 0 * ones(numele,1); % Fiber direction variables

% Combine design variables
xval = [x; theta];

% Bounds for densities and fiber directions
xmin_x = 1e-4 * ones(numele,1); % Lower bound for densities
xmax_x = 1 * ones(numele,1); % Upper bound for densities
xmin_theta = -pi * ones(numele,1); % Lower bound for fiber directions
xmax_theta =  pi * ones(numele,1); % Upper bound for fiber directions


% INITIALIZE MMA OPTIMIZER
%Reference from: https://www.top3d.app/tutorials/3d-topology-optimization-using-method-of-moving-asymptotes-top3dmma

m     = 1;                          % The number of general constraints.
n     = numel(xval);                % The number of design variables x_j.

xmin  = [xmin_x; xmin_theta];       % Column vector with the lower bounds for the variables x_j.
xmax  = [xmax_x; xmax_theta];       % Column vector with the upper bounds for the variables x_j.

xold1 = xval;        % xval, one iteration ago (provided that iter>1).
xold2 = xval;        % xval, two iterations ago (provided that iter>2).

low   = xmin;  %ones(n,1);    % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = xmax;  %ones(n,1);    % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).

a0    = 1;                % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 10000*ones(m,1);  % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.

% Filter design variable 
xphy=xval; % 

%(xphy is used only in FE_analysis, objective_function, and final ploting
% but is not not used in the MMA. In the MMA unfiltered x i used)

% Optimization loop
change = 1; iter = 0;
while change > 1e-3 && iter < 50
    iter = iter + 1;

    % FE Analysis
    [U] = FE_analysis(xphy, penal, numnode, numele, gs, edofMat, coords, conn, freedofs, F, matprop);

    % Objective function and sensitivities
    [c, dc_dx, dc_theta] = objective_function(U, xphy, penal, numele, gs, edofMat, coords, conn, matprop);
    
    % Volume constraint and sensitivities
    [v, dv_dx, dv_theta] = volume_constraint(xphy, volfrac, numele, ve);
    
    % filtering of sensitivites
    dc_dx = H*(dc_dx./Hs);
    dv_dx = H*(dv_dx./Hs);

    dc_theta = H*(dc_theta./Hs);
    dv_theta = H*(dv_theta./Hs);

    % Combine sensitivities
    df0dx = [dc_dx; dc_theta]; % Objective function sensitivitie
    dfdx = [dv_dx;  dv_theta]'; % Volume constraint sensitivities (transpose)
    
    % Initial values for MMA
    f0val = c; % Initial objective function value
    fval = v; % Initial volume constraint value
    
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
    xphy(1:numele) =     (H*xval(1:numele))./Hs; % density 
    xphy(numele+1:end) = (H*xval(numele+1:end))./Hs; % theta
    % xval=xphy; 

    % Print results
    change = max(abs(xval(1:numele) - xold1(1:numele)));
    fprintf('It %d: Obj = %f, V const = %f, Change = %f\n', iter, c, v, change);
    
    % Plot design
    colormap(jet);
    patch('Faces',conn','Vertices',coords','FaceVertexCData',xphy(1:numele),...
          'FaceColor','flat','EdgeColor',[0,0,0]); axis equal off; colorbar
    caxis([0 1]);  drawnow;
end
warning on

%plot angles
% Update design variables
x = xphy(1:numele);
theta = xphy(numele+1:end);

figure
patch('Faces',conn','Vertices',coords','FaceVertexCData',x,...
      'FaceColor','flat','EdgeColor',[0 0 0]); colorbar; %axis equal off;
axis equal; hold on 
len=0.05; 
u_len=len*cos(theta); 
v_len=len*sin(theta); 

x_cen=zeros(numele,1); y_cen=x_cen; 
for i=1:numele
    x_cen(i)=mean(coords(1, conn(:,i)));
    y_cen(i)=mean(coords(2, conn(:,i)));
end 

ind=find(xphy(1:numele)>0.5);
quiver(x_cen(ind),y_cen(ind),u_len(ind),v_len(ind),'Color','r','LineWidth',1); 


