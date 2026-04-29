% MMA Optimization for Topology Optimization for composites
% By Zahur Ullah 23/5/2025
% Compliance minimisation with volume and stress constraints 
% For stress we use Tsai-Wu failure criterion
% It is ptimising both density and theta
tic
clear; clc; close all; warning off

% Parameters
volfrac = 0.4; % Volume fraction
pl = 3; % Penalization factor
rmin = 1.5; % Filter radius
p = 8; % P-norm for the global strength constraint
s_pl=0.8; % penality for Stress 
max_itr=500; 


%Material properties composites (from Guowei Ma)
matprop.E1=39e3; % Young's modulus in fiber direction
matprop.E2=8.4e3; % Young's modulus perpendicular to fiber direction
matprop.nu12=0.26; % Major Poisson's ratio
matprop.nu21=matprop.nu12*matprop.E2/matprop.E1; % Minor Poisson's ratio
matprop.G12=4.2e3; % Shear modulus


% %Material properties composites (from Guowei Ma)
% matprop.E1=1; % Young's modulus in fiber direction
% matprop.E2=1; % Young's modulus perpendicular to fiber direction
% matprop.nu12=0.3; % Major Poisson's ratio
% matprop.nu21=matprop.nu12*matprop.E2/matprop.E1; % Minor Poisson's ratio
% matprop.G12=matprop.E1/(2*(1+matprop.nu12)); % Shear modulus
Xt=1062; Xc=610; Yt=31; Yc=118; S=72; 

F1=1/Xt - 1/Xc; 
F2=1/Yt - 1/Yc; 
F11=1/(Xt*Xc); 
F22=1/(Yt*Yc); 
F12=-(1/2)*sqrt(F11*F22); 
F66=1/S^2; 

matprop.F1=F1; 
matprop.F2=F2; 
matprop.F11=F11; 
matprop.F22=F22; 
matprop.F12=F12; 
matprop.F66=F66; 

%[coords, conn, edofMat, numnode, numele, freedofs, F, H]= problem_setup_hook(rmin);
[coords, conn, edofMat, numnode, numele, freedofs, F, H]= problem_setup_Lbrac60(rmin);
%[coords, conn, edofMat, numnode, numele, freedofs, F, H]= problem_setup_mbb(rmin);

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
% x = 1 * ones(numele,1); % Density variables
%theta = 0 * ones(numele,1); % Fiber direction variables
%%
% Principal stress initialisation for theta
fprintf('Calculating initial fiber directions via isotropic solve...\n');
% Temporarily store original properties
matprop_true = matprop;
% Set material properties to isotropic
matprop.E1 = 1000; matprop.E2 = 1000;
matprop.nu12 = 0.3; matprop.nu21 = 0.3;
matprop.G12 = matprop.E1 / (2*(1 + matprop.nu12));
% Set a uniform density and dummy theta for the solve
x_iso = volfrac * ones(numele,1);
theta_iso = zeros(numele,1); 
xval_iso = [x_iso; theta_iso];
% Run ONE FE analysis (penal = 1 for linear solve)
[U_iso, K_iso] = FE_analysis(xval_iso, 1.0, numnode, numele, gs, edofMat, coords, conn, freedofs, F, matprop);
% Calculate stress directions
[~, stressDirection] = calculatePrincipalStress(U_iso, numele, gs, edofMat, coords, conn, matprop);
% Assign result to actual theta and restore properties
theta = stressDirection;
% xval(numele+1:end) = theta;  
matprop = matprop_true; % Restore real composite properties
fprintf('Initial theta assigned. Starting optimisation loop...\n');
%%
% Combine design variables
xval = [x; theta];

% Bounds for densities and fiber directions
xmin_x = 1e-4 * ones(numele,1); % Lower bound for densities
xmax_x = 1    * ones(numele,1); % Upper bound for densities
xmin_theta = -pi * ones(numele,1); % Lower bound for fiber directions
xmax_theta =  pi * ones(numele,1); % Upper bound for fiber directions

% INITIALIZE MMA OPTIMIZER
%Reference from: https://www.top3d.app/tutorials/3d-topology-optimization-using-method-of-moving-asymptotes-top3dmma
m     = 2;     % The number of general constraints.
n     = numel(xval); % The number of design variables x_j.

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

cI_1=1; % Initialize correction factor for previous iteration
alpha_I = 0.5; 
gpn_I_1 = 1; 
ge_I_1 = ones(numele, 1);

% Optimization loop
change = 1; iter = 0;

data_Itr_gmax_c=zeros(max_itr,3); 
data_Itr_x_ge  =zeros(numele,3,iter);

while change > 1e-5 && iter < max_itr 
    iter = iter + 1;

    % FE Analysis
    [K, U] = FE_analysis(xphy, pl, numnode, numele, gs, edofMat, coords, ...
                         conn, freedofs, F, matprop);

    % Objective function and sensitivities
    [c, dc_dx, dc_theta] = obj_func_sens(U, xphy, pl, numele, gs, edofMat, ...
                                         coords, conn, matprop);
    
    % Volume constraint and sensitivities
    [v_con, dv_dx, dv_theta] = volume_const(xphy, volfrac, numele, ve);

    % %increase p in each iteration (to reduce ge at the end<0)

    % if iter<30
    %     p=p0 
    % else 
    %     p0=p0+1; 
    %     p=min(40, p0 + 1)
    % end 

    % Global stess constraint
    [gpn, gmax_con, ge, dgmax_x, dgmax_theta, cI] = stress_const(xphy, coords, conn, ...
        matprop, pl, p, s_pl, gs, gs1, numele, U, edofMat, K, numnode,freedofs,cI_1,alpha_I,gpn_I_1,ge_I_1,H,Hs);
             
    % update values
    cI_1=cI; gpn_I_1=gpn; ge_I_1=ge; 

    % filtering of sensitivites
    dc_dx = H*(dc_dx./Hs);
    dc_theta = H*(dc_theta./Hs);
    dv_dx = H*(dv_dx./Hs);
    dv_theta = H*(dv_theta./Hs);
    dgmax_x = H*(dgmax_x./Hs);
    dgmax_theta = H*(dgmax_theta./Hs);

    % objective function and constraint values for MMA
    f0val = c; % objective function value
    % fval = [v_con]; % volume constraint value only
    fval = [gmax_con; v_con]; % volume and stress constraint value

    % Combine sensitivities
    df0dx = [dc_dx; dc_theta]; % Objective function sensitivitie
    
    % only volume constraint
    % dfdx = [dv_dx'   dv_theta']; 

    %volume and stress constraints
    dfdx = [dgmax_x' dgmax_theta'
            dv_dx'   dv_theta']; % Volume and stress constraint sensitivities (transpose)

    % MMA update
    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low1, upp1] = mmasub(m, numel(xval), iter, xval, xmin,...
        xmax, xold1, xold2, f0val, df0dx, fval, dfdx, low, upp, a0, a, c_MMA, d);
   
    low=low1;
    upp=upp1;
    
    % Update old values
    xold2 = xold1; 
    xold1 = xval;

    %current values of the design variables
    xval = xmma;

    %filter design variables both density and 
    xphy(1:numele) =     (H*xval(1:numele))./Hs; % density 
    xphy(numele+1:end) = (H*xval(numele+1:end))./Hs; % theta

    % Print results
    change = max(abs(xval(1:numele) - xold1(1:numele)));

    %save data to file 
    data_Itr_gmax_c(iter,:)=[iter  max(ge)  c]; 
    data_Itr_x_ge(:,:,iter)=[xphy(1:numele) xphy(numele+1:end) ge];

    fprintf('It %d: Obj = %f, v_con = %f, , stress_con = %f, Change = %f\n', iter, c, v_con, gmax_con, change);
    
    % Plot design
    subplot(2,1,1)
    colormap(jet);
    patch('Faces',conn','Vertices',coords','FaceVertexCData',xphy(1:numele),...
          'FaceColor','flat','EdgeColor',[0,0,0]); axis equal off; colorbar;  
    caxis([0 1]); drawnow;

    subplot(2,1,2)
    colormap(jet);
    patch('Faces',conn','Vertices',coords','FaceVertexCData',ge,...
          'FaceColor','flat','EdgeColor',[0,0,0]); axis equal off; colorbar; 
    caxis([min(ge) max(ge)]); drawnow;
end
warning on

sim_time=toc;
save('L_bracket_div100', 'conn', 'coords', 'data_Itr_gmax_c', 'data_Itr_x_ge','sim_time');  




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

ind=find(xphy(1:numele)>0.3);
quiver(x_cen(ind),y_cen(ind),u_len(ind),v_len(ind),'Color','r','LineWidth',2); 





