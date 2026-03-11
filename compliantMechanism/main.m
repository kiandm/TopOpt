% main_MMA.m 
tic; clear; clc;
%addpath('Mains/')
addpath('Functions/')
dbstop if error
% Parameters-----------------------------------------------
volfrac = 0.25; penal = 3.0; rmin = 3; E0 = 1.0; 
Emin = 1e-4; maxit = 180; tol = 1e-3; nu = 0.3; 
k_in = 1; k_out = 0.001;
beta = 2; eta = 0.5;
De = E0/(1-nu^2) * [1, nu,     0; 
                  nu,  1,     0; 
                   0,  0, (1-nu)/2]; % Material matrix
% 1. Import mesh
[nn,nodes,nel,enodes,ndof,edof] = createMesh();
% 2. Define BC and load nodesets
[BCs_xy, LC,force_in, disp_out] = loadConditions();
disp(nodes(LC, :))
% Show mesh and loads
%plotGeometry(nodes, enodes, BCs_xy, LC, force_in, disp_out);
% 3. Precompute KE0 per element
KE0 = cell(nel,1);
gp = [-1, 1]/sqrt(3); w = [1, 1];
for e = 1:nel
    n  = enodes(e,:);
    xe = nodes(n,2); ye = nodes(n,3);
    KE0{e} = assembleElementStiffnessQ4(xe, ye, De, gp, w);
end
% 4. Precompute element centroids for filtering
ecent = zeros(nel,2);
for e = 1:nel
    n  = enodes(e,:);
    xy = [nodes(n,2), nodes(n,3)];
    ecent(e,:) = mean(xy,1);
end
% 5. Build filter structure
[Wi, Wj, Ww, Wsum] = SensitivityFilter(ecent, rmin);  % sparse weights
% 6. Initialise design
x = volfrac * ones(nel,1);

% 7. Initialise MMA OPTIMIZER
%Reference from: https://www.top3d.app/tutorials/3d-topology-optimization-using-method-of-moving-asymptotes-top3dmma
m     = 1;                          % The number of general constraints.
n     = nel;                        % The number of design variables x_j.
xmin = zeros(n,1);
xmax = ones(n,1);
xval = x;

%xmin = max(0, xval - 0.2);
%xmax = min(1, xval + 0.2);

xold1 = xval;                       % xval, one iteration ago (provided that iter>1).
xold2 = xval;                       % xval, two iterations ago (provided that iter>2).
low   = xmin;  %ones(n,1);          % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = xmax;  %ones(n,1);          % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                          % The constants a_0 in the term a_0*z.
a     = zeros(m,1);                 % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 10000*ones(m,1);            % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);                 % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
% Filter design variable 
%xphy=xval; % 
x = applySensitivityFilter(x, x, Wi, Wj, Ww, Wsum); % density filter
x_phy = heavisideProjection(x, beta, eta);
%(xphy is used only in FE_analysis, objective_function, and final ploting
% but is not not used in the MMA. In the MMA unfiltered x i used)

% 8. Optimization loop
iterationHistory = zeros(maxit, 4);  % Preallocate iteration history array
for it = 1:maxit
    fprintf('--- Iteration %d start ---\n', it);
    xold = x;
    % ----- Assemble K(x) and F -------------
    x = applySensitivityFilter(x, x, Wi, Wj, Ww, Wsum);
    x_phy = heavisideProjection(x, beta, eta);
    fprintf('Assembling system...\n');
    [K, F] = assembleSystem(nodes, enodes, ndof, edof, KE0, x_phy, penal, E0, Emin);
    fprintf('Assembly done.\n');
    % Adding spring BCs
    fprintf('Adding springs...\n');
    for i = 1:length(force_in)
        nid = force_in(i);
        k = 1;
        dof_x = 2*nid-1;
        K(dof_x,dof_x) = K(dof_x,dof_x) + k;
    end
    for i = 1:length(disp_out)
        nid = disp_out(i);
        k = 0.001;
        dof_x = 2*nid-1;
        K(dof_x,dof_x) = K(dof_x,dof_x) + k;
    end
    fprintf('Springs done.\n');
    % ----- Apply loads ---------------------
    F_in = zeros(ndof,1);
    Pin = 1.0;
    F_in(2*force_in(:)-1) = Pin;
    fprintf('Solving U_in...\n');
    U_in = solveSystem(ndof,K,F_in,BCs_xy,LC);
    fprintf('U_in done.\n');


    F_out = zeros(ndof,1);
    F_out(2*disp_out(:)-1) = 1.0;
    fprintf('Solving U_out...\n');
    U_out = solveSystem(ndof, K, F_out, BCs_xy, LC);
    fprintf('U_out done.\n');
    % ----- Solve ----------------------------
    %U = solveSystem(ndof, K, F, BCs_xy, LC);

    fprintf("Solving U...\n");
    U = solveSystem(ndof, K, F, BCs_xy, LC);
    fprintf("Solved U.\n");
    %fprintf('rank(K) = %d (ndof = %d)\n', rank(full(K)), ndof);

    % ----- Compliance and sensitivities -----
    fprintf('Computing sensitivities...\n');
    [c, dc] = complianceandsensitivities(U_in, U_out, F_out, edof, KE0, x_phy, penal, E0, Emin);
    fprintf('Sensitivities done. c = %e\n', c);
    % ----- Sensitivity filter ---------------
    fprintf('Applying filter...\n');
    %dcf = applySensitivityFilter(dc, x, Wi, Wj, Ww, Wsum);
    [~, dxphy] = heavisideProjection(x, beta, eta);
    dc_projected = dc .* dxphy;
    dcf = applySensitivityFilter(dc_projected, x, Wi, Wj, Ww, Wsum);
    fprintf('Filter done.\n');
    % ----- MMA routine ----------------------
    %x = OC_update(x, dcf, volfrac);
    
    % Initial values for MMA
    f0val = c; % Initial objective function value
    df0dx = dcf;
    fval = mean(x) - volfrac; % Initial volume constraint value
    dfdx = (1/nel) * ones(1,n);
    % MMA call
    fprintf('Running MMA...\n');
    move = 0.2;
    xmin = max(0, xval - move);
    xmax = min(1, xval + move);
    [xnew,~,~,~,~,~,~,~,~,low,upp] = mmasub( ...
              m, n, it, xval, xmin, xmax, xold1, xold2, ...
                f0val, df0dx, fval, dfdx, low, upp, a0, a, c_MMA, d );
    xold2 = xold1;
    xold1 = xval;
    xval  = xnew;
    x = xnew;     % updated design
    fprintf('MMA done.\n');
    if mod(it, 30) == 0
        beta = min(beta * 2, 64);
        fprintf('Beta updated to %d\n', beta);
    end
    fprintf('--- Iteration %d complete ---\n', it);
    % ----- Convergence, logging -------------
    change = max(abs(x - xold));
    fprintf('It. %4d | Obj (c): %12.5e | Vol: %6.4f | chg: %8.4e\n', ...
            it, c, mean(x), change);
    
    % ----- Iteration history ----------------
    iterationHistory(it, :) = [it, c, mean(x), change];
    % ----- Quick plot -----------------------
    fprintf('Running plot...\n');
    if mod(it,5)==1 || change<tol
        plotDensityMesh(nodes, enodes, x);
        drawnow;
    end
    fprintf('plot done.\n');
    if change < tol
        break;
    end
end

% Measure of non-discreteness, M
nGrey = sum(x > 0.05 & x < 0.95);
M = nGrey / n;
disp(M)    

% plot iteration convergence history
figure;
plot(iterationHistory(1:it, 1), iterationHistory(1:it, 2), '-o');
xlabel('Iteration');
ylabel('Objective Function');
title('Convergence History');
grid on;
% Plot final design density
figure;
plotDensityMesh(nodes, enodes, x);
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Final Design Density');
colorbar;
grid on;

toc

