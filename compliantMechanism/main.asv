% main.m 
tic; clear; clc;
addpath('Functions/')
dbstop if error
% Parameters-----------------------------------------------
volfrac = 0.25; % Target volume fraction
penal = 3.0;   % SIMP penalisation factor
rmin = 3;      % Radius for filter
E0 = 1.0;      % Initial Young's modulus value
Emin = 1e-4;   % Minimum stiffness (Prevents 0 stiffness elements)

maxit = 300;   % Iteration limit
tol = 1e-3;    % Stopping criteria for change in density
nu = 0.3;      % Poisson's ratio

k_in = 1;      % Input spring constant
k_out = 0.001; % Output spring constant
beta = 0.5;      % Heaviside curvature parameter
eta = 0.5;
De = E0/(1-nu^2) * [1, nu,     0; 
                  nu,  1,     0; 
                   0,  0, (1-nu)/2]; % Material matrix

% Initialise video writer
v = VideoWriter('topology_optimisation_7.mp4', 'MPEG-4');
v.FrameRate = 5;  % frames per second
v.Quality = 95;
open(v);
fig = figure('Visible','off');
%% 
% 1. Import mesh
[nn,nodes,nel,enodes,ndof,edof] = createMesh();

% 2. Define BC and load nodesets
[BCs_xy, LC,force_in, disp_out] = loadConditions();

% Show mesh and loads
%plotGeometry(nodes, enodes, BCs_xy, LC, force_in, disp_out);

% 3. Precompute KE0 per element and element centroids for filtering
KE0 = cell(nel,1);
gp = [-1, 1]/sqrt(3); w = [1, 1];   % Gauss quadrature
ecent = zeros(nel,2);
for e = 1:nel
    n  = enodes(e,:);
    xe = nodes(n,2); ye = nodes(n,3);
    KE0{e} = assembleElementStiffnessQ4(xe, ye, De, gp, w);
    n  = enodes(e,:);
    xy = [nodes(n,2), nodes(n,3)];
    ecent(e,:) = mean(xy,1);
end

% 4. Build filter structure
[Wi, Wj, Ww, Wsum] = SensitivityFilter(ecent, rmin);  % sparse weights

% 5. Initialise design
x = volfrac * ones(nel,1);
c_best = +inf;
%% 
% 6. Initialise MMA OPTIMIZER
%Reference from: https://www.top3d.app/tutorials/3d-topology-optimization-using-method-of-moving-asymptotes-top3dmma
m     = 1;                          % The number of general constraints.
n     = nel;                        % The number of design variables x_j.
xmin = zeros(n,1);
xmax = ones(n,1);
xval = x;
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
%(xphy is used only in FE_analysis, objective_function, and final ploting
% but is not not used in the MMA. In the MMA unfiltered x i used)
%%
% 7. Optimization loop
iterationHistory = zeros(maxit, 4);  % Preallocate iteration history array
for it = 1:maxit
    xold = x;
    
    % Apply filters to get x_tilde then x_phy
    x_tilde = applyDensityFilter(x, Wi, Wj, Ww, Wsum);
    x_phy = heavisideProjection(x_tilde, beta, eta);
    M = sum(4 * x_phy .* (1 - x_phy))/n * 100;

    % Density clipping
    %x_phy = max(x_phy, 1e-3);
    
    % Assemble
    [K] = assembleSystem(nodes, enodes, ndof, edof, KE0, x_phy, penal, E0, Emin);

    % Adding spring BCs
    for i = 1:length(force_in)
        nid = force_in(i);
        k = k_in;
        dof_x = 2*nid-1;
        K(dof_x,dof_x) = K(dof_x,dof_x) + k;
    end

    for i = 1:length(disp_out)
        nid = disp_out(i);
        k = k_out;
        dof_x = 2*nid-1;
        K(dof_x,dof_x) = K(dof_x,dof_x) + k;
    end

    % ----- Apply loads & solve ---------------------
    F_in = zeros(ndof,1);
    F_in(2*force_in(:)-1) = 1.0;

    F_out = zeros(ndof,1);
    F_out(2*disp_out(:)-1) = 1.0;
    
    U_in  = solveSystem(ndof, K, F_in, BCs_xy, LC);
    U_out = solveSystem(ndof, K, F_out, BCs_xy, LC);

    % ----- Compliance and sensitivities -----
    [c, dc] = complianceandsensitivities(U_in, U_out, F_out, edof, KE0, x_phy, penal, E0, Emin);

    % ----- Sensitivity filter ---------------
    [~, dxphy] = heavisideProjection(x_tilde, beta, eta);
    dc_projected = dc .* dxphy;
    dcf = applySensitivityFilter(dc_projected, Wi, Wj, Ww, Wsum);
    
    % ----- MMA routine ----------------------
    %x = OC_update(x, dcf, volfrac);
    
    % Initial values for MMA
    f0val = c;                     % Initial objective function value
    df0dx = dcf;
    fval = mean(x_phy) - volfrac;  % Initial volume constraint value
    %dfdx = (1/nel) * ones(1,n);
    dfdx_raw      = (1/nel) * dxphy;   % dxphy already computed above
    dfdx_filtered = applySensitivityFilter(dfdx_raw, Wi, Wj, Ww, Wsum);
    dfdx          = dfdx_filtered';

    % MMA call
    %move = 0.2;
    move = min(0.2, 0.5/(beta+1));
    xmin = max(0, xval - move);
    xmax = min(1, xval + move);

    [xnew,~,~,~,~,~,~,~,~,low,upp] = mmasub( ...
              m, n, it, xval, xmin, xmax, xold1, xold2, ...
                f0val, df0dx, fval, dfdx, low, upp, a0, a, c_MMA, d );

    xold2 = xold1;
    xold1 = xval;
    xval  = xnew;
    x = xnew;     % updated design
    
    % Beta continuation scheme
    if mod(it, 50) == 0 && beta < 64
        beta = min(beta * 2, 64);
        % Reset MMA memory to avoid oscillation after beta jump
        if M > 20
            xold1 = xval;
            xold2 = xval;
            low   = max(0, xval - move);
            upp   = min(1, xval + move);
        end
        fprintf('Beta updated to %.1f, MMA reset at iteration %d\n', beta, it);
    end

    % ----- Convergence, logging -------------
    change = max(abs(x - xold));
    M = sum(4 * x_phy .* (1 - x_phy))/n * 100;
    fprintf('It. %4d | Obj (c): %12.5e | Vol: %6.4f | chg: %8.4e| M: %6.4f\n', ...
            it, c, mean(x_phy), change, M);

    % Track best design
    if c < c_best
        c_best = c;
        x_best = x;
        x_phy_best = x_phy;
        it_best = it;
    end

    % ----- Iteration history ----------------
    iterationHistory(it, :) = [it, c, mean(x_phy), change];

    % ----- Quick plot -----------------------
    %if mod(it,5)==1 || change<tol
    %    plotDensityMesh(nodes, enodes, x_phy);
    %    drawnow;
    %end

    % Capture frame every iteration
    if mod(it, 5) == 0
        clf(fig);
        plotDensityMesh(nodes, enodes, x_phy);
        title(sprintf('It: %d | c: %.4f | Vol: %.3f | beta: %.1f', it, c, mean(x_phy), beta));
        colorbar;
        %drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end

    if change < tol
        break;
    end

end
close(v);
close(fig);
fprintf('Video saved\n');
%%

% plot iteration convergence history + volume constraint
figure;
yyaxis left
plot(iterationHistory(1:it, 1), iterationHistory(1:it, 2), '-o', 'MarkerSize', 3);
ylabel('Objective Function (c)');
yyaxis right
plot(iterationHistory(1:it, 1), iterationHistory(1:it, 3), '-o', 'MarkerSize', 3);
yline(volfrac, '--', sprintf('Target = %.2f', volfrac), 'LabelHorizontalAlignment', 'left');
ylabel('Volume Fraction');
xlabel('Iteration');
title('Convergence History');
legend('Objective Function', 'Volume Fraction', 'Location', 'best');
grid on;

% Plot best design density
x_phy_best_filtered = heavisideProjection(x_phy_best, beta, eta);
figure;
plotDensityMesh(nodes, enodes, x_phy_best_filtered);
xlabel('X Coordinate');
ylabel('Y Coordinate');
title(sprintf('Filtered Best Design at Iteration %d (c=%.4f)', it_best, c_best));
colorbar;
grid on;

toc

