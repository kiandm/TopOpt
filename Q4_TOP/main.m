% main.m 
tic; clear; clc;
addpath('Mains/')
addpath('Functions/')
% Parameters-----------------------------------------------
volfrac = 0.5; penal = 3.0; rmin = 2.0; E0 = 1.0; 
Emin = 1e-9; maxit = 200; tol = 1e-3; nu = 0.3; 
De = E0/(1-nu^2) * [1, nu,     0; 
                  nu,  1,     0; 
                   0,  0, (1-nu)/2]; % Material matrix
% 1. Import mesh
[nn,nodes,nel,enodes,ndof,edof] = createMesh();
% 2. Define BC and load nodesets
[BCs_xy,loads_xy] = loadConditions();
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
% 7. Optimisation loop
iterationHistory = zeros(maxit, 4);  % Preallocate iteration history array (saved 3 seconds)
for it = 1:maxit
    xold = x;
    % ----- Assemble K(x) and F --------------
    [K, F] = assembleSystem(nodes, enodes, ndof, edof, KE0, x, penal, E0, Emin);
    % ----- Apply loads ----------------------
    P = 1.0;
    F(2*loads_xy(:)) = F(2*loads_xy(:)) - P;
    % ----- Solve ----------------------------
    U = solveSystem(ndof, K, F, BCs_xy);
    % ----- Compliance and sensitivities -----
    [c, dc] = complianceandsensitivities(U, edof, KE0, x, penal, E0, Emin);
    % ----- Sensitivity filter ---------------
    dcf = applySensitivityFilter(dc, x, Wi, Wj, Ww, Wsum);
    % ----- OC update ------------------------
    x = OC_update(x, dcf, volfrac);
    % x = fmincon()
    % ----- Convergence, logging -------------
    change = max(abs(x - xold));
    fprintf('It. %4d | Obj (c): %12.5e | Vol: %6.4f | chg: %8.4e\n', ...
            it, c, mean(x), change);
    % ----- Iteration history ----------------
    iterationHistory(it, :) = [it, c, mean(x), change];
    % ----- Quick plot -----------------------
    %if mod(it,5)==1 || change<tol
        %plotDensityMesh(nodes, enodes, x);
        %drawnow;
    %end
    if change < tol
        break;
    end
end
% 8. Plot iteration convergence history
figure;
plot(iterationHistory(1:it, 1), iterationHistory(1:it, 2), '-o');
xlabel('Iteration');
ylabel('Objective Function (Compliance)');
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

