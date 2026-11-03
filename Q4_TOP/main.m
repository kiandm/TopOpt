% main.m 
tic; clear; clc;
%addpath('Mains/')
addpath('Functions/')
% Parameters-----------------------------------------------
volfrac = 0.5; penal = 3.0; 
rmin_mm = 2.0;   % filter radius [mm], physical
unit2mm = 1.0;   % mesh-file units -> mm (1 = already mm, 1000 = metres, 25.4 = inches)
E0 = 1.0; Emin = 1e-9; maxit = 100; tol = 1e-3; nu = 0.3; 
De = E0/(1-nu^2) * [1, nu,     0; 
                  nu,  1,     0; 
                   0,  0, (1-nu)/2]; % Material matrix
% 1. Import mesh
[nn,nodes,nel,enodes,ndof,edof] = createMesh();
nodes(:,2:3) = nodes(:,2:3) * unit2mm;   % all geometry in mm from here on
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
Xe = nodes(:,2); Xe = Xe(enodes);        % nel x 4
Ye = nodes(:,3); Ye = Ye(enodes);
ecent = [mean(Xe,2), mean(Ye,2)];                                  % [mm]
Ae    = 0.5*abs(sum(Xe.*Ye(:,[2 3 4 1]) - Xe(:,[2 3 4 1]).*Ye, 2)); % shoelace area [mm^2]
Atot  = sum(Ae);  he = sqrt(Ae);
% 5. Build filter structure
[W, WA] = SensitivityFilter(ecent, rmin_mm, Ae);
nNbr = full(sum(W > 0, 2)) - 1;
fprintf('Filter: rmin=%.3g mm | h=%.3g..%.3g mm | neighbours min/mean = %d/%.1f\n', ...
        rmin_mm, min(he), max(he), min(nNbr), mean(nNbr));
if min(nNbr) == 0 || rmin_mm < 1.5*max(he)   % rule of thumb
    warning('rmin is small relative to the element size: filter may be ineffective.');
end% 6. Initialise design
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
    dcf = applySensitivityFilter(dc, x, W, WA, Ae);
    % ----- OC update ------------------------
    x = OC_update(x, dcf, Ae, volfrac);
    % x = fmincon()
    % ----- Convergence, logging -------------
    change = max(abs(x - xold));
    fprintf('It. %4d | Obj (c): %12.5e | Vol: %6.4f | chg: %8.4e\n', ...
            it, c, sum(Ae.*x)/Atot, change);
    % ----- Iteration history ----------------
    iterationHistory(it, :) = [it, c, sum(Ae.*x)/Atot, change];
    % ----- Quick plot -----------------------
    %if mod(it,5)==1 || change<tol
        %plotDensityMesh(nodes, enodes, x);
        %drawnow;
    %end
    if change < tol
        break;
    end
end
% Measure of non-discreteness, M
nGrey = sum(x > 0.05 & x < 0.95);
nel = length(x);
M = nGrey / nel;
disp(M)  
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

