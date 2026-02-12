function main_FMINCON
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

x0 = volfrac * ones(nel,1);

P = 1.0;  % load magnitude

% Objective
function [cf, gradf] = obj_fun(xvar)
    [Kf, Ff] = assembleSystem(nodes, enodes, ndof, edof, KE0, xvar, penal, E0, Emin);
    Ff(2*loads_xy(:)) = Ff(2*loads_xy(:)) - P;
    Uf = solveSystem(ndof, Kf, Ff, BCs_xy);
    [cf, dcf_local] = complianceandsensitivities(Uf, edof, KE0, xvar, penal, E0, Emin);
    gradf = applySensitivityFilter(dcf_local, xvar, Wi, Wj, Ww, Wsum);
end

% Constraints – prefer linear inequality for volume:
A = (1/numel(x0)) * ones(1,numel(x0));  b = volfrac;
Aeq = []; beq = [];
lb = zeros(size(x0));  ub = ones(size(x0));

% Options
opts = optimoptions('fmincon', 'Algorithm','interior-point', ...
    'MaxIter',100,...
    'TolX',1e-12, ...
    'TolFun', 1e-12, ... 
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'Display','iter');

[x, fval, exitflag, output] = fmincon(@obj_fun, x0, A, b, Aeq, beq, lb, ub, [], opts);

% Post-process
figure; plotDensityMesh(nodes, enodes, x); title('Final Design Density'); colorbar;
% plot iteration convergence history
figure;
plot(iterationHistory(1:it, 1), iterationHistory(1:it, 2), '-o');
xlabel('Iteration');
ylabel('Objective Function (Compliance)');
title('Convergence History');
grid on;
toc
end
