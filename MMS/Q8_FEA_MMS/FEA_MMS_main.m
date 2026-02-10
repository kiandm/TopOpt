% FEA_MMS_Main.m 
tic
clear; clc;
addpath('Functions/')
% Parameters
E = 1; nu = 0.3; test = 256;  % Mesh refinement level
De = E/(1-nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2]; % Material matrix
% 1. Create mesh
[nodes, enodes, edof, nelx, nely, nn, ndof] = createMesh(test);
% 2. Define MMS functions
mms = mmsFunctions(E, nu);
% 3. Assemble system
[K, F] = assembleSystem(nodes, enodes, edof, De, mms);
% 4. Apply boundary conditions and solve
[U, U_exact] = solveSystem(test, nodes, K, F, mms);
% 5. Compute errors
errors = computeError(nodes, enodes, edof, U, U_exact, De, mms);
% 6. Display results
fprintf('nelx: %.3f\n', test);
fprintf('Relative displacement L2 error:  %.6e\n', errors.rel_l2_error);
fprintf('Relative stress L2 error:  %.6e\n', errors.rel_stress_error);
toc

