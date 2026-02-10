% FEA_MMS_Main.m 
tic; clear; clc;
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
% 6. Initialize design
x = volfrac * ones(nel,1);
% 7. Optimization loop
for it = 1:maxit
    xold = x;
    % ----- Assemble K(x) and F -------------
    [K, F] = assembleSystem(nodes, enodes, ndof, edof, KE0, x, penal, E0, Emin);
    % ----- Apply loads ---------------------
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
    % ----- Convergence, logging -------------
    change = max(abs(x - xold));
    fprintf('It. %4d | Obj (c): %12.5e | Vol: %6.4f | chg: %8.4e\n', ...
            it, c, mean(x), change);
    % ----- Iteration history ----------------
    iterationHistory(it, :) = [it, c, mean(x), change];
    % ----- Quick plot -----------------------
    if mod(it,5)==1 || change<tol
        plotDensityMesh(nodes, enodes, x);
        drawnow;
    end
    if change < tol
        break;
    end
end
% plot iteration convergence history
figure;
plot(iterationHistory(1:it, 1), iterationHistory(1:it, 2), '-o');
xlabel('Iteration');
ylabel('Objective Function (Compliance)');
title('Convergence History');
grid on;
toc

function [nn,nodes,nel,enodes,ndof,edof] = createMesh()
% Imports list of node coordinates, as well as element construction from
% Cubit, then constructs connectivity.
addpath('C:\Users\sfxz32\OneDrive - Durham University\Desktop\Cubit')
nodes_file    = readmatrix("lbeam_nodes.dat");      % [id, x, y]
elements_file = readmatrix("lbeam_elements.dat");   % [eid, n1, n2, n3, n4]
% Node coordinates & DoFs
nn = size(nodes_file, 1);           % Number of nodes
nodes = zeros(nn, 5);               % [id, x, y, dofUx, dofUy]
nodes(:,1) = nodes_file(:,1);       % id
nodes(:,2) = nodes_file(:,2);       % x
nodes(:,3) = nodes_file(:,3);       % y
nodes(:,4) = 2*nodes(:,1) - 1;      % Ux DOF = 2*id - 1
nodes(:,5) = 2*nodes(:,1);          % Uy DOF = 2*id
% Element connectivity
nel = size(elements_file, 1);
enodes = elements_file(:, 2:5);     % take columns n1 n2 n3 n4
edof = zeros(nel, 8);
for e = 1:nel
    n = enodes(e,:);
    edof(e,:) = [2*n(1)-1, 2*n(1), ...
                 2*n(2)-1, 2*n(2), ...
                 2*n(3)-1, 2*n(3), ...
                 2*n(4)-1, 2*n(4)];
end
% Total DOFs
ndof = 2 * nn;
end

function [BCs_xy,loads_xy] = loadConditions()
% Defines two nodesets from Cubit ABAQUS file, containing nodes for BCs,
% and nodes for loads.

BCs_xy = [
       2,     182,     218,     219,     220,     221,...
     222,     223,     224,     225,     226,     227,...
     228,     229,     230,     231,     232,     233,...
     234,     235,     236,     237,     238,     239,...
     240
     ];

loads_xy = [122];

end

function KE = assembleElementStiffnessQ4(xe, ye, De, gp, w)
    KE = zeros(8,8);
    for ixi = 1:2
        xi = gp(ixi); wx = w(ixi);
        for ieta = 1:2
            eta = gp(ieta); wy = w(ieta);
            [N, dNdxi] = shapeFcnQ4(xi, eta);
            % Jacobian & inverse
            J = [dNdxi(1,:)*xe, dNdxi(1,:)*ye;
                 dNdxi(2,:)*xe, dNdxi(2,:)*ye];
            detJ = det(J);          % determinant
            invJ = J \ eye(2);
            dNdx = invJ * dNdxi;    % 2x4
            % B matrix (3x8)
            B = zeros(3,8);
            for a = 1:4
                B(:,2*a-1:2*a) = [ dNdx(1,a), 0;
                                   0,         dNdx(2,a);
                                   dNdx(2,a), dNdx(1,a) ];
            end
            KE = KE + (B' * De * B) * abs(detJ) * (wx * wy);
        end
    end
end

function [Wi, Wj, Ww, Wsum] = SensitivityFilter(ecent, rmin)
% Returns sparse weight triplets (i,j,wij) with fac = rmin - dist >= 0
% and Wsum(i) = sumj wij for normalization.
    nel = size(ecent,1);
    rij = []; cij = []; wij = [];
    for i = 1:nel
        pi = ecent(i,:);
        for j = 1:nel
            pj = ecent(j,:);
            dist = norm(pi - pj);
            fac  = rmin - dist;
            if fac > 0
                rij(end+1,1) = i; 
                cij(end+1,1) = j; 
                wij(end+1,1) = fac; 
            end
        end
    end
    Wi = rij; Wj = cij; Ww = wij;
    % sum of weights per row i
    Wsum = accumarray(Wi, Ww, [nel,1], @sum);
end

function [K, F] = assembleSystem(nodes, enodes, ndof, edof, KE0, x, penal, E0, Emin)
    K = sparse(ndof, ndof);
    F = sparse(ndof, 1);
    

    % Scale each element KE0 by Escale(x)
    Escale = Emin + (E0 - Emin) * (x.^penal);  % nel x 1

    for e = 1:size(enodes, 1)
        dofs = edof(e,:);
        K(dofs, dofs) = K(dofs, dofs) + Escale(e) * KE0{e};
    end

end

function [U] = solveSystem(ndof, K, F, BCs_xy)

    fixeddofs = sort([2*BCs_xy-1; 2*BCs_xy]);
    freedofs = setdiff(1:ndof, fixeddofs);
    fixed_values = zeros(numel(fixeddofs), 1);

    % Solve
    Kff = K(freedofs, freedofs);
    Kfc = K(freedofs, fixeddofs);
    Ff = F(freedofs);
    
    U = zeros(ndof, 1);
    U(fixeddofs) = fixed_values;
    U(freedofs) = Kff \ (Ff - Kfc * fixed_values);
    
end

function [c, dc] = complianceandsensitivities(U, edof, KE0, x, penal, E0, Emin)
    nel = numel(x);
    dc  = zeros(nel,1);
    c   = 0.0;
    for e = 1:nel
        dofs = edof(e,:);
        Ue   = U(dofs);
        se   = Ue' * KE0{e} * Ue;               % strain energy with unit E
        Esc = Emin + (E0 - Emin) * x(e)^penal;  % element E
        c    = c + Esc * se;
        dc(e)= -(E0 - Emin) * penal * x(e)^(penal-1) * se;
    end
end

function dcf = applySensitivityFilter(dc, x, Wi, Wj, Ww, Wsum)
% dĉi = [sumj wij * xj * dcj] / [xi * sumj wij]
    nel = numel(x);
    num = zeros(nel,1);
    for k = 1:numel(Ww)
        i = Wi(k); j = Wj(k); w = Ww(k);
        num(i) = num(i) + w * x(j) * dc(j);
    end
    dcf = num ./ (x .* Wsum);
end

function xnew = OC_update(x, dc, volfrac)
    l1 = 0; l2 = 1e9; move = 0.2;
    xnew = x;
    while (l2 - l1) / (l2 + l1) > 1e-6
        lmid = 0.5*(l1 + l2);
        xcand = max(0.001, max(x - move, min(1.0, min(x + move, x .* sqrt(-dc ./ lmid)))));
        if mean(xcand) - volfrac > 0
            l1 = lmid;
        else
            l2 = lmid;
        end
        xnew = xcand;
    end
end

function plotDensityMesh(nodes, enodes, x)
    % simple patch coloring by x
    clf;
    hold on;
    colormap(gray); % black=1 (solid), white=0 (void) if we invert later
    for e = 1:size(enodes,1)
        n  = enodes(e,:);
        xy = [nodes(n,2), nodes(n,3)];
        c = 1 - max(0,min(1,x(e)));
        c = [c c c];
        patch('Faces', [1 2 3 4], ...
              'Vertices', xy, ...
              'FaceColor', c, ...
              'EdgeColor', [0.7 0.7 0.7], ...
              'LineWidth', 0.25);
    end
    axis equal tight off;
    title(sprintf('Density (vol=%.3f)', mean(x)));
    hold off;
end


