
function topopt_modular()
% TOPOLOGY OPTIMIZATION (SIMP) — Modular template built on your FEA style
% -------------------------------------------------------------
% Goal: Minimize compliance subject to a volume fraction constraint.
% Mesh: Structured grid of Q4 elements with unit spacing.
% Material: Plane stress, unit thickness.
% -------------------------------------------------------------

clear; clc;

%% ------------------------ User Parameters ------------------------
nelx   = 60;               % #elements in x
nely   = 20;               % #elements in y
volfrac= 0.5;              % target volume fraction (0..1)
penal  = 3.0;              % SIMP penalty (pushes x toward 0/1)
rmin   = 1.5;              % filter radius (mesh-independence)
move   = 0.2;              % OC move limit per iteration (stability)
E      = 1.0;              % Young's modulus
nu     = 0.3;              % Poisson's ratio

% Plane stress constitutive matrix
De = E/(1-nu^2) * [ 1, nu, 0;
                    nu, 1, 0;
                    0,  0, (1-nu)/2 ];

% Gauss points (2x2) for Q4 KE precomputation
gp = [-1, 1]/sqrt(3);
w  = [1, 1];

%% ------------------------ Build Mesh & DOFs ------------------------
% Build nodes, element connectivity, DOF arrays, and node_id helper
[nodes, enodes, edof, node_id, nn, ndof] = buildMesh(nelx, nely);

% Precompute solid element stiffness KE (unit square Q4, plane stress)
KE = elementKe_Q4(De, gp, w);

%% ------------------------ TopOpt Initialization ------------------------
% Design variable: element density x(ely,elx)
x      = volfrac * ones(nely, nelx); % Start from uniform design
change = 1.0;                        % Iteration change measure
iter   = 0;

%% ------------------------ Optimization Loop ------------------------
while change > 0.01
    iter = iter + 1;
    xold = x;

    % --- Assemble global K with SIMP scaling (x^penal) ---
    K = assembleGlobalK(nelx, nely, ndof, edof, x, penal, KE);

    % --- Define loads and boundary conditions  ---
    [F, fixeddofs, right_mid] = defineLoadsAndBC(nelx, nely, ndof, node_id);

    % --- Solve FE: K_ff * U_f = F_f ---
    U = solveFE(K, F, fixeddofs);

    % --- Compliance and sensitivities (dc) ---
    [c, dc] = complianceAndSensitivity(nelx, nely, edof, U, KE, x, penal);

    % --- Filter sensitivities (cone filter) ---
    dcf = sensitivityFilter(dc, x, rmin);

    % --- Update densities via Optimality Criteria (OC) ---
    x = ocUpdate(x, dcf, volfrac, move);

    % --- Measure change & report ---
    change = max(abs(x(:) - xold(:)));
    fprintf('Iter %3d: Obj (compliance) = %10.4f | Vol = %6.3f | Change = %6.3f\n', ...
            iter, c, mean(x(:)), change);

    % --- Visualize current design ---
    visualizeDensity(x);
end

fprintf('Optimization finished. Final compliance: %.6f, volume: %.3f\n', c, mean(x(:)));

end % ====== topopt_modular() ======


%% ================== Mesh Builder ==================
function [nodes, enodes, edof, node_id, nn, ndof] = buildMesh(nelx, nely)
% Builds nodes table, element-node connectivity, element-DOF connectivity,
% and a node_id helper aligned with your numbering (y outer, x inner).

nnx  = nelx + 1; 
nny  = nely + 1;
nn   = nnx * nny;    % number of nodes
ndof = 2 * nn;       % Ux, Uy per node

% Node table: [id, x, y, dofUx, dofUy]
nodes = zeros(nn, 5);
node  = 0;
for y = 0:nely
    for x = 0:nelx
        node = node + 1;
        nodes(node,1) = node;
        nodes(node,2) = x;          % unit spacing: x coordinate
        nodes(node,3) = y;          % unit spacing: y coordinate
        nodes(node,4) = 2*node - 1; % Ux DOF
        nodes(node,5) = 2*node;     % Uy DOF
    end
end

% Node id helper: (x,y) → node ID
node_id = @(x,y) (y*(nelx+1) + x + 1);

% Element-node connectivity (each row: [n1 n2 n3 n4] = [bl br tr tl])
enodes = zeros(nelx*nely, 4);
elem   = 0;
for ex = 1:nelx
    for ey = 1:nely
        elem = elem + 1;
        x0 = ex-1; y0 = ey-1;
        n1 = node_id(x0,   y0  ); % bottom-left
        n2 = node_id(x0+1, y0  ); % bottom-right
        n3 = node_id(x0+1, y0+1); % top-right
        n4 = node_id(x0,   y0+1); % top-left
        enodes(elem,:) = [n1, n2, n3, n4];
    end
end

% Element-DOF connectivity (each row: 8 DOFs for Q4)
edof = zeros(nelx*nely, 8);
for e = 1:(nelx*nely)
    n = enodes(e,:);
    edof(e,:) = [2*n(1)-1, 2*n(1), 2*n(2)-1, 2*n(2), ...
                 2*n(3)-1, 2*n(3), 2*n(4)-1, 2*n(4)];
end
end


%% ================== Element KE (Q4, plane stress) ==================
function KE = elementKe_Q4(De, gp, w)
% Precomputes the solid element stiffness KE for a unit square Q4 element.
% Uses 2x2 Gauss quadrature. Assumes axis-aligned uniform mesh (unit size).

KE = zeros(8,8);

% Nodal coords of a unit element (local instance; global elements share size)
xe = [0; 1; 1; 0];
ye = [0; 0; 1; 1];

for ixi = 1:2
    xi = gp(ixi); wx = w(ixi);
    for ieta = 1:2
        eta = gp(ieta); wy = w(ieta);

        % dN/dxi, dN/deta (rows: [xi; eta], cols: N1..N4)
        dNdxi = [ (eta-1)/4,  (1-eta)/4,  (eta+1)/4, -(eta+1)/4;
                  (xi -1)/4, -(xi+1)/4,  (xi+1)/4,   (1 -xi)/4 ];

        % Jacobian J from (xi,eta) to (x,y)
        J = [dNdxi(1,:)*xe, dNdxi(1,:)*ye;
             dNdxi(2,:)*xe, dNdxi(2,:)*ye];

        % Derivatives in physical coordinates: [dN/dx; dN/dy]
        dNdx = J \ dNdxi;

        % Build B (3x8)
        B = zeros(3,8);
        for i = 1:4
            B(:,2*i-1:2*i) = [ dNdx(1,i), 0;
                               0,          dNdx(2,i);
                               dNdx(2,i),  dNdx(1,i) ];
        end

        % Accumulate KE
        KE = KE + (B' * De * B) * det(J) * (wx * wy);
    end
end
end


%% ================== Assembly (SIMP scaling) ==================
function K = assembleGlobalK(nelx, nely, ndof, edof, x, penal, KE)
% Builds sparse global K from precomputed solid KE scaled by x^penal.

K = sparse(ndof, ndof);
elem = 0;
for ex = 1:nelx
    for ey = 1:nely
        elem = elem + 1;
        dofs = edof(elem,:);
        scale = x(ey,ex)^penal; % SIMP stiffness scaling
        K(dofs,dofs) = K(dofs,dofs) + scale * KE;
    end
end
end


%% ================== Loads & Boundary Conditions ==================
function [F, fixeddofs, right_mid] = defineLoadsAndBC(nelx, nely, ndof, node_id)
% Your MBB-like supports & load (kept consistent with your code style).

% Corner nodes (for clarity; not directly used below)
left_bottom  = node_id(0, 0);
left_top     = node_id(0, nely);
right_bottom = node_id(nelx, 0);
right_top    = node_id(nelx, nely);

% Fix Uy at both left corners and Ux at both left corners (your choice)
fixeddofs = [ 2*left_bottom-1;   % Ux(left-bottom)
              2*left_bottom;     % Uy(left-bottom)
              2*left_top;        % Uy(left-top)
              2*left_top-1 ];    % Ux(left-top)

% Force: downward unit load at mid of right edge
right_mid = node_id(nelx, nely/2);

F = sparse(ndof,1);
F(2*right_mid) = -1; % Apply -1 in Uy at right_mid
end


%% ================== FE Solve ==================
function U = solveFE(K, F, fixeddofs)
% Partitions DOFs into free/fixed and solves the reduced system.

ndof    = size(K,1);
alldofs = (1:ndof)';
freedofs= setdiff(alldofs, fixeddofs);

U = zeros(ndof,1);
U(freedofs) = K(freedofs,freedofs) \ F(freedofs);
% Fixed DOFs remain zero (prescribed displacement = 0)

end


%% ================== Compliance & Sensitivities ==================
function [c, dc] = complianceAndSensitivity(nelx, nely, edof, U, KE, x, penal)
% Computes compliance c and sensitivities dc (per element).
% c = sum x^penal * Ue' * KE * Ue
% dc = -penal * x^(penal-1) * (Ue' * KE * Ue)

c  = 0.0;
dc = zeros(nely, nelx);

elem = 0;
for ex = 1:nelx
    for ey = 1:nely
        elem = elem + 1;
        dofs = edof(elem,:);
        Ue   = U(dofs);
        ce   = Ue' * KE * Ue;              % strain energy for solid KE
        c    = c + x(ey,ex)^penal * ce;    % compliance contribution
        dc(ey,ex) = -penal * x(ey,ex)^(penal-1) * ce; % sensitivity
    end
end
end


%% ================== Sensitivity Filter (cone) ==================
function dcf = sensitivityFilter(dc, x, rmin)
% Mesh-independence cone filter:
% dcf(i,j) = (sum_{k,l in neighborhood} w * x(k,l) * dc(k,l)) / (x(i,j) * sum_w)
% where w = max(0, rmin - distance)

[nely, nelx] = size(x);
dcf = zeros(nely, nelx);

for i = 1:nelx
    for j = 1:nely
        wsum = 0.0;
        val  = 0.0;
        i1 = max(i - floor(rmin), 1);
        i2 = min(i + floor(rmin), nelx);
        j1 = max(j - floor(rmin), 1);
        j2 = min(j + floor(rmin), nely);
        for k = i1:i2
            for l = j1:j2
                dist = sqrt((i-k)^2 + (j-l)^2);
                w = max(0, rmin - dist);
                wsum = wsum + w;
                val  = val + w * x(l,k) * dc(l,k);
            end
        end
        dcf(j,i) = val / max(1e-12, (x(j,i) * wsum)); % safe divide
    end
end
end


%% ================== OC Update ==================
function xnew = ocUpdate(x, dc, volfrac, move)
% Optimality Criteria update with bisection on lambda to meet volume constraint.
% xnew = clamp( x * sqrt(-dc / lambda), with move limits and bounds [0.001, 1] )

l1 = 0; l2 = 1e12;        % bounds for lambda
xnew = x;                 % initialize

while (l2 - l1) > 1e-4
    lmid = 0.5 * (l2 + l1);

    % Candidate update: x * sqrt(-dc / lmid)
    % Apply move limits and bounds
    xcand = x .* sqrt( max(0, -dc) / max(1e-12, lmid) );
    xcand = max(0.001, min(1.0, min(x + move, max(x - move, xcand))));

    % Volume check: mean(xcand) ? volfrac
    if mean(xcand(:)) - volfrac > 0
        % Too much material → increase lambda to shrink x
        l1 = lmid;
    else
        % Too little material → decrease lambda to grow x
        l2 = lmid;
    end
    xnew = xcand;
end
end


%% ================== Visualization ==================
function visualizeDensity(x)
% Simple grayscale density plot (black ~ solid, white ~ void)
colormap(gray);
imagesc(1 - x);           % invert for dark-solid
axis equal; axis tight; axis off;
drawnow;
end
