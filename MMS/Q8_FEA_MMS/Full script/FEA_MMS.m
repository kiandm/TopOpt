clear; 
clc;

% Material properties
E  = 1;
nu = 0.3;
De = E/(1-nu^2) * [1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2]; %De for plane stress

test = 4;

% Mesh parameters
nelx = test; 
nely = test;
nnx  = nelx + 1; 
nny = nely + 1;                   % Number of nodes in x and y directions
nn   = nnx * nny;                 % Total number of nodes
nel  = nelx * nely;               % Number of elements
ndof = 2 * nn;                    % Number of DoFs
%% 

%------------------------------------

% Assigning node coordinates & DOFs
nodes = zeros(nn, 5);             % [id, x, y, dofUx, dofUy]
node = 0;
for y = 0:nely
    for x = 0:nelx
        node = node + 1;
        nodes(node,1) = node;
        nodes(node,2) = x/nelx;          % x coordinate
        nodes(node,3) = y/nely;          % y coordinate
        nodes(node,4) = 2*node - 1;      % Ux DOF
        nodes(node,5) = 2*node;          % Uy DOF
    end
end
%Defining node_id function with inputs (x,y) and expression (y*(nelx+1) + x + 1)
node_id = @(x,y) (y*(nelx+1) + x + 1);

%%
%------------------------------------
% This could also be wrong. Need to re-derive
 
% Defining exact solutions for Method of Manufactured Solutions (MMS)
% Assuming the displacement field u is trigonometric
ux_exact = @(x,y) sin(2*pi*x) .* sin(2*pi*y);
uy_exact = @(x,y) sin(2*pi*x) .* sin(2*pi*y);

% Strains are given by derivatives of the displacement field
strain_xx = @(x,y) 2*pi*cos(2*pi*x) .* sin(2*pi*y);
strain_yy = @(x,y) 2*pi*sin(2*pi*x) .* cos(2*pi*y);
strain_xy = @(x,y) 2*pi * (cos(2*pi*x) .* sin(2*pi*y) + sin(2*pi*x) .* cos(2*pi*y)); % gamma_xy

% Then calculating corresponding stresses (eqn 5.7 plane stress)
Estar = E/(1 - nu^2);
sigma_xx = @(x,y) 2*pi*Estar * (cos(2*pi*x).*sin(2*pi*y) + nu*sin(2*pi*x).*cos(2*pi*y));
sigma_yy = @(x,y) 2*pi*Estar * (nu*cos(2*pi*x).*sin(2*pi*y) + sin(2*pi*x).*cos(2*pi*y));
% sigma_xy = @(x,y) 2*pi*Estar * ((1 - nu)/2) * (sin(2*pi*x) .* cos(2*pi*y) + cos(2*pi*x) .* sin(2*pi*y));
sigma_xy = @(x,y) pi*Estar*(1-nu) * (cos(2*pi*x).*sin(2*pi*y) + sin(2*pi*x).*cos(2*pi*y));

% Stress derivatives for plane stress
dsigmaxx_dx = @(x,y)  4*pi^2*Estar * (-sin(2*pi*x).*sin(2*pi*y) + nu*cos(2*pi*x).*cos(2*pi*y));
dsigmaxy_dy = @(x,y)  2*pi^2*Estar * ((1-nu) * ((cos(2*pi*x).*cos(2*pi*y) - sin(2*pi*x).*sin(2*pi*y))));
dsigmayy_dy = @(x,y)  4*pi^2*Estar * (nu*cos(2*pi*x).*cos(2*pi*y) - sin(2*pi*x).*sin(2*pi*y));
dsigmaxy_dx = @(x,y)  2*pi^2*Estar * ((1-nu) * ((cos(2*pi*x).*cos(2*pi*y) - sin(2*pi*x).*sin(2*pi*y))));

% Correct body forces (from equilibrium: ∇·σ + b = 0)
body_force_x = @(x,y) -4*pi^2*Estar * (((nu+1)/2)*cos(2*pi*x).*cos(2*pi*y) - ((3-nu)/2)*sin(2*pi*x).*sin(2*pi*y));
body_force_y = @(x,y) -4*pi^2*Estar * (((nu+1)/2)*cos(2*pi*x).*cos(2*pi*y) - ((3-nu)/2)*sin(2*pi*x).*sin(2*pi*y));
%------------------------------------
%%
%------------------------------------

% Element connectivity in terms of node IDs
% enodes defines each row to be an element with the associated node IDs
enodes = zeros(nel, 4); % [n1, n2, n3, n4] as bl, br, tr, tl
elem = 0; %counter to be updated in for loop
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
%%
% DOF connectivity per element
% Defining edof that stores the 8 DoF IDs associated with each element (4
% nodes -> 8 DoF)
edof = zeros(nel, 8); 
for e = 1:nel
    n = enodes(e,:);
    % For all 8 columns in row e, assign them the DoFs 2ni-1 and 2ni
    edof(e,:) = [2*n(1)-1, 2*n(1), 2*n(2)-1, 2*n(2), 2*n(3)-1, 2*n(3), 2*n(4)-1, 2*n(4)];
end
%%
%------------------------------------

% Gauss quadrature (2x2) for bilinear elements
% Define Gauss points at +/-(1/sqrt(3)) with weights w=1
gp = [-1, 1]/sqrt(3);
w  = [1, 1];

% Initialise global stiffness matrix of size ndof x ndof where ndof = 2 * no. of nodes
K = sparse(ndof, ndof);

% Load definition
%top_mid = node_id(nelx/2, nely);
F = sparse(ndof,1);
%F(2*top_mid) = -1;
%%
% Assemble
% Begin by iterating through each element and retrieving the associated 
% DoFs, node numbers, and x,y coords
for e = 1:nel
    dofs = edof(e,:);
    n    = enodes(e,:);
    xe   = nodes(n,2);
    ye   = nodes(n,3);
    % Initialise element stiffness matrix
    KE = zeros(8,8);
    % Initialise element force vector
    FE = zeros(8,1);
    % For each of the 2 gauss points and weights define the natural coords
    for ixi = 1:2
        xi = gp(ixi); wx = w(ixi);
        for ieta = 1:2
            eta = gp(ieta); wy = w(ieta);
            
            % Define the derivatives of the 4 bilinear shape functions
            N = [ (1-xi)*(1-eta)/4, (1+xi)*(1-eta)/4, (1+xi)*(1+eta)/4, (1-xi)*(1+eta)/4 ];

            % dN/dxi and dN/deta (rows: [xi; eta], cols: N1..N4)
            dNdxi = [ (eta-1)/4,  (1-eta)/4,  (eta+1)/4, -(eta+1)/4;
                      (xi -1)/4, -(xi+1)/4,  (xi+1)/4,   (1 -xi)/4 ];

            % Jacobian (2x2)
            J = [dNdxi(1,:)*xe, dNdxi(1,:)*ye;
                 dNdxi(2,:)*xe, dNdxi(2,:)*ye];
            detJ = abs(det(J));

            % dN/dx, dN/dy
            dNdx = J \ dNdxi; % 2x4

            % B matrix (3x8)
            % B=[B1 B2 B3 B4] where each Bi is 3x2 -> 3x8
            B = zeros(3,8); % Initialise
            for i = 1:4    % Loop over the 4 nodes of the element
                % Fill columns for node i: columns (2*i-1) and (2*i) 
                % correspond to Ux_i and Uy_i so for all rows, each 2 columns 
                B(:,2*i-1:2*i) = [ dNdx(1,i),            0;
                                   0,            dNdx(2,i);
                                   dNdx(2,i),  dNdx(1,i) ];
            end

            KE = KE + (B' * De * B) * detJ * (wx * wy);

            % Adding the body forces for MMS
            % Real coordinates of gauss points
            x_gp = N * xe;
            y_gp = N * ye;

            bx = body_force_x(x_gp, y_gp);
            by = body_force_y(x_gp, y_gp);

            % Distributing forces to the nodes
            for i = 1:4
                FE(2*i-1) = FE(2*i-1) + N(i) * bx * detJ * wx * wy;
                FE(2*i)   = FE(2*i)   + N(i) * by * detJ * wx * wy;
            end
        end
    end

    % Assemble into global K
    K(dofs, dofs) = K(dofs, dofs) + KE; % Find the block of K associated 
                                        % with (dofs,dofs) and add KE to it
    F(dofs) = F(dofs) + FE;
end
%%
%------------------------------------

% Boundary conditions
% Retrieve the coordinates of the boundaries of the structure
boundary_nodes = [];
for x = 0:nelx
    boundary_nodes = [boundary_nodes; node_id(x,0); node_id(x,nely)];
end
for y = 1:nely-1  % Avoid duplicating corners
    boundary_nodes = [boundary_nodes; node_id(0,y); node_id(nelx,y)];
end
boundary_nodes = unique(boundary_nodes);
%%
% Initialise a column vector to store the displacements given by the exact
% solution, then populate with calculations for each node
U_exact = zeros(ndof,1);
for i = 1:nn
    % For each node (x,y) calculate the exact displacement
    %Might need to use gauss points instead of node points
    x = nodes(i,2);
    y = nodes(i,3);
    U_exact(2*i-1) = ux_exact(x, y);
    U_exact(2*i)   = uy_exact(x, y);
end
%%
% Create column vector with all DoFs, then create vector containing all free
% DoFs (i.e. difference between vector with all DoFs and fixeddofs)
% Create DOF list
fixeddofs = sort([2*boundary_nodes-1; 2*boundary_nodes]);

%------------------------------------

% freedofs is the set of DoFs I want to solve for
alldofs  = (1:ndof)';
freedofs = setdiff(alldofs, fixeddofs);
fixed_values = U_exact(fixeddofs);

% Reduced stiffness matrix extracts the submatrix corresponding to free DoFs only
% Then extracting the submatrix corresponding to fixed DoFs only
Kff = K(freedofs, freedofs);
Kfc = K(freedofs, fixeddofs);

% Reduced force vector extracts the forces acting on the free DOFs
Ff  = F(freedofs);
%%
%------------------------------------

% Solve
U = zeros(ndof,1); % Initialise displacement vector 
U(fixeddofs) = fixed_values;
U(freedofs) = Kff \ (Ff - Kfc * fixed_values); % Solve K_ff * U_f = F_f

% Compliance
compliance = full(F' * U);
fprintf('Compliance = %.6e\n', compliance);
%%
%------------------------------------

% Then calculate the error vector between the simulation and exact solution
% Might be calculating the error wrong - look at this
error_U = U - U_exact;

% Absolute displacement error at each node 
disp_error = zeros(nn, 1);
for i = 1:nn
    % Obtaining the FE displacements and exact displacements to compare
    ux_fe = U(2*i-1);
    uy_fe = U(2*i);
    ux_ex = U_exact(2*i-1);
    uy_ex = U_exact(2*i);
    disp_error(i) = sqrt((ux_fe - ux_ex)^2 + (uy_fe - uy_ex)^2);
end   
avg_disp_error = mean(disp_error);
%%
%------------------------------------

% Initialize L2 displacement and stress error norm
l2_error_squared = 0;
l2_exact_squared = 0;
stress_l2_error_squared = 0;
stress_exact_squared = 0;
total_area = 0;
% Loop over elements
for e = 1:nel
    dofs = edof(e,:);
    n = enodes(e,:);
    xe = nodes(n,2);
    ye = nodes(n,3);
    
    % Get FEA solution at element nodes
    Ue = U(dofs);  % [u1x, u1y, u2x, u2y, u3x, u3y, u4x, u4y]'
    U_exact_e = U_exact(dofs);
    
    % Loop over Gauss points
    for ixi = 1:2
        xi = gp(ixi); wx = w(ixi);
        for ieta = 1:2
            eta = gp(ieta); wy = w(ieta);
            
            % Shape functions
            N = [ (1-xi)*(1-eta)/4, (1+xi)*(1-eta)/4, ...
                  (1+xi)*(1+eta)/4, (1-xi)*(1+eta)/4 ];
            
            % Shape function derivatives for Jacobian
            dNdxi = [ (eta-1)/4,  (1-eta)/4,  (eta+1)/4, -(eta+1)/4;
                      (xi-1)/4,  -(xi+1)/4,   (xi+1)/4,   (1-xi)/4 ];
            
            % Jacobian
            J = [dNdxi(1,:)*xe, dNdxi(1,:)*ye;
                 dNdxi(2,:)*xe, dNdxi(2,:)*ye];
            detJ = abs(det(J));  % Use absolute value
            total_area = total_area + detJ * wx * wy;
            % dN/dx, dN/dy
            dNdx = J \ dNdxi;

            % B matrix (for stress calculation)
            B = zeros(3,8);
            for i = 1:4
                B(:,2*i-1:2*i) = [ dNdx(1,i),            0;
                                   0,            dNdx(2,i);
                                   dNdx(2,i),  dNdx(1,i) ];
            end
            
            % Physical coordinates
            x_gp = N * xe;
            y_gp = N * ye;
            
            % FEA solution at Gauss point
            ux_fe_gp = N * Ue(1:2:8);  % N1*ux1 + N2*ux2 + N3*ux3 + N4*ux4
            uy_fe_gp = N * Ue(2:2:8);
            
            % Exact solution at Gauss point (evaluate analytic function)
            ux_ex_gp = ux_exact(x_gp, y_gp);
            uy_ex_gp = uy_exact(x_gp, y_gp);

            % Stress calculations also
            strain_fe = B * Ue;
            stress_fe = De * strain_fe;
            stress_ex = [sigma_xx(x_gp, y_gp);
                         sigma_yy(x_gp, y_gp);
                         sigma_xy(x_gp, y_gp)];
            
            % Error at Gauss point
            error_x = ux_fe_gp - ux_ex_gp;
            error_y = uy_fe_gp - uy_ex_gp;

            stress_error = stress_fe - stress_ex;
            stress_error_norm_sq = sum(stress_error.^2);
            stress_exact_norm_sq = sum(stress_ex.^2);
            
            % Accumulate stress L2 error
            stress_l2_error_squared = stress_l2_error_squared + ...
                stress_error_norm_sq * detJ * wx * wy;
            
            stress_exact_squared = stress_exact_squared + ...
                stress_exact_norm_sq * detJ * wx * wy;

            
            % Add to L2 norms (square of error integrated)
            l2_error_squared = l2_error_squared + ...
                (error_x^2 + error_y^2) * detJ * wx * wy;
            
            l2_exact_squared = l2_exact_squared + ...
                (ux_ex_gp^2 + uy_ex_gp^2) * detJ * wx * wy;
        end
    end
end

% Compute relative L2 error
l2_error = sqrt(l2_error_squared);
l2_exact = sqrt(l2_exact_squared);
rel_l2_error = l2_error / l2_exact;

stress_l2_error = sqrt(stress_l2_error_squared);
stress_exact_norm = sqrt(stress_exact_squared);
rel_stress_l2_error = stress_l2_error / stress_exact_norm;
%%
%------------------------------------
nels = [2 4 8 16 32 64 128 256];
uErr = [1.19E+00 2.34E-01 6.24E-02 1.59E-02 4.01E-03 1.00E-03 2.51E-04 6.27E-05];
sigErr = [9.55E-01 4.35E-01 2.24E-01 1.13E-01 5.66E-02 2.83E-02 1.42E-02 7.09E-03];

figure(1); clf;
loglog(nels,uErr,'k-o','MarkerSize',8,'MarkerFaceColor','w'); hold on;
set(gca,'FontSize',16);
%legend('4 noded element','8 noded element');
xlabel('1/h');
ylabel('error')

figure(1); clf;
loglog(nels,sigErr,'k-o','MarkerSize',8,'MarkerFaceColor','w'); hold on;
set(gca,'FontSize',16);
%legend('4 noded element','8 noded element');
xlabel('1/h');
ylabel('error')

uRate   = mean((log(uErr(2:end))-log(uErr(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1))));
sigRate = mean((log(sigErr(2:end))-log(sigErr(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1))));
%------------------------------------
%%
% Outputs
fprintf('nelx: %.3f\n', test);
fprintf('Relative displacement L2 error:  %.6e\n', rel_l2_error);
fprintf('Relative stress L2 error:  %.6e\n', rel_stress_l2_error);

%% ===================== VISUALIZATION =====================
% % -----------------------------
% % Helper settings for plotting
% % -----------------------------
% deform_scale = 0.1;     % scale factor for deformed shape (tweak for visibility)
% line_color   = [0.3 0.3 0.3];
% face_alpha   = 0.6;     % transparency of element faces (0..1)
% 
% % Extract nodal displacement components for convenience
% Ux = U(1:2:end);  % Ux at each node (odd DOFs)
% Uy = U(2:2:end);  % Uy at each node (even DOFs)
% Umag = sqrt(Ux.^2 + Uy.^2);
% 
% % Original coordinates
% X = nodes(:,2);
% Y = nodes(:,3);
% 
% % Deformed coordinates
% Xd = X + deform_scale * Ux;
% Yd = Y + deform_scale * Uy;
% 
% % -----------------------------------------
% % 1) Undeformed vs. deformed mesh (overlay)
% % -----------------------------------------
% figure('Name','Undeformed vs Deformed Mesh','Color','w'); hold on; axis equal;
% title('Undeformed (gray wireframe) and Deformed (blue) mesh');
% xlabel('x'); ylabel('y'); grid on;
% 
% % Plot undeformed mesh as wireframe
% for e = 1:size(enodes,1)
%     n = enodes(e,:);
%     patch('Faces', 1:4, ...
%           'Vertices', [X(n), Y(n)], ...
%           'FaceColor', 'none', ...
%           'EdgeColor', line_color, ...
%           'LineWidth', 0.75);
% end
% 
% % Plot deformed mesh as light-filled quads
% for e = 1:size(enodes,1)
%     n = enodes(e,:);
%     patch('Faces', 1:4, ...
%           'Vertices', [Xd(n), Yd(n)], ...
%           'FaceColor', [0.2 0.6 1.0], ...
%           'FaceAlpha', face_alpha, ...
%           'EdgeColor', 'none');
% end
% 
% % Mark supports (fixed DOFs) and the load node
% if exist('fixeddofs','var')
%     fixed_nodes = unique(ceil(fixeddofs/2)); % map DOFs to node IDs
%     scatter(X(fixed_nodes), Y(fixed_nodes), 40, 'r', 'filled', 'DisplayName','Supports');
% end
% if exist('right_mid','var')
%     scatter(X(right_mid), Y(right_mid), 60, 'k', 'filled', 'DisplayName','Load node');
% end
