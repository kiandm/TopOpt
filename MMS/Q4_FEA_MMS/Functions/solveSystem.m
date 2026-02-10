function [U, U_exact] = solveSystem(test, nodes, K, F, mms)
    ndof = size(K, 1);
    nn = size(nodes, 1);
    
    % Compute exact solution at nodes
    U_exact = zeros(ndof, 1);
    for i = 1:nn
        x = nodes(i,2); 
        y = nodes(i,3);
        U_exact(2*i-1) = mms.ux_exact(x, y);
        U_exact(2*i) = mms.uy_exact(x, y);
    end
    
    % Boundary conditions (all boundaries fixed)
    %nelx = max(nodes(:,2)) * (size(nodes,1)/(max(nodes(:,3))+1) - 1);
    %nely = max(nodes(:,3)) * (size(nodes,1)/(max(nodes(:,2))+1) - 1);
    nelx = test;
    nely = test;
    node_id = @(x,y) (y*(nelx+1) + x + 1);
    
    boundary_nodes = [];
    for x = 0:nelx
        boundary_nodes = [boundary_nodes; node_id(x,0); node_id(x,nely)];
    end
    for y = 1:nely-1
        boundary_nodes = [boundary_nodes; node_id(0,y); node_id(nelx,y)];
    end
    boundary_nodes = unique(boundary_nodes);
    
    fixeddofs = sort([2*boundary_nodes-1; 2*boundary_nodes]);
    freedofs = setdiff(1:ndof, fixeddofs);
    fixed_values = U_exact(fixeddofs);
    
    % Solve
    Kff = K(freedofs, freedofs);
    Kfc = K(freedofs, fixeddofs);
    Ff = F(freedofs);
    
    U = zeros(ndof, 1);
    U(fixeddofs) = fixed_values;
    U(freedofs) = Kff \ (Ff - Kfc * fixed_values);
end