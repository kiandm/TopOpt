function [U, U_exact] = solveSystem(~, nodes, K, F, mms)
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
    x = nodes(:,2);    y = nodes(:,3);
    tol = 1e-12;
    boundary_nodes = find( ...       
    abs(x) < tol | abs(x-1) < tol | ...       
    abs(y) < tol | abs(y-1) < tol );    
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