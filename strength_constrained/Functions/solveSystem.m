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