function [U] = solveSystem(ndof, K, F, BCs_xy, LC)
    
    BCs_xy = BCs_xy(:);
    LC     = LC(:);

    BCs_y_dofs = 2 * BCs_xy;
    LC_dofs = [2 * LC - 1; 2 * LC];

    fixeddofs = sort([BCs_y_dofs; LC_dofs]);
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