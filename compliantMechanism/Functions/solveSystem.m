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

    fprintf('Starting linear solve, size Kff = %d x %d\n', size(Kff));
    fprintf('rcond(Kff) = %e\n', rcond(full(Kff)));
    tic;
    U(freedofs) = Kff \ (Ff - Kfc * fixed_values);
    fprintf('Solve completed in %.2f seconds\n', toc);
    
    if rcond(full(Kff)) < 1e-12
        warning("Kff is nearly singular: rcond = %g", rcond(Kff));
    end

end