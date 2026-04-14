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
