function [K, F] = assembleSystem(nodes, enodes, ndof, edof, KE0, x, penal, E0, Emin)
    K = sparse(ndof, ndof);
    F = sparse(ndof, 1);
    

    % Scale each element KE0 by Escale(x)
    % Escale = Emin + (E0 - Emin) * (x.^penal);  % nel x 1

    % for e = 1:size(enodes, 1)
    %     dofs = edof(e,:);
    %     K(dofs, dofs) = K(dofs, dofs) + Escale(e) * KE0{e};
    % end
    nel = size(edof,1);
    Escale = Emin + (E0-Emin)*x.^penal;
    KEall  = cell2mat(cellfun(@(k) k(:), KE0', 'UniformOutput', false));   % 64 x nel
    sK = reshape(KEall .* Escale', 64*nel, 1);
    iK = reshape(kron(edof, ones(8,1))', 64*nel, 1);
    jK = reshape(kron(edof, ones(1,8))', 64*nel, 1);
    K  = sparse(iK, jK, sK, ndof, ndof);  F = sparse(ndof, 1);

end
