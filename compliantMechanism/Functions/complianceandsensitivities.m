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
