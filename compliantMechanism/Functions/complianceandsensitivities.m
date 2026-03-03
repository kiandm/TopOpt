function [c, dc] = complianceandsensitivities(U_in, U_out, F_out, edof, KE0, x, penal, E0, Emin)
    nel = numel(x);
    dc  = zeros(nel,1);
    c   = 0.0;
    for e = 1:nel
        dofs = edof(e,:);
        
        Ue_in   = U_in(dofs);
        Ue_out   = U_out(dofs);
        
        dKe = penal*x(e)^(penal-1)*(KE0{e}*(E0-Emin));

        %se   = Ue' * KE0{e} * Ue;               % strain energy with unit E
        %Esc = Emin + (E0 - Emin) * x(e)^penal;  % element E
        c    = F_out' * U_in;
        dc(e)= - Ue_out' * dKe * Ue_in;
    end
end
