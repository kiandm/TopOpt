function [c, dc, dtheta] = complianceandsensitivities(U, edof, KE0, x, penal, E0, Emin, C, theta, nodes, enodes, gp, w)
    nel = numel(x);
    dc  = zeros(nel,1);
    dtheta = zeros(nel,1);
    c   = 0.0;
    for e = 1:nel
        dofs = edof(e,:);
        Ue   = U(dofs);
        se   = Ue' * KE0{e} * Ue;               % strain energy with unit E
        Esc = Emin + (E0 - Emin) * x(e)^penal;  % element E
        c    = c + Esc * se;
        dc(e)= -(E0 - Emin) * penal * x(e)^(penal-1) * se;

        n  = enodes(e,:);
        xe = nodes(n,2); ye = nodes(n,3);
        dCrot = dRotateElasticMatrix(C, theta(e));   % pass C not De
        dKe_dtheta = assembleElementStiffnessWithFixedC(xe, ye, dCrot, gp, w);
        dtheta(e) = -(Emin + x(e)^penal * (E0-Emin)) * (Ue' * dKe_dtheta * Ue);
    end
end
