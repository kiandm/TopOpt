function [c, dc, dtheta] = complianceandsensitivities(E1,E2,nu12,nu21,G12,U, edof, KE0, x, penal, E0, Emin, De, theta,nodes, enodes, gp, w)
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
        dDe_dtheta = dRotateElasticMatrix(De, theta(e));
        dKe_dtheta = assembleElementStiffnessQ4(e,nel,xe,ye,dDe_dtheta,gp,w, E1,E2,nu12,nu21,G12,theta);

        dtheta(e) = -(Emin + x(e)^penal * (E0-Emin)) * (Ue' * dKe_dtheta * Ue);
    end
end
