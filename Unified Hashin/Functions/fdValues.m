function [c, g] = fdValues(xv, H, Hs, beta, eta, penal, numnode, numele, gs, edofMat, coords, conn, freedofs, F, matprop, dphix_ref, dphiy_ref, strength, opts)
% Values only (no sensitivities) of compliance c and failure constraint(s) g at design variables xv, via the same
% filter -> FE -> Hashin chain as the FD script. Used by the directional (Taylor) test.
xphy = apply_filter_and_projection(xv, numele, H, Hs, beta, eta);
[U, ~, KE0, dK] = FE_analysis(xphy, penal, numnode, numele, gs, edofMat, coords, conn, freedofs, F, matprop, dphix_ref, dphiy_ref);
c = objective_function(U, xphy, penal, numele, gs, edofMat, coords, conn, matprop, dphix_ref, dphiy_ref);
g = Hashin(U, dK, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs, dphix_ref, dphiy_ref, opts);
end