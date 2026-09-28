function filename = saveOptResults(xphy, U, numele, numnode, gs, edofMat, coords, conn, ...
    matprop, strength, penal, dphix_ref, dphiy_ref, freedofs, F, FailIdx, g_hs, vonMises, iter, M, filename)
% Bundles the converged state of an optimisation run - the design, the
% displacement field it produced, the mesh/problem data needed to recompute
% stresses from them, and the Hashin results already computed for this design
% - into a single struct and saves it to <filename>.mat. This lets the design
% be re-checked against a different failure criterion later (see
% checkTsaiWu.m) without re-running the optimisation.
%
% filename is optional; if omitted, a timestamped name is generated so
% repeated runs don't overwrite each other.
%
% Load the result back with:  r = load(filename);

if nargin < 21 || isempty(filename)
    filename = sprintf('opt_results_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
end

results.xphy      = xphy;       % converged design: [density(1:numele); theta(numele+1:end)]
results.U         = U;          % converged global displacement field
results.numele    = numele;
results.numnode   = numnode;
results.gs        = gs;         % Gauss point table (gauss_domain.m)
results.edofMat   = edofMat;
results.coords    = coords;
results.conn      = conn;
results.matprop   = matprop;    % E1, E2, nu12, nu21, G12
results.strength  = strength;   % Xt, Xc, Yt, Yc, S
results.penal     = penal;
results.dphix_ref = dphix_ref;
results.dphiy_ref = dphiy_ref;
results.freedofs  = freedofs;
results.F         = F;          % applied load vector
results.FailIdx   = FailIdx;    % converged per-element, per-mode Hashin index (numele x 4)
results.g_hs      = g_hs;       % converged Hashin p-norm constraints [g_ft;g_fc;g_mt;g_mc]
results.vonMises  = vonMises;   % converged per-element von Mises stress
results.iter      = iter;       % iterations taken to converge
results.M         = M;          % final greyness measure

save(filename, '-struct', 'results');
fprintf('Saved optimisation results to %s\n', filename);
end