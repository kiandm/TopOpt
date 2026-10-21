function [g, dg_dx, dg_dtheta, FailIdx, vonMises] = Hashin( ...
    U, dK, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs, dphix_ref, dphiy_ref, opts)
% HASHIN  Failure-constraint dispatcher.  opts.agg selects how the Hashin indices are aggregated:
%   'separated'  four mode indices, each its own p-norm over the elements -> g is 4 x 1, FailIdx is numele x 4
%   'smoothmax'  Dong-type: smooth max of fibre/matrix indices, ONE p-norm -> g is 1 x 1, FailIdx is numele x 1
% All numerical constants come from opts (see setOpts.m).
if nargin < 16 || isempty(opts)
    error('Hashin:opts', 'Hashin needs the options struct as its last argument, e.g. opts = setOpts();');
end
switch opts.agg
    case 'separated'
        [g, dg_dx, dg_dtheta, FailIdx, vonMises] = Hashin_separated( ...
            U, dK, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs, dphix_ref, dphiy_ref, opts);
    case 'smoothmax'
        [g, dg_dx, dg_dtheta, FailIdx, vonMises] = Hashin_smoothmax( ...
            U, dK, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs, dphix_ref, dphiy_ref, opts);
    otherwise
        error('Hashin:agg', 'Unknown opts.agg "%s".', opts.agg);
end
end