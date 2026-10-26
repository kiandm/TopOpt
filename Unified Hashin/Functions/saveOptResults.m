function filename = saveOptResults(xphy, U, numele, numnode, gs, edofMat, coords, conn, ...
    matprop, strength, penal, dphix_ref, dphiy_ref, freedofs, F, FailIdx, g_hs, vonMises, iter, M, ...
    opts, exitflag, history, c, V, runtime)
% Bundles the final state of one run - the design, the displacement field it produced, the mesh/problem
% data needed to recompute stresses from them, the Hashin results, AND everything needed to repeat the run
% (options, git commit) - into one struct, and saves it to Results/<filename>.
%
% The file name carries a timestamp first (checkTsaiWu.m picks the newest file by sorting on it), then the
% options that matter and the commit.
%
% Load the result back with:  r = load(fullfile('Results', filename));

[commit, dirty] = getGitInfo();

seedtag = '';  if ~isempty(opts.seed), seedtag = sprintf('_s%d', opts.seed); end
projtag = '';  if strcmp(opts.projection, 'none'), projtag = '_noproj'; end   % projection off: say so in the name
filename = sprintf('opt_results_%s_%s_P%g_c%g%s%s_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'), ...
                   opts.agg, opts.p, opts.c_max, projtag, seedtag, commit);
if dirty == 1, filename = strrep(filename, '.mat', '_dirty.mat'); end   % uncommitted edits: not reproducible from the hash

filepath = fullfile(pwd, 'Results');
if ~exist(filepath, 'dir'), mkdir(filepath); end

results.xphy      = xphy;       % final design: [density(1:numele); theta(numele+1:end)]
results.U         = U;          % displacement field of that design
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
results.FailIdx   = FailIdx;    % element-averaged stress-ratio index: numele x 4 (separated) or numele x 1 (smoothmax)
results.g_hs      = g_hs;       % final p-norm constraint(s): [g_ft;g_fc;g_mt;g_mc] or g_hs
results.vonMises  = vonMises;   % per-element von Mises stress
results.iter      = iter;       % iterations run
results.M         = M;          % final greyness measure [%]
% ---- added for reproducibility and comparison ----
results.opts      = opts;       % every option of the run (setOpts.m)
results.exitflag  = exitflag;   % 'converged' (stopping rule met) or 'maxiter'
results.history   = history;    % iterations x 10: [iter, c, V, change, g(1:4), M, seconds]; g padded with NaN for smoothmax
results.c         = c;          % final compliance
results.V         = V;          % final volume fraction
results.runtime   = runtime;    % wall-clock seconds
results.commit    = commit;     % git commit hash of the code
results.dirty     = dirty;      % 1 = tracked files had uncommitted changes when the run started
results.date      = datestr(now);

save(fullfile(filepath, filename), '-struct', 'results');
fprintf('Saved optimisation results to %s\n', fullfile('Results', filename));
end