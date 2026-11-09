function T = rebuildThetaSweep(agg, ndiv, folder, expected)
% rebuildThetaSweep  Rebuild a theta-sweep table from the results files the sweep already saved, without re-running anything.
%   T = rebuildThetaSweep('smoothmax', 200)      every Results/opt_results_*.mat with opts.agg = 'smoothmax' and opts.ndiv = 200
%   plotThetaSweep(T, maxiter=2000)
% Same columns as thetaSweep.m. A smoothmax run has one constraint and one index: they go in g_ft and maxFI_ft, the rest are NaN.
% If an angle appears more than once the newest file is used, so a later smoothmax test at this ndiv would replace it:
% narrow the line marked FILTER if the folder holds other runs.
if nargin < 3, folder = 'Results'; end
if nargin < 4, expected = 5:5:175; end
vars = {'theta_deg','theta_rad','g_ft','g_fc','g_mt','g_mc','maxFI_ft','maxFI_fc','maxFI_mt','maxFI_mc', ...
        'g_c','c','V','M','iter','dV_win','C0','runtime'};
d = dir(fullfile(folder, 'opt_results_*.mat'));
D = nan(0, numel(vars));  flag = {};  file = {};  stamp = [];
for k = 1:numel(d)
    f = fullfile(d(k).folder, d(k).name);
    try
        o = load(f, 'opts');  o = o.opts;
        if ~strcmp(o.agg, agg) || o.ndiv ~= ndiv || ~strcmp(o.problem, 'Lbrac60'), continue; end   % FILTER
        S = load(f, 'xphy', 'numele', 'FailIdx', 'g_hs', 'c', 'V', 'M', 'iter', 'history', 'runtime', 'exitflag');
        solid = S.xphy(1:S.numele) > 0.3;                          % same mask as runCase.m
        g  = nan(1,4);  g(1:numel(S.g_hs)) = S.g_hs(:)';
        mx = max(S.FailIdx(solid,:), [], 1);  FI = nan(1,4);  FI(1:numel(mx)) = mx;
        w  = min(o.win, S.iter);  vl = S.history(S.iter-w+1:S.iter, 3);
        D(end+1,:) = [round(rad2deg(o.theta_init), 6), o.theta_init, g, FI, S.c/o.c_max - 1, S.c, S.V, S.M, S.iter, ...
                      (max(vl) - min(vl))/vl(end), S.history(1,2), S.runtime];
        flag{end+1,1} = S.exitflag;  file{end+1,1} = d(k).name;  stamp(end+1,1) = d(k).datenum;
    catch ME
        fprintf('skipped %s: %s\n', d(k).name, ME.message);       % e.g. an old file without opts
    end
end
if isempty(D), error('rebuildThetaSweep:none', 'No %s results with ndiv = %d in %s.', agg, ndiv, folder); end
[~, ord] = sortrows([D(:,1), stamp]);                              % by angle, then by time
D = D(ord,:);  flag = flag(ord);  file = file(ord);
last = [D(2:end,1) ~= D(1:end-1,1); true];                         % keep the newest file of each angle
D = D(last,:);  flag = flag(last);  file = file(last);
T = [array2table(D, 'VariableNames', vars), table(flag, file, 'VariableNames', {'exitflag', 'file'})];
missing = setdiff(expected, T.theta_deg);
if ~isempty(missing), warning('rebuildThetaSweep:missing', 'No results file for theta_0 = %s deg.', mat2str(missing(:).')); end
out = fullfile(folder, sprintf('theta_sweep_%s_%s.csv', agg, datestr(now, 'yyyymmdd_HHMMSS')));
writetable(T, out);
fprintf('%d angles rebuilt -> %s\n', height(T), out);
end