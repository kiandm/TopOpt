addpath('Functions');  if ~exist('Results','dir'), mkdir('Results'); end
%% Sweep settings (edit here)
aggs    = {'separated', 'smoothmax'};                      % mode-separated (4 constraints) and mode-collapsed (one smooth-max constraint)
th_deg  = 5:5:175;                                         % not 0 or 180: those start exactly on the angle bounds
ndiv    = 150;                                             % elements along the 100 mm side, a multiple of 5 (3600 elements at 75, 14400 at 150)
maxiter = 2000;                                            % the default in setOpts.m is 1000
% Run normally: every (aggregation, angle) pair one after another; one CSV per aggregation, Results/theta_sweep_<agg>_<time>.csv.
% Run as a Slurm job array (hamilton_thetaSweep.sh): task k runs ONE pair and writes Results/theta_sweep_part_<agg>_th<deg>.csv;
% mergeThetaSweep.m joins the parts of one aggregation afterwards. Tasks 1..35 are aggs{1} at th_deg(1..35), tasks 36..70 are aggs{2}.
nT     = numel(th_deg);  nJobs = numel(aggs) * nT;
task   = str2double(getenv('SLURM_ARRAY_TASK_ID'));        % NaN when the variable is not set
isTask = ~isnan(task);
if isTask
    if task < 1 || task > nJobs || task ~= round(task), error('thetaSweep:task', 'SLURM_ARRAY_TASK_ID = %g, but there are %d runs: use --array=1-%d', task, nJobs, nJobs); end
    jobs = task;
else
    jobs = 1:nJobs;
end
vars = {'theta_deg','theta_rad','g_ft','g_fc','g_mt','g_mc','maxFI_ft','maxFI_fc','maxFI_mt','maxFI_mc', ...
        'g_c','c','V','M','iter','dV_win','C0','runtime'};
n = numel(jobs);  D = nan(n, numel(vars));  flag = repmat({'error'}, n, 1);  file = repmat({''}, n, 1);  aggOf = cell(n, 1);
for i = 1:n
    aggOf{i} = aggs{ceil(jobs(i) / nT)};
    th = th_deg(mod(jobs(i) - 1, nT) + 1);
    o = struct('agg',aggOf{i},'projection','heaviside','rmin_phys',4,'ndiv',ndiv,'maxiter',maxiter,'plot',false,'save',true);
    o.theta_init = deg2rad(th);
    o.tag = sprintf('th%03d', th);                         % in the results file name: parallel runs must not overwrite each other
    try
        r  = runCase(o);
        g  = nan(1,4);  g(1:numel(r.g_hs))   = r.g_hs;     % smoothmax has ONE constraint and ONE index: they go in the *_ft columns, the other modes stay NaN
        FI = nan(1,4);  FI(1:numel(r.maxFI)) = r.maxFI;
        D(i,:)  = [th, o.theta_init, g, FI, r.g_c, r.c, r.V, r.M, r.iter, r.dV_win, r.C0, r.runtime];
        flag{i} = r.exitflag;
        file{i} = r.file;
    catch ME
        warning('%s, theta_init = %g deg failed: %s', aggOf{i}, th, ME.message);
        D(i,1:2) = [th, o.theta_init];
    end
    fprintf('%s, theta_init %g deg: %s\n', aggOf{i}, th, flag{i});  % main_heaviside calls clc, so print after each run
end
T = [array2table(D, 'VariableNames', vars), table(flag, file, 'VariableNames', {'exitflag','file'})];
disp(T);
stamp = datestr(now, 'yyyymmdd_HHMMSS');
ua = unique(aggOf, 'stable');
for ia = 1:numel(ua)
    Ta = T(strcmp(aggOf, ua{ia}), :);
    if isTask
        out = fullfile('Results', sprintf('theta_sweep_part_%s_th%03d.csv', ua{ia}, Ta.theta_deg));
    else
        out = fullfile('Results', sprintf('theta_sweep_%s_%s.csv', ua{ia}, stamp));
    end
    writetable(Ta, out);
end