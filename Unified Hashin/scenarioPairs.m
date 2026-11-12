% scenarioPairs.m  Separated vs smooth-max on the four-mode scenario: cantilever, fibres fixed along the beam axis, compliance limit
% c_max = c_rel * (compliance of the solid start) so that the failure modes, not stiffness, govern the design.
% The fibre field is fixed, so slightly different loads act as replicate instances. Same start and options for both aggregations.
%   serial :  scenarioPairs                         every (aggregation, load) run one after another, then the summary table
%   Slurm  :  array task k runs ONE run (tasks 1..n = aggs{1}, n+1..2n = aggs{2}); afterwards run  summarisePairs(ndiv)
addpath('Functions');  if ~exist('Results','dir'), mkdir('Results'); end
%% settings (edit here)
aggs  = {'separated', 'smoothmax'};
loads = 0.76:0.01:0.90;                  % load_scale instances (0.74 - 0.86 was the feasible window in my ndiv = 50 runs)
c_rel = 10;                              % c_max = c_rel * compliance of the solid starting design
ndiv  = 50;                              % 50 reproduces my numbers (2 mm load patch); then repeat with 100
base  = struct('problem','cantilever','ndiv',ndiv,'rmin_phys',4,'theta_init',0,'theta_fixed',true, ...
               'maxiter',1500,'beta_step',25,'win',20,'min_iter',0,'plot',false,'save',false);
%% which runs
nL = numel(loads);  nJobs = numel(aggs) * nL;
task = str2double(getenv('SLURM_ARRAY_TASK_ID'));       % NaN when not running as a Slurm array
if isnan(task), jobs = 1:nJobs; else, jobs = task; end
if any(jobs < 1 | jobs > nJobs | jobs ~= round(jobs)), error('scenarioPairs:task', 'task %g, but there are %d runs: use --array=1-%d', task, nJobs, nJobs); end
for j = jobs
    ia = ceil(j / nL);  il = mod(j - 1, nL) + 1;
    o = base;  o.agg = aggs{ia};  o.load_scale = loads(il);
    o0 = o;  o0.maxiter = 2;  o0.win = 1;  o0.c_max = 1000;      % 2-iteration run: iteration 1 is the solid design
    r0 = runCase(o0);  o.c_max = c_rel * r0.C0;
    r = runCase(o);
    r.load = loads(il);  r.c_rel = c_rel;
    save(fullfile('Results', sprintf('scenarioPairs_%s_s%.2f_ndiv%d.mat', aggs{ia}, loads(il), ndiv)), 'r');
    fprintf('%-10s load %.2f: %-9s V = %.4f  (%d iterations)\n', aggs{ia}, loads(il), r.exitflag, r.V, r.iter);
end
if isnan(task), summarisePairs(ndiv, aggs); end