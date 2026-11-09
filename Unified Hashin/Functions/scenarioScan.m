% scenarioScan.m  Which failure modes are active at the same time?  Run from this folder.
% For each (load_scale, c_rel) the optimiser is run with the options in 'base'; c_max = c_rel * (compliance of the solid starting design).
% A failure mode is ACTIVE when its shadow price (MMA multiplier, volume fraction saved per unit relaxation of g_k) exceeds tol_lam.
% g_k ~ 0 alone is not enough: a constraint can sit at 0 with a price of 1e-4 while compliance does the work.
% Shadow prices are only meaningful in a CONVERGED run (at maxiter they saturate at hundreds when a constraint is violated).
addpath('Functions');  if ~exist('Results','dir'), mkdir('Results'); end
base      = struct('problem','cantilever','agg','separated','ndiv',50,'rmin_phys',4,'theta_init',0,'theta_fixed',true, ...
                   'maxiter',1000,'beta_step',25,'win',20,'min_iter',0,'plot',false,'save',false);
loads     = [0.35 0.5 0.65 0.8];       % load_scale
c_rels    = [10 30];                   % c_max / compliance of the solid starting design
tol_lam   = 1e-3;
modes     = {'ft','fc','mt','mc'};
rows = {};
for s = loads
    for cr = c_rels
        o = base;  o.load_scale = s;
        o0 = o;  o0.maxiter = 2;  o0.win = 1;  o0.c_max = 1000;   % 2-iteration run: iteration 1 is the solid design
        r0 = runCase(o0);  C0 = r0.C0;                  % compliance of the solid starting design (iteration 1)
        o.c_max = cr * C0;
        r = runCase(o);
        active = modes(r.lam(2:5) > tol_lam);
        fprintf('s=%.2f c_rel=%2d %-9s V=%.3f c/c_max=%.2f  g=[%s]  price=[%s]  active: %s   worst Gauss pt %.2f\n', s, cr, r.exitflag, r.V, (r.g_c + 1), ...
                sprintf('%7.3f', r.gS), sprintf('%8.1e', r.lam), strjoin(active, '+'), r.gp.gpMaxAll);
        rows(end+1,:) = {s, cr, r.exitflag, r.V, (r.g_c + 1), r.gS, r.lam, strjoin(active,'+'), r.gp.gpMaxAll};   %#ok<SAGROW>
    end
end
save(fullfile('Results', ['scenarioScan_' datestr(now,'yyyymmdd_HHMMSS') '.mat']), 'rows', 'base');