function opts = setOpts(user)
% SETOPTS  Single source of truth for every option of a run.
%   opts = setOpts()                       all defaults
%   opts = setOpts(struct('p',16))         defaults with p overridden
% An unknown field name is an ERROR (a typo must never silently run the default),
% and so is an option value that is not implemented yet.

if nargin < 1 || isempty(user), user = struct(); end

% ---- problem -------------------------------------------------------------------------
d.problem     = 'Lbrac60';   % 'Lbrac60' | 'mbb' | 'cantilever'  (Functions/problem_setup_<problem>.m)
d.ndiv        = 100;         % elements along the 100 mm side (Lbrac60 only)
d.rmin_phys   = 4;           % density / angle filter radius [mm]
d.c_max       = 1941.74;     % compliance limit  F'*U
% ---- material model ------------------------------------------------------------------
d.penal       = 3.0;         % SIMP penalty
d.q           = 0.8;         % stress relaxation exponent
% ---- failure constraint --------------------------------------------------------------
d.agg         = 'separated'; % 'separated': 4 Hashin modes, one p-norm each   |   'smoothmax': Dong-type, one p-norm
d.p           = 8;           % p-norm exponent P, applied to the (linear) stress-ratio index
d.gate        = 'sigmoid';   % separated only. 'sigmoid' (a cubic gate is added in Phase 2)
d.k_gate      = 50;          % sigmoid sharpness [1/MPa]
d.eps_max     = 1e-3;        % smoothmax only: smooth-max regularisation
d.eps_sq      = 1e-12;       % keeps sqrt(I) differentiable at I = 0
d.delta_frac  = 0.05;        % smoothmax only: cubic Heaviside half-width = delta_frac*min(Xt,Xc) (fibre), min(Yt,Yc) (matrix)
d.pcorr       = 'off';       % adaptive p-norm correction. 'off' (added in Phase 2)
% ---- filters and projection ----------------------------------------------------------
d.thetafilter = 'double';    % 'double' = pi-periodic double-angle filter ('continuity' is added in Phase 2)
d.eta         = 0.5;         % Heaviside threshold
d.beta0       = 1;           % initial Heaviside sharpness
d.beta_max    = 32;          % final Heaviside sharpness
d.beta_factor = 1.5;         % beta is multiplied by this ...
d.beta_step   = 50;          % ... every beta_step iterations
% ---- optimiser -----------------------------------------------------------------------
d.maxiter     = 1000;
d.crow        = [];          % compliance-row scale for MMA. [] = automatic: norm(df0dx)/norm(dc/c_max) at iteration 1
d.fscale_coef = 0.02;        % objective scale: fscale = fscale_coef*numele
d.theta_init  = pi/2;        % initial fibre angle [rad]
d.theta_pert  = 0;           % [rad] seeded random perturbation of the initial angle, +-theta_pert/2 (needs a seed)
d.seed        = [];          % [] = none
% ---- stopping rule: volume plateau after beta_max is reached, with the constraints satisfied -----
d.win         = 20;          % plateau window [iterations]
d.tol_plateau = 1e-3;        % relative volume change over the window
d.tol_g       = 1e-2;        % largest constraint value allowed at exit
% ---- output --------------------------------------------------------------------------
d.plot        = true;
d.save        = true;

% ---- merge the user's values -----------------------------------------------------------
names = fieldnames(user);
for k = 1:numel(names)
    if ~isfield(d, names{k})
        error('setOpts:unknownOption', 'Unknown option "%s".  Valid options:\n  %s', names{k}, strjoin(fieldnames(d)', ', '));
    end
    d.(names{k}) = user.(names{k});
end
opts = d;

% ---- validation --------------------------------------------------------------------------
oneOf(opts, 'problem',     {'Lbrac60', 'mbb', 'cantilever'});
oneOf(opts, 'agg',         {'separated', 'smoothmax'});
oneOf(opts, 'gate',        {'sigmoid'}, '  (the cubic gate is added in Phase 2)');
oneOf(opts, 'pcorr',       {'off'},     '  (the adaptive correction is added in Phase 2)');
oneOf(opts, 'thetafilter', {'double'},  '  (the continuity filter is added in Phase 2)');
if opts.p < 1,                          error('setOpts:p', 'opts.p must be >= 1.'); end
if opts.beta_max < opts.beta0,          error('setOpts:beta', 'opts.beta_max must be >= opts.beta0.'); end
if opts.maxiter < 1 || opts.maxiter ~= round(opts.maxiter), error('setOpts:maxiter', 'opts.maxiter must be a positive integer.'); end
if opts.win >= opts.maxiter,            error('setOpts:win', 'opts.win must be smaller than opts.maxiter.'); end
if ~isempty(opts.crow) && opts.crow <= 0, error('setOpts:crow', 'opts.crow must be positive, or [] for automatic.'); end
if opts.theta_pert < 0 || opts.theta_pert > pi, error('setOpts:theta_pert', 'opts.theta_pert must lie in [0, pi].'); end
if opts.theta_pert > 0 && isempty(opts.seed), error('setOpts:seed', 'opts.theta_pert > 0 needs opts.seed, so that the run can be repeated.'); end
if ~strcmp(opts.problem, 'Lbrac60') && opts.ndiv ~= 100
    error('setOpts:ndiv', 'opts.ndiv is only supported for problem "Lbrac60".');
end
end

function oneOf(opts, field, allowed, note)
if nargin < 4, note = ''; end
if ~any(strcmp(opts.(field), allowed))
    v = opts.(field); if ~ischar(v), v = '<not a string>'; end
    error('setOpts:value', 'opts.%s = "%s" is not available.  Allowed: %s%s', field, v, strjoin(allowed, ', '), note);
end
end