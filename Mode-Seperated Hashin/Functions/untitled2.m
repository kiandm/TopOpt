function demoAggregationOrder()
% Demonstrates why mode-separated and mode-collapsed Hashin can disagree on
% the FEASIBILITY of the exact same design, using just two elements - no
% mesh or FE solve needed. This is a different mechanism from
% demoModeConflict.m (which showed a trade-off AT a single point): this one
% is about the ORDER of aggregation ACROSS different elements.
%
% MECHANISM: mode-separated runs one p-norm PER MODE, across all elements,
% then requires every mode's p-norm to be <=1 independently. Mode-collapsed
% instead combines fibre+matrix at EACH element first (Dong et al. Eq 15-16),
% then runs a SINGLE p-norm across elements on that already-combined value.
%
% Consider two elements that are each exactly as damaged as the other, but
% in DIFFERENT modes: element A is fibre-tension-critical, element B is
% matrix-tension-critical, both at the same fraction of their own limit and
% negligible in every other mode. Mode-separated correctly sees two
% independent constraints, each safely below 1 (the fibre p-norm is
% dominated by A alone - B's contribution there is negligible - and vice
% versa for the matrix p-norm). Mode-collapsed instead produces two elements
% that EACH combine to ~1 (since each is critical in its own mode), then
% p-norms those two near-1 values TOGETHER - and p-norming two
% simultaneously-near-the-limit values pushes the aggregate ABOVE 1, even
% though no single element, in any single mode, actually exceeds its own
% allowable anywhere.

strength.Xt = 1062; strength.Xc = 610; strength.Yt = 31; strength.Yc = 118; strength.S = 72;
Xt = strength.Xt; Yt = strength.Yt;
k_gate = 50; p = 8;
delta_f = 0.05*min(strength.Xt, strength.Xc);
delta_m = 0.05*min(strength.Yt, strength.Yc);
eps_max = 1e-3;

% ---- Headline snapshot: both elements at 95% of their own critical mode ----
frac = 0.95;