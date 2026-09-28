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

% Crosses threshold at 92.993% allowable stress fraction with
% mode-collapsed vs 99.756% mode-separated

strength.Xt = 1062; strength.Xc = 610; strength.Yt = 31; strength.Yc = 118; strength.S = 72;
Xt = strength.Xt; Yt = strength.Yt;
k_gate = 50; p = 8;
delta_f = 0.05*min(strength.Xt, strength.Xc);
delta_m = 0.05*min(strength.Yt, strength.Yc);
eps_max = 1e-3;

% ---- Headline snapshot: both elements at 95% of their own critical mode ----
frac = 0.95;
sA = [frac*Xt, 2.0, 5.0];    % element A: fibre-tension critical, safe elsewhere, material-axis stress [s1, s2, t12], s1 is close to Xt
sB = [10.0, frac*Yt, 5.0];   % element B: matrix-tension critical, safe elsewhere, material-axis stress [s1, s2, t12], s2 is close to Yt
% Both stresses are at 95% allowable stress

fprintf('Element A (fibre-critical):  s1=%.1f, s2=%.1f, t12=%.1f\n', sA(1), sA(2), sA(3));
fprintf('Element B (matrix-critical): s1=%.1f, s2=%.1f, t12=%.1f\n\n', sB(1), sB(2), sB(3));

[IftA, ImtA] = separatedFtMt(sA, strength, k_gate);
[IftB, ImtB] = separatedFtMt(sB, strength, k_gate);
g_ft = (IftA^p + IftB^p)^(1/p) - 1;
g_mt = (ImtA^p + ImtB^p)^(1/p) - 1;

HcA = collapsedIndex(sA, strength, delta_f, delta_m, eps_max);
HcB = collapsedIndex(sB, strength, delta_f, delta_m, eps_max);
g_h = (HcA^p + HcB^p)^(1/p) - 1;

fprintf('---- mode-separated (this project''s Hashin.m) ----\n');
fprintf('  g_ft = %+.4f  (%s)\n', g_ft, passFail(g_ft));
fprintf('  g_mt = %+.4f  (%s)\n\n', g_mt, passFail(g_mt));
fprintf('---- mode-collapsed (Dong et al. 2025, single combined constraint) ----\n');
fprintf('  g_h  = %+.4f  (%s)\n\n', g_h, passFail(g_h));
if g_h > 0 && g_ft <= 0 && g_mt <= 0
    fprintf(['Neither element individually breaches ANY mode, yet the collapsed\n' ...
             'formulation reports the design as infeasible. Separating the modes\n' ...
             'avoids this p-norm "double counting" across elements that are each\n' ...
             'critical in a different mode.\n\n']);
end

% ---- Sweep: how far apart are the two formulations across criticality? ----
fracs = 0:0.02:1.0;
g_sep_sweep = zeros(size(fracs));
g_h_sweep   = zeros(size(fracs));
for k = 1:numel(fracs)
    sA_k = [fracs(k)*Xt, 2.0, 5.0];
    sB_k = [10.0, fracs(k)*Yt, 5.0];

    [IftA_k, ImtA_k] = separatedFtMt(sA_k, strength, k_gate);
    [IftB_k, ImtB_k] = separatedFtMt(sB_k, strength, k_gate);
    g_ft_k = (IftA_k^p + IftB_k^p)^(1/p) - 1;
    g_mt_k = (ImtA_k^p + ImtB_k^p)^(1/p) - 1;
    g_sep_sweep(k) = max(g_ft_k, g_mt_k);   % binding (worst) separated constraint

    HcA_k = collapsedIndex(sA_k, strength, delta_f, delta_m, eps_max);
    HcB_k = collapsedIndex(sB_k, strength, delta_f, delta_m, eps_max);
    g_h_sweep(k) = (HcA_k^p + HcB_k^p)^(1/p) - 1;
end

figure;
plot(fracs*100, g_sep_sweep, 'b-', 'LineWidth', 1.5); hold on;
plot(fracs*100, g_h_sweep, 'g-.', 'LineWidth', 1.5);
yline(0, 'k:');
xlabel('Each element''s stress, as a % of its own critical mode''s allowable');
ylabel('Constraint value g (<=0 is feasible)');
legend('worst of g_{ft}, g_{mt} (separated)', 'g_h (collapsed)', 'feasibility boundary', 'Location', 'best');
title('Constraint value against percentage of critical stress, showing difference in feasibility');
grid on;
end

function [I_ft_gated, I_mt_gated] = separatedFtMt(sig, strength, k_gate)
s1 = sig(1); s2 = sig(2); t12 = sig(3);
Xt = strength.Xt; Yt = strength.Yt; S = strength.S;
g1 = 1/(1+exp(-k_gate*s1));
g2 = 1/(1+exp(-k_gate*s2));
I_ft_gated = g1 * ((s1/Xt)^2 + (t12/S)^2);
I_mt_gated = g2 * ((s2/Yt)^2 + (t12/S)^2);
end

function H_gp = collapsedIndex(sig, strength, delta_f, delta_m, eps_max)
s1 = sig(1); s2 = sig(2); t12 = sig(3);
S = strength.S;
Hf = heavisideBlend(strength.Xt, strength.Xc, s1, delta_f);
Hm = heavisideBlend(strength.Yt, strength.Yc, s2, delta_m);
FI_f = sqrt((s1/Hf)^2 + (t12/S)^2);
FI_m = sqrt((s2/Hm)^2 + (t12/S)^2);
diff_fm = FI_f - FI_m;
H_gp = 0.5*((FI_f+FI_m) + sqrt(diff_fm^2+eps_max) - sqrt(eps_max));
end

function h = heavisideBlend(a, b, c, delta)
% Eq. 12 of Dong et al. 2025 - see Hashin.asv for the same function.
if c > delta
    h = a;
elseif c < -delta
    h = b;
else
    h = 0.75*(a-b)*(c/delta - c^3/(3*delta^3)) + (a+b)/2;
end
end

function s = passFail(g)
if g <= 0
    s = 'safe';
else
    s = 'VIOLATED';
end
end
