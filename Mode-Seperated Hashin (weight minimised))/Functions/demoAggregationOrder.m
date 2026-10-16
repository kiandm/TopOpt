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

% Crosses threshold at 97.07% allowable stress fraction with mode-collapsed (index squared, like-for-like with the
% separated I = (s/X)^2 + (t/S)^2) vs 99.76% mode-separated.  (With the un-squared collapsed index it is 92.99%.)

strength.Xt = 1062; strength.Xc = 610; strength.Yt = 31; strength.Yc = 118; strength.S = 72;
Xt = strength.Xt; Yt = strength.Yt;
k_gate = 50; p = 8;
delta_f = 0.05*min(strength.Xt, strength.Xc);
delta_m = 0.05*min(strength.Yt, strength.Yc);
eps_max = 1e-3;

% ---- Headline snapshot: both elements at 98% of their own critical mode ----
frac = 0.98;
sA = [frac*Xt, 2.0, 5.0];    % element A: fibre-tension critical, safe elsewhere, material-axis stress [s1, s2, t12], s1 is close to Xt
sB = [10.0, frac*Yt, 5.0];   % element B: matrix-tension critical, safe elsewhere, material-axis stress [s1, s2, t12], s2 is close to Yt
% Both stresses are at 98% allowable stress

fprintf('Element A (fibre-critical):  s1=%.1f, s2=%.1f, t12=%.1f\n', sA(1), sA(2), sA(3));
fprintf('Element B (matrix-critical): s1=%.1f, s2=%.1f, t12=%.1f\n\n', sB(1), sB(2), sB(3));

[IftA, ImtA] = separatedFtMt(sA, strength, k_gate);
[IftB, ImtB] = separatedFtMt(sB, strength, k_gate);
g_ft = (IftA^p + IftB^p)^(1/p) - 1;
g_mt = (ImtA^p + ImtB^p)^(1/p) - 1;

HcA = collapsedIndex(sA, strength, delta_f, delta_m, eps_max);
HcB = collapsedIndex(sB, strength, delta_f, delta_m, eps_max);
g_h = (HcA^(2*p) + HcB^(2*p))^(1/p) - 1;

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
fracs = 0.5:0.001:1.0;
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
    g_h_sweep(k) = (HcA_k^(2*p) + HcB_k^(2*p))^(1/p) - 1;
end
% ---- Where each curve crosses g = 0 (linear interpolation between sweep points) ----
x_sep = zeroCross(fracs*100, g_sep_sweep);
x_col = zeroCross(fracs*100, g_h_sweep);
fprintf('---- g = 0 crossings (stress as %% of each element''s own critical allowable) ----\n');
fprintf('  mode-separated (worst of g_ft, g_mt): x = %.2f %%\n', x_sep);
fprintf('  mode-collapsed (g_h, squared index):  x = %.2f %%\n', x_col);
if isnan(x_sep) || isnan(x_col)
    fprintf('  No crossing inside the swept range (%.0f-%.0f %%).\n\n', fracs(1)*100, fracs(end)*100);
else
    fprintf('  Separated approach: +%.2f percentage points, i.e. %.2f %% improvement over collapsed\n', x_sep - x_col, 100*(x_sep - x_col)/x_col);
    fprintf('  (stress scales linearly with load, so this is also the extra load the design can carry)\n\n');
end

% ---- Does the gap survive when both formulations use the SAME exponent on the stress-ratio index? ----
%   mode-separated  ||I||_{P/2} <= 1  with I = FI^2   (Hashin.m's p = 8 on I is P = 16)
%   mode-collapsed  ||Hc||_P    <= 1  with Hc linear  (Dong et al. as published: p = 8 is P = 8)
Ps = [4 8 16 32 64];
x_sep_P = zeros(size(Ps));   x_col_P = zeros(size(Ps));
for kk = 1:numel(Ps)
    x_sep_P(kk) = crossingPct('separated', Ps(kk), strength, k_gate, delta_f, delta_m, eps_max);
    x_col_P(kk) = crossingPct('collapsed', Ps(kk), strength, k_gate, delta_f, delta_m, eps_max);
end
fprintf('---- g = 0 crossing (stress as %% of allowable) against exponent P on the stress-ratio index ----\n');
fprintf('     P    separated   collapsed   gap [pts]   improvement [%%]\n');
for kk = 1:numel(Ps)
    fprintf('  %4d    %8.2f    %8.2f    %7.2f    %8.2f\n', Ps(kk), x_sep_P(kk), x_col_P(kk), ...
            x_sep_P(kk)-x_col_P(kk), 100*(x_sep_P(kk)-x_col_P(kk))/x_col_P(kk));
end
fprintf(['  Hashin.m with p = 8 on I corresponds to P = 16; Dong''s linear index with p = 8 corresponds to P = 8.\n' ...
         '  As published, the comparison is therefore the P = 16 separated row against the P = 8 collapsed row.\n\n']);
figure;
semilogx(Ps, x_sep_P, 'b-o', 'LineWidth', 1.5); hold on;
semilogx(Ps, x_col_P, 'g-s', 'LineWidth', 1.5);
yline(100, 'k:');
xlabel('Exponent P applied to the stress-ratio index'); ylabel('Stress at which g = 0 [% of allowable]');
legend('separated', 'collapsed', '100 %', 'Location', 'best'); grid on;
title('Feasibility threshold against p-norm exponent');

figure;
plot(fracs*100, g_sep_sweep, 'b-', 'LineWidth', 1.5); hold on;
plot(fracs*100, g_h_sweep, 'g-.', 'LineWidth', 1.5);
yline(0, 'k:');
xlabel('Each element''s stress, as a % of its own critical mode''s allowable');
ylabel('Constraint value g (<=0 is feasible)');
legend('worst of g_{ft}, g_{mt} (separated)', 'g_h (collapsed, squared index)', 'feasibility boundary', 'Location', 'best');
title('Constraint value against percentage of critical stress, showing difference in feasibility');
grid on;

% Printing when each line crosses y=0, and the percentage change between
% this crossing for g_sep and g_h


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

function x0 = zeroCross(x, y)
% first upward crossing of y through 0, linearly interpolated; NaN if none
i = find(y(1:end-1) <= 0 & y(2:end) > 0, 1);
if isempty(i)
    x0 = NaN;
else
    x0 = x(i) - y(i)*(x(i+1) - x(i))/(y(i+1) - y(i));
end
end


function x = crossingPct(formulation, P, strength, k_gate, delta_f, delta_m, eps_max)
% stress (percent of each element's own critical allowable) at which the two-element constraint reaches g = 0
x = 100 * fzero(@(fr) gAtFrac(fr, formulation, P, strength, k_gate, delta_f, delta_m, eps_max), [0.3 1.2]);
end

function g = gAtFrac(fr, formulation, P, strength, k_gate, delta_f, delta_m, eps_max)
sA = [fr*strength.Xt, 2.0, 5.0];     % element A: fibre-tension critical
sB = [10.0, fr*strength.Yt, 5.0];    % element B: matrix-tension critical
switch formulation
    case 'separated'                 % ||I||_{P/2} per mode, I = FI^2, binding mode
        [IA_f, IA_m] = separatedFtMt(sA, strength, k_gate);
        [IB_f, IB_m] = separatedFtMt(sB, strength, k_gate);
        q = P/2;
        g = max((IA_f^q + IB_f^q)^(1/q), (IA_m^q + IB_m^q)^(1/q)) - 1;
    case 'collapsed'                 % ||Hc||_P with the linear (un-squared) index
        HA = collapsedIndex(sA, strength, delta_f, delta_m, eps_max);
        HB = collapsedIndex(sB, strength, delta_f, delta_m, eps_max);
        g = (HA^P + HB^P)^(1/P) - 1;
end
end

function s = passFail(g)
if g <= 0
    s = 'safe';
else
    s = 'VIOLATED';
end
end
