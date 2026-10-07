function demoModeConflict()
% Engineers a single material-point stress state where fibre-tension and
% matrix-tension are in genuine conflict for the fibre-angle design variable,
% and compares how this project's mode-separated Hashin treats it against
% Dong et al. 2025's mode-collapsed (smooth-max) formulation. No mesh or FE
% solve needed - this operates directly on stresses, the way Hashin.m does
% internally at a single Gauss point.
%
% MECHANISM: for a fixed global stress state (sx,sy,txy), the material-axis
% stresses after rotating by fibre angle theta satisfy s1(theta)+s2(theta) =
% sx+sy for every theta - the trace of the stress tensor is invariant under
% rotation. So wherever both s1 and s2 are tensile, rotating theta can only
% trade fibre stress against matrix stress; it can never reduce both at once.
%
% Because this material's Yt is ~34x smaller than Xt (31 vs 1062 MPa), a
% stress split that is EQUALLY damaging to both modes puts nearly all the
% raw stress into s1 and very little into s2 - and moving away from that
% split degrades the matrix index far faster than it improves the fibre
% index. That produces a sharp corner in the true (hard-max) envelope right
% where I_ft = I_mt, which is the point engineered below.

strength.Xt = 1062; strength.Xc = 610; strength.Yt = 31; strength.Yc = 118; strength.S = 72;
Xt = strength.Xt; Yt = strength.Yt; S = strength.S;

% Engineer sx,sy (at theta=0, material axes = global axes) so I_ft = I_mt exactly.
ratio   = Xt/Yt;                 % how much more stress the fibre direction tolerates
S_total = 550;                   % sx+sy, MPa: the rotation-invariant stress "budget"
s2_0 = S_total/(ratio+1);
s1_0 = S_total - s2_0;
sx = s1_0; sy = s2_0; txy = 0;

fprintf('sx=%.1f, sy=%.1f MPa (sx+sy=%.1f, invariant under fibre rotation)\n', sx, sy, sx+sy);
fprintf('At theta=0: s1=%.1f, s2=%.1f -> I_ft=%.4f, I_mt=%.4f (tied by construction)\n\n', ...
    s1_0, s2_0, (s1_0/Xt)^2, (s2_0/Yt)^2);

k_gate  = 50;
delta_f = 0.05*min(strength.Xt, strength.Xc);
delta_m = 0.05*min(strength.Yt, strength.Yc);
eps_max = 1e-3;

theta_deg = -30:0.5:30;
theta = deg2rad(theta_deg);
c = cos(theta); s = sin(theta);
s1  = c.^2*sx + s.^2*sy + 2*c.*s*txy;
s2  = s.^2*sx + c.^2*sy - 2*c.*s*txy;
t12 = -c.*s*sx + c.*s*sy + (c.^2-s.^2)*txy;

% Sigmoid gates (Hashin.m convention) - confirm they stay saturated at ~1
% across this sweep, so fibre-compression/matrix-compression stay ~0 and
% the I_ft/I_mt comparison below is the whole story.
g1 = 1./(1+exp(-k_gate*s1));
g2 = 1./(1+exp(-k_gate*s2));
fprintf('Gate range across +-30 deg: g1 in [%.6f, %.6f], g2 in [%.6f, %.6f]\n', ...
    min(g1), max(g1), min(g2), max(g2));
fprintf('(both stay ~1 throughout, so fibre-compression/matrix-compression stay ~0)\n\n');

% ---- mode-separated (this project's Hashin.m) ----
I_ft = (s1/Xt).^2 + (t12/S).^2;
I_mt = (s2/Yt).^2 + (t12/S).^2;

% ---- mode-collapsed (Dong et al. 2025, Eq 12-16, revived from Hashin.asv) ----
H_collapsed = zeros(size(theta));
for k = 1:numel(theta)
    Hf = heavisideBlend(strength.Xt, strength.Xc, s1(k), delta_f);
    Hm = heavisideBlend(strength.Yt, strength.Yc, s2(k), delta_m);
    FI_f = sqrt((s1(k)/Hf)^2 + (t12(k)/S)^2);
    FI_m = sqrt((s2(k)/Hm)^2 + (t12(k)/S)^2);
    diff_fm = FI_f - FI_m;
    H_collapsed(k) = 0.5*((FI_f+FI_m) + sqrt(diff_fm^2+eps_max) - sqrt(eps_max));
end
H_collapsed_sq = H_collapsed.^2;   % squared, for a fair comparison against I_ft/I_mt

figure;
plot(theta_deg, I_ft, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_deg, I_mt, 'r-', 'LineWidth', 1.5);
plot(theta_deg, max(I_ft, I_mt), 'k--', 'LineWidth', 1);
plot(theta_deg, H_collapsed_sq, 'g-.', 'LineWidth', 1.5);
yline(1, 'k:');
xlabel('Fibre angle \theta (deg, rotated away from the engineered tie point)');
ylabel('Failure index');
legend('I_{ft} (separated)', 'I_{mt} (separated)', 'max(I_{ft},I_{mt}) (true, hard)', ...
       'H_{Hs} (Dong et al. smooth-max)', 'failure = 1', 'Location', 'best');
title('Fibre/matrix tension conflict at a single material point');
grid on;
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