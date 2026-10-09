function [TW, g_tw, vonMises] = checkTsaiWu(filename, doPlot)
% Post-hoc Tsai-Wu check on a design that was optimised under a different
% criterion (mode-separated Hashin, via saveOptResults.m). This recomputes
% element stresses from the SAVED displacement field U and design xphy - no
% FE re-solve, no sensitivities - and evaluates the classical 2D Tsai-Wu
% failure index, using the same F1,F2,F11,F22,F66,F12 coefficients and the
% same stress-relaxation exponent q as TsaiWu.m/Hashin.m, so the result is
% directly comparable to the saved Hashin FailIdx. This is exactly the kind
% of cross-criterion comparison Dong et al. 2025 run in their own results
% (e.g. checking a Hashin-optimised design's Tsai-Wu index, and vice versa).
%
% Inputs:
%   filename  path to a .mat file written by saveOptResults.m
%   doPlot    optional, default true - draws a Hashin-vs-Tsai-Wu comparison
%
% Outputs:
%   TW        numele x 1, per-element Tsai-Wu index (>=1 means predicted failure)
%   g_tw      scalar, p-norm aggregate constraint value (Hp - 1); <=0 is safe,
%             same convention as g_hs elsewhere in this project
%   vonMises  numele x 1, per-element von Mises stress, for cross-reference
% 
% Execute this file with "[TW, g_tw, vonMises] = checkTsaiWu(mfilename);" in
% the command window.
addpath('Results\')
if nargin < 2
    doPlot = true;
end
if nargin < 1 || isempty(filename)                 % default: newest result file in this folder
    d = dir('Results/opt_results_*.mat');  names = sort({d.name});  filename = names{end};      % file names carry a timestamp -> last = newest
end
r = load(filename);
xphy = r.xphy; U = r.U; numele = r.numele; gs = r.gs; edofMat = r.edofMat;
matprop = r.matprop; strength = r.strength;
dphix_ref = r.dphix_ref; dphiy_ref = r.dphiy_ref;

q = 0.8; % stress relaxation exponent - matches Hashin.m/TsaiWu.m so this index
         % is on the same footing as the saved Hashin FailIdx it's compared against

% Strength allowables -> Tsai-Wu coefficients (identical to TsaiWu.m)
Xt = strength.Xt; Xc = strength.Xc; Yt = strength.Yt; Yc = strength.Yc; S = strength.S;
F1  = 1/Xt - 1/Xc; F2  = 1/Yt - 1/Yc;
F11 = 1/(Xt*Xc);   F22 = 1/(Yt*Yc);
F66 = 1/S^2;       F12 = -0.5*sqrt(F11*F22);

% Material stiffness in material axes - identical to Hashin.m/TsaiWu.m
E1 = matprop.E1; E2 = matprop.E2; nu12 = matprop.nu12; nu21 = matprop.nu21; G12 = matprop.G12;
C0 = [ E1/(1-nu12*nu21), nu21*E1/(1-nu12*nu21),   0;
       nu12*E2/(1-nu12*nu21), E2/(1-nu12*nu21),    0;
       0,                                    0, G12];

ndof = size(edofMat,2);
TW = zeros(numele,1); vonMises = zeros(numele,1);

% Element loop: recompute stresses from the saved U/xphy and evaluate Tsai-Wu
for e = 1:numele
    Ue = U(edofMat(e,:)); xdens = xphy(e); theta = xphy(numele + e);
    c = cos(theta); s = sin(theta);
    T_eps = [ c^2, s^2,  c*s;              % global strain -> material-axis strain
              s^2, c^2, -c*s;
             -2*c*s, 2*c*s, c^2-s^2 ];
    Tinv = [c^2,  s^2, -2*c*s;             % material-axis stress -> global stress
            s^2,  c^2,  2*c*s;
            c*s, -c*s,  c^2-s^2];

    TW_e = 0; gcount = (e-1)*4;
    for i = 1:2
        for j = 1:2   % 2x2 Gauss quadrature -> 4 integration points per element
            gcount = gcount + 1;
            wt  = gs(6,gcount);
            jac = gs(7,gcount);

            gp_idx = gcount - (e-1)*4;
            dNdx = [dphix_ref(:,gp_idx)'; dphiy_ref(:,gp_idx)'];

            B = zeros(3,ndof);
            B(1,1:2:end) = dNdx(1,:);  B(2,2:2:end) = dNdx(2,:);
            B(3,1:2:end) = dNdx(2,:);  B(3,2:2:end) = dNdx(1,:);

            eps_l = T_eps * (B * Ue);
            sig = xdens^q * (C0 * eps_l);   % density-relaxed material-axis stress
            s1 = sig(1); s2 = sig(2); t12 = sig(3);

            sig_global = Tinv * sig;
            sx = sig_global(1); sy = sig_global(2); txy = sig_global(3);
            vonMises(e) = vonMises(e) + sqrt(sx^2 - sx*sy + sy^2 + 3*txy^2) / 4;

            TW_gp = F1*s1 + F2*s2 + F11*s1^2 + F22*s2^2 + F66*t12^2 + 2*F12*s1*s2;
            TW_e = TW_e + TW_gp * wt * jac;   % integrate over the element area
        end
    end
    TW(e) = TW_e;
end

% p-norm aggregate, same convention as Hashin.m/TsaiWu.m: g <= 0 is safe
p = 8;
TWp  = (sum(TW.^p))^(1/p);
g_tw = TWp - 1;

% Summary, restricted to elements that actually carry material
solid  = xphy(1:numele) >= 0.3;
nSolid = sum(solid);
nFail  = sum(TW >= 1 & solid);
fprintf('Tsai-Wu check on %s:\n', filename);
fprintf('  max(TW) over solid elements = %.3f, g_tw (p-norm - 1) = %.3f\n', ...
    max(TW(solid)), g_tw);
fprintf('  %d / %d solid elements exceed TW = 1 (%.1f%%)\n', nFail, nSolid, 100*nFail/max(nSolid,1));
if nFail == 0
    fprintf('  -> Design satisfies Tsai-Wu everywhere it satisfies Hashin.\n');
else
    fprintf('  -> Design VIOLATES Tsai-Wu in %d elements despite being Hashin-optimised.\n', nFail);
end

if doPlot
    if isfield(r, 'FailIdx'), FI = r.FailIdx; else, FI = r.TW; end   % mode-separated or collapsed results
    mask = ~solid;
    hs_f = max(FI, [], 2);   hs_f(mask) = NaN;                       % worst Hashin sub-mode, per element
    tw_f = TW;               tw_f(mask) = NaN;
    d_f  = hs_f - tw_f;                                              % negative = Tsai-Wu is the more critical criterion
    fields = {hs_f, tw_f, d_f};
    ttl    = {'max Hashin index (as optimised)', 'Tsai-Wu index (post-hoc)', 'Hashin - Tsai-Wu'};
    figure; set(gcf, 'Color', 'white');
    for k = 1:3
        subplot(1,3,k);
        patch('Faces', r.conn', 'Vertices', r.coords', 'FaceVertexCData', fields{k}, ...
              'FaceColor', 'flat', 'EdgeColor', 'none');
        axis equal off; colorbar; title(ttl{k});
        if k < 3
            clim([0 1.2]);
        else                                                         % diverging scale centred on 0
            lim = max(abs(d_f));  clim([-lim lim]);
            colormap(gca, [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ones(128,1) linspace(1,0,128)' linspace(1,0,128)']);
        end
    end
    % fprintf('  max(Hashin - TW) over solid elements = %.3f\n', max(d_f));
    fprintf('  max(Hashin - TW) over solid elements = %.3f\n', max(d_f, [], 'all'));
end
end