function [g_h, dgh_dx, dgh_dtheta, FailIdx, HashinIdx_gp, vonMises] = Hashin( ...
    U, dK, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs, dphix_ref, dphiy_ref)
% Classical (non-mode-separated) Hashin stress constraint with adjoint sensitivities
% Each of the four classical Hashin sub-modes - fibre tension (ft), fibre
% compression (fc), matrix tension (mt), matrix compression (mc) - is
% aggregated independently over the mesh via a p-norm, giving four
% separate constraints g = [g_ft; g_fc; g_mt; g_mc] for the optimiser
% instead of one combined index.

% Strength allowables
Xt = strength.Xt; Xc = strength.Xc;
Yt = strength.Yt; Yc = strength.Yc;
S  = strength.S;   % NOTE: in-plane shear strength reused as transverse (S23) allowable
                   % for the matrix-compression term - standard simplification

% Smoothing parameters (Dong et al. 2025, Eqs. 12 & 16)
delta_f = 0.05*min(Xt,Xc);  % Heaviside transition half-width, fibre mode (tune)
delta_m = 0.05*min(Yt,Yc);  % Heaviside transition half-width, matrix mode (tune)
eps_max = 1e-3;             % smooth-max regularization (tune)

% Material stiffness (material coordinates)
E1   = matprop.E1; E2   = matprop.E2; nu12 = matprop.nu12;
nu21 = matprop.nu21; G12  = matprop.G12;
C0 = [ E1/(1-nu12*nu21), nu21*E1/(1-nu12*nu21),   0;
       nu12*E2/(1-nu12*nu21), E2/(1-nu12*nu21),    0;
       0,                                    0, G12];

% Initialisation
HashinIdx = zeros(numele,1); dHdx    = zeros(numele,1);
dHdth     = zeros(numele,1); vonMises = zeros(numele,1);
ndof = size(edofMat,2);
dKE0th = zeros(ndof,ndof,numele);
fadj_elem = zeros(numele,ndof);
q = 0.8; % stress interpolation exponent 

% Element loop
for e = 1:numele
    Ue = U(edofMat(e,:)); xdens = xphy(e); theta = xphy(numele + e);
    c = cos(theta); s = sin(theta);
    T_eps = [ c^2, s^2,  c*s;
              s^2, c^2, -c*s;
             -2*c*s, 2*c*s, c^2-s^2 ];
    Tinv = [c^2,  s^2, -2*c*s;
            s^2,  c^2,  2*c*s;
            c*s, -c*s,  c^2-s^2];
    dTinv_dth = [-sin(2*theta),  sin(2*theta), -2*cos(2*theta);
                  sin(2*theta), -sin(2*theta),  2*cos(2*theta);
                  cos(2*theta), -cos(2*theta), -2*sin(2*theta)];
    dCxy_dth = dTinv_dth * C0 * Tinv' + Tinv * C0 * dTinv_dth'; 
    dT_eps_dth = [-sin(2*theta),  sin(2*theta),   cos(2*theta);
                           sin(2*theta), -sin(2*theta),  -cos(2*theta);
                          -2*cos(2*theta), 2*cos(2*theta), -2*sin(2*theta)];
    H_e   = 0.0; dHdx_e = 0.0;
    fadj_e = zeros(ndof,1); dH_dth_e = 0.0;
    dKE_dth   = zeros(8,8); gcount = (e-1)*4;
    for i = 1:2
        for j = 1:2
            gcount = gcount + 1;
            wt  = gs(6,gcount);
            jac = gs(7,gcount);
            gp_idx = gcount - (e-1)*4; % local Gauss point index within this element (1..4)
            dNdx = [dphix_ref(:,gp_idx)'; dphiy_ref(:,gp_idx)'];
            % B-matrix
            B = zeros(3,ndof);
            B(1,1:2:end) = dNdx(1,:);  B(2,2:2:end) = dNdx(2,:);
            B(3,1:2:end) = dNdx(2,:);  B(3,2:2:end) = dNdx(1,:);
            dKE_dth = dKE_dth + (xdens^penal) * jac * wt * B' * dCxy_dth * B;
            % Strain/stress (material axes)
            eps_l = T_eps * (B * Ue);
            sig_unscaled = C0 * eps_l;
            sig = xdens^q * sig_unscaled;
            s1 = sig(1); s2 = sig(2); t12 = sig(3);
            sig_global = Tinv * sig;
            sx = sig_global(1); sy = sig_global(2); txy = sig_global(3);
            vm_gp = sqrt(sx^2 - sx*sy + sy^2 + 3*txy^2);
            vonMises(e) = vonMises(e) + vm_gp / 4;
            %% ---- Hashin failure index at this Gauss point (Dong et al. 2025, Eqs. 12-16) ----
            % Hf = heavisideBlend(Xt, Xc, s1, delta_f);   % Eq. 12, blends the ALLOWABLE, not the index
            % Hm = heavisideBlend(Yt, Yc, s2, delta_m);     
            s1c = min(max(s1, -delta_f), delta_f);
            Hf = 0.75*(Xt-Xc)*(s1c/delta_f - s1c^3/(3*delta_f^3)) + (Xt+Xc)/2;
            s2c = min(max(s2, -delta_m), delta_m);
            Hm = 0.75*(Yt-Yc)*(s2c/delta_m - s2c^3/(3*delta_m^3)) + (Yt+Yc)/2;
            
            FI_f = sqrt((s1/Hf)^2 + (t12/S)^2);         % Eq. 13
            FI_m = sqrt((s2/Hm)^2 + (t12/S)^2);         % Eq. 14           
            diff_fm = FI_f - FI_m;
            root_fm = sqrt(diff_fm^2 + eps_max);
            H_gp = 0.5*((FI_f+FI_m) + root_fm - sqrt(eps_max));   % Eq. 16 (corrected)          
            H_e = H_e + H_gp * wt * jac;
            
            %% ---- derivatives of H_gp w.r.t. (s1, s2, t12) ----
            %dHf_ds1 = dHeavisideBlend(Xt, Xc, s1, delta_f);
            %dHm_ds2 = dHeavisideBlend(Yt, Yc, s2, delta_m); 
            dHf_ds1 = 0.75*(Xt-Xc)*(1/delta_f - s1c^2/delta_f^3);
            dHm_ds2 = 0.75*(Yt-Yc)*(1/delta_m - s2c^2/delta_m^3);
            dFIf_ds1  = (1/FI_f) * (s1/Hf) * (1/Hf - s1/Hf^2*dHf_ds1);   % Eq. 34
            dFIf_dt12 = (1/FI_f) * (t12/S^2);
            dFIm_ds2  = (1/FI_m) * (s2/Hm) * (1/Hm - s2/Hm^2*dHm_ds2);   % Eq. 35
            dFIm_dt12 = (1/FI_m) * (t12/S^2);           
            dmax_dFIf = 0.5 + 0.5*diff_fm/root_fm;    % Eq. 32
            dmax_dFIm = 0.5 - 0.5*diff_fm/root_fm;    % Eq. 33            
            psi = [ dmax_dFIf * dFIf_ds1;
                    dmax_dFIm * dFIm_ds2;
                    dmax_dFIf * dFIf_dt12 + dmax_dFIm * dFIm_dt12 ];
            %% sensitivities
            % Density
            dsig_dx = q * xdens^(q-1) * sig_unscaled;
            dHdx_gp = psi' * dsig_dx;
            dHdx_e  = dHdx_e + dHdx_gp * wt * jac;
            dH_dth_gp = psi' * (xdens^q * C0 * dT_eps_dth * (B * Ue));
            dH_dth_e  = dH_dth_e + dH_dth_gp * wt * jac;
            % Adjoint RHS contribution
            fadj_gp = B' * T_eps' * C0' * psi * xdens^q;
            fadj_e = fadj_e + fadj_gp * wt * jac;
        end
    end
    HashinIdx(e) = H_e;
    dHdx(e)      = dHdx_e;
    dHdth(e)     = dH_dth_e;
    % fadj_elem{e} = fadj_e;
    fadj_elem(e,:) = fadj_e';
    % dKE0th{e}    = dKE_dth;
    dKE0th(:,:,e)  = dKE_dth;
end
% p-norm aggregation over elements 
p = 8; % 8 16 32
Hp  = (sum(HashinIdx.^p))^(1/p);
g_h = Hp - 1;
fac = (HashinIdx.^(p-1)) / (Hp^(p-1));
HashinIdx_gp = HashinIdx; 
% assemble adjoint RHS
fadj = accumarray(edofMat(:), reshape(fac .* fadj_elem, [], 1), [size(U,1) 1]);
% adjoint solve
lambda = zeros(size(U));
lambda(freedofs) = dK \ fadj(freedofs);
% final sensitivities
Uall = U(edofMat)';       % ndof x numele
Lall = lambda(edofMat)';  % ndof x numele
KU  = reshape(pagemtimes(KE0,    reshape(Uall,ndof,1,numele)), ndof, numele);
dKU = reshape(pagemtimes(dKE0th, reshape(Uall,ndof,1,numele)), ndof, numele);
quadK  = sum(Lall .* KU,  1)';   % le' * Ke0    * Ue, one value per element
quaddK = sum(Lall .* dKU, 1)';   % le' * dKE0th * Ue, one value per element
xdens_all = xphy(1:numele);
dgh_dx     = fac .* dHdx  - (penal ./ xdens_all) .* quadK;
dgh_dtheta = fac .* dHdth - quaddK;
end

% function h = heavisideBlend(a, b, c, delta)
% % Eq. 12: smooth, compact-support cubic blend between allowable a (c>delta)
% % and b (c<-delta); C1-continuous at the band edges (verified symbolically).
% if c > delta
%     h = a;
% elseif c < -delta
%     h = b;
% else
%     h = 0.75*(a-b)*(c/delta - c^3/(3*delta^3)) + (a+b)/2;
% end
% end
% 
% function dh = dHeavisideBlend(a, b, c, delta)
% if c > delta || c < -delta
%     dh = 0;
% else
%     dh = 0.75*(a-b)*(1/delta - c^2/delta^3);
% end
% end