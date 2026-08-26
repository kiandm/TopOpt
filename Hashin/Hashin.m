function [g_h, dgh_dx, dgh_dtheta, HashinIdx, HashinIdx_gp, vonMises] = Hashin( ...
    U, dK, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs, dphix_ref, dphiy_ref)
% Classical (non-mode-separated) Hashin stress constraint with adjoint sensitivities
% Fibre and matrix modes combined per Gauss point via KS smooth-max;
% tension/compression within each mode gated via sigmoid (differentiable Macaulay).
% Mode separation (keeping fibre/matrix as distinct constraints) is a later extension.

% Strength allowables
Xt = strength.Xt; Xc = strength.Xc;
Yt = strength.Yt; Yc = strength.Yc;
S  = strength.S;   % NOTE: in-plane shear strength reused as transverse (S23) allowable
                    % for the matrix-compression term - standard simplification, flag in methods.

% Sharpness parameters 
k_gate = 50;   % sigmoid sharpness for tension/compression gating within a mode
k_ks   = 20;   % KS aggregation sharpness for combining fibre & matrix modes

% Material stiffness (material coordinates)
E1   = matprop.E1; E2   = matprop.E2; nu12 = matprop.nu12;
nu21 = matprop.nu21; G12  = matprop.G12;
C0 = [ E1/(1-nu12*nu21), nu21*E1/(1-nu12*nu21),   0;
       nu12*E2/(1-nu12*nu21), E2/(1-nu12*nu21),    0;
       0,                                    0, G12];

% Initialisation
HashinIdx = zeros(numele,1); dHdx    = zeros(numele,1);
dHdth     = zeros(numele,1); vonMises = zeros(numele,1);
dKE0th = cell(numele,1); ndof = size(edofMat,2);
%gp   = [-1 1]/sqrt(3);

fadj_elem = cell(numele,1);
q = 0.8; % stress interpolation exponent 

% Element loop
for e = 1:numele

    %xe = coords(1,conn(:,e))'; ye = coords(2,conn(:,e))';
    Ue = U(edofMat(e,:)); xdens = xphy(e); theta = xphy(numele + e);
    c = cos(theta); s = sin(theta);
    T_eps = [ c^2, s^2,  c*s;
              s^2, c^2, -c*s;
             -2*c*s, 2*c*s, c^2-s^2 ];
    H_e   = 0.0; dHdx_e = 0.0;
    fadj_e = zeros(ndof,1); dH_dth_e = 0.0;
    dKE_dth   = zeros(8,8); gcount = (e-1)*4;

    for i = 1:2
        for j = 1:2

            gcount = gcount + 1;
            %xi = gp(i); eta = gp(j);
            wt  = gs(6,gcount);
            jac = gs(7,gcount);

            % % Shape derivatives
            % dNdxi = 0.25 * ...
            %     [-(1-eta),  (1-eta),  (1+eta), -(1+eta);
            %      -(1-xi),  -(1+xi),   (1+xi),   (1-xi)];
            % J = [dNdxi(1,:)*xe, dNdxi(1,:)*ye;
            %      dNdxi(2,:)*xe, dNdxi(2,:)*ye];
            % dNdx = J \ dNdxi;
            gp_idx = gcount - (e-1)*4; % local Gauss point index within this element (1..4)
            dNdx = [dphix_ref(:,gp_idx)'; dphiy_ref(:,gp_idx)'];


            % B-matrix
            B = zeros(3,ndof);
            for a = 1:4
                B(:,2*a-1:2*a) = ...
                    [dNdx(1,a), 0;
                     0, dNdx(2,a);
                     dNdx(2,a), dNdx(1,a)];
            end

            Tinv = [c^2,  s^2, -2*c*s;
                    s^2,  c^2,  2*c*s;
                    c*s, -c*s,  c^2-s^2];
            dTinv_dth = [-sin(2*theta),  sin(2*theta), -2*cos(2*theta);
                          sin(2*theta), -sin(2*theta),  2*cos(2*theta);
                          cos(2*theta), -cos(2*theta), -2*sin(2*theta)];
            dCxy_dth = dTinv_dth * C0 * Tinv' + Tinv * C0 * dTinv_dth';
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

            %% ---- Hashin failure index at this Gauss point ----
            % Sub-mode indices
            I_ft = (s1/Xt)^2 + (t12/S)^2;   % fibre tension
            I_fc = (s1/Xc)^2;               % fibre compression
            I_mt = (s2/Yt)^2 + (t12/S)^2;   % matrix tension
            I_mc = (s2/(2*S))^2 + ((Yc/(2*S))^2 - 1)*(s2/Yc) + (t12/S)^2; % matrix compression

            % Smooth tension/compression gates
            g1 = 1/(1+exp(-k_gate*s1));     % ~1 in tension, ~0 in compression
            g2 = 1/(1+exp(-k_gate*s2));

            I_fiber  = g1*I_ft + (1-g1)*I_fc;
            I_matrix = g2*I_mt + (1-g2)*I_mc;

            % KS smooth-max combine (numerically stable form)
            a_ks = k_ks*I_fiber; b_ks = k_ks*I_matrix;
            m_ks = max(a_ks, b_ks);
            H_gp = (m_ks + log(exp(a_ks-m_ks) + exp(b_ks-m_ks))) / k_ks;
            w_f  = exp(a_ks-m_ks) / (exp(a_ks-m_ks) + exp(b_ks-m_ks)); % softmax weight, fibre
            w_m  = 1 - w_f;                                            % softmax weight, matrix

            H_e = H_e + H_gp * wt * jac;

            %% ---- derivatives of H_gp w.r.t. (s1, s2, t12) ----
            dIft_ds1 = 2*s1/Xt^2;   dIft_dt12 = 2*t12/S^2;
            dIfc_ds1 = 2*s1/Xc^2;   dIfc_dt12 = 0;
            dImt_ds2 = 2*s2/Yt^2;   dImt_dt12 = 2*t12/S^2;
            dImc_ds2 = s2/(2*S^2) + ((Yc/(2*S))^2 - 1)/Yc; dImc_dt12 = 2*t12/S^2;

            dIfiber_ds1  = k_gate*g1*(1-g1)*(I_ft - I_fc) + g1*dIft_ds1 + (1-g1)*dIfc_ds1;
            dIfiber_dt12 = g1*dIft_dt12 + (1-g1)*dIfc_dt12;

            dImatrix_ds2  = k_gate*g2*(1-g2)*(I_mt - I_mc) + g2*dImt_ds2 + (1-g2)*dImc_ds2;
            dImatrix_dt12 = g2*dImt_dt12 + (1-g2)*dImc_dt12;

            psi = [ w_f * dIfiber_ds1;
                    w_m * dImatrix_ds2;
                    w_f * dIfiber_dt12 + w_m * dImatrix_dt12 ];

            %% sensitivities
            % Density
            dsig_dx = q * xdens^(q-1) * sig_unscaled;
            dHdx_gp = psi' * dsig_dx;
            dHdx_e  = dHdx_e + dHdx_gp * wt * jac;

            % Theta
            dT_eps_dth = [-sin(2*theta),  sin(2*theta),   cos(2*theta);
                           sin(2*theta), -sin(2*theta),  -cos(2*theta);
                          -2*cos(2*theta), 2*cos(2*theta), -2*sin(2*theta)];
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
    fadj_elem{e} = fadj_e;
    dKE0th{e}    = dKE_dth;
end

% p-norm aggregation over elements 
p = 8; % 8 16 32
Hp  = (sum(HashinIdx.^p))^(1/p);
g_h = Hp - 1;
fac = (HashinIdx.^(p-1)) / (Hp^(p-1));

HashinIdx_gp = HashinIdx; 

% assemble adjoint RHS
fadj = zeros(size(U));
for e = 1:numele
    fadj(edofMat(e,:)) = fadj(edofMat(e,:)) + fac(e) * fadj_elem{e};
end

% adjoint solve
lambda = zeros(size(U));
% lambda(freedofs) = K(freedofs,freedofs) \ fadj(freedofs); % (OLD)
lambda(freedofs) = dK \ fadj(freedofs); % (NEW)


% final sensitivities
dgh_dx = zeros(numele,1);
dgh_dtheta = zeros(numele,1);
for e = 1:numele
    xdens = xphy(e);
    Ue = U(edofMat(e,:));
    le = lambda(edofMat(e,:));
    Ke0 = KE0{e};
    dgh_dx(e) = fac(e)*dHdx(e) ...
               - (penal/xdens) * (le' * Ke0 * Ue);
    dgh_dtheta(e) = fac(e) * dHdth(e) ...
                   - (le' * dKE0th{e} * Ue);
end
end
