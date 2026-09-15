function [g, dg_dx, dg_dtheta, FailIdx, vonMises] = Hashin( ...
    U, dK, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs, dphix_ref, dphiy_ref)
% Mode-separated Hashin stress constraints with adjoint sensitivities.
% Each of the four classical Hashin sub-modes - fibre tension (ft), fibre
% compression (fc), matrix tension (mt), matrix compression (mc) - is
% aggregated independently over the mesh via a p-norm, giving four
% separate constraints g = [g_ft; g_fc; g_mt; g_mc] for the optimiser
% instead of one combined index.
%
% Within a mode, tension/compression selection uses a smooth sigmoid
% gate on the driving stress component (s1 for fibre modes, s2 for
% matrix modes) rather than a hard branch, so a mode's gated index
% decays smoothly to ~0 outside its own physical stress regime while
% staying differentiable everywhere.

% Strength allowables
Xt = strength.Xt; Xc = strength.Xc;
Yt = strength.Yt; Yc = strength.Yc;
S  = strength.S;   % NOTE: in-plane shear strength reused as transverse (S23) allowable
                    % for the matrix-compression term - standard simplification, flag in methods.

k_gate = 50;   % sigmoid sharpness for tension/compression gating within a mode (tune)
q      = 0.8;  % stress interpolation exponent

% Material stiffness (material coordinates)
E1   = matprop.E1; E2   = matprop.E2; nu12 = matprop.nu12;
nu21 = matprop.nu21; G12  = matprop.G12;
C0 = [ E1/(1-nu12*nu21), nu21*E1/(1-nu12*nu21),   0;
       nu12*E2/(1-nu12*nu21), E2/(1-nu12*nu21),    0;
       0,                                    0, G12];

nmode = 4; % [1]=fibre tension [2]=fibre compression [3]=matrix tension [4]=matrix compression

% Initialisation
HashinIdx = zeros(numele,nmode); dHdx = zeros(numele,nmode);
dHdth     = zeros(numele,nmode); vonMises = zeros(numele,1);
ndof = size(edofMat,2);
dKE0th = zeros(ndof,ndof,numele);
fadj_elem = zeros(numele,ndof,nmode);

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

    H_e = zeros(nmode,1); dHdx_e = zeros(nmode,1); dH_dth_e = zeros(nmode,1);
    fadj_e = zeros(ndof,nmode);
    dKE_dth = zeros(ndof,ndof); gcount = (e-1)*4;

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

            %% ---- classical Hashin sub-mode indices, smoothly gated by stress sign ----
            g1 = 1/(1+exp(-k_gate*s1));   % ~1 in fibre tension, ~0 in fibre compression
            g2 = 1/(1+exp(-k_gate*s2));   % ~1 in matrix tension, ~0 in matrix compression

            I_ft = (s1/Xt)^2 + (t12/S)^2;                                    % fibre tension
            I_fc = (s1/Xc)^2;                                                % fibre compression
            I_mt = (s2/Yt)^2 + (t12/S)^2;                                    % matrix tension
            I_mc = (s2/(2*S))^2 + ((Yc/(2*S))^2 - 1)*(s2/Yc) + (t12/S)^2;    % matrix compression

            Igate = [g1*I_ft; (1-g1)*I_fc; g2*I_mt; (1-g2)*I_mc];
            H_e = H_e + Igate * wt * jac;

            %% ---- derivatives of each gated index w.r.t. (s1, s2, t12) ----
            dg1_ds1 = k_gate*g1*(1-g1);
            dg2_ds2 = k_gate*g2*(1-g2);

            dIft_ds1 = 2*s1/Xt^2;  dIft_dt12 = 2*t12/S^2;
            dIfc_ds1 = 2*s1/Xc^2;
            dImt_ds2 = 2*s2/Yt^2;  dImt_dt12 = 2*t12/S^2;
            dImc_ds2 = s2/(2*S^2) + ((Yc/(2*S))^2 - 1)/Yc;  dImc_dt12 = 2*t12/S^2;

            % Psi columns: d(gated index)/d[s1;s2;t12], one column per mode [ft fc mt mc]
            Psi = zeros(3,nmode);
            Psi(1,1) = dg1_ds1*I_ft + g1*dIft_ds1;         Psi(3,1) = g1*dIft_dt12;
            Psi(1,2) = -dg1_ds1*I_fc + (1-g1)*dIfc_ds1;
            Psi(2,3) = dg2_ds2*I_mt + g2*dImt_ds2;         Psi(3,3) = g2*dImt_dt12;
            Psi(2,4) = -dg2_ds2*I_mc + (1-g2)*dImc_ds2;    Psi(3,4) = (1-g2)*dImc_dt12;

            %% sensitivities
            % Density
            dsig_dx = q * xdens^(q-1) * sig_unscaled;
            dHdx_gp = Psi' * dsig_dx;              % nmode x 1
            dHdx_e  = dHdx_e + dHdx_gp * wt * jac;

            % Theta
            dH_dth_gp = Psi' * (xdens^q * C0 * dT_eps_dth * (B * Ue));  % nmode x 1
            dH_dth_e  = dH_dth_e + dH_dth_gp * wt * jac;

            % Adjoint RHS contribution
            fadj_gp = B' * T_eps' * C0' * Psi * xdens^q;   % ndof x nmode
            fadj_e = fadj_e + fadj_gp * wt * jac;
        end
    end

    HashinIdx(e,:) = H_e';
    dHdx(e,:)      = dHdx_e';
    dHdth(e,:)     = dH_dth_e';
    fadj_elem(e,:,:) = reshape(fadj_e, [1,ndof,nmode]);
    dKE0th(:,:,e)    = dKE_dth;
end

% p-norm aggregation over elements, independently per mode
p = 8; % 8 16 32
Hp  = (sum(HashinIdx.^p, 1)).^(1/p);      % 1 x nmode
g   = Hp' - 1;                            % nmode x 1
fac = (HashinIdx.^(p-1)) ./ (Hp.^(p-1));  % numele x nmode

FailIdx = HashinIdx; % per-element, per-mode aggregated index (for diagnostics/plotting)

% assemble adjoint RHS, one column per mode
fadj = zeros(size(U,1), nmode);
for k = 1:nmode
    contrib = fac(:,k) .* fadj_elem(:,:,k);   % numele x ndof
    fadj(:,k) = accumarray(edofMat(:), reshape(contrib, [], 1), [size(U,1) 1]);
end

% adjoint solve - single multi-RHS solve reusing the factorised stiffness
lambda = zeros(size(U,1), nmode);
lambda(freedofs,:) = dK \ fadj(freedofs,:);

% final sensitivities
Uall = U(edofMat)';       % ndof x numele
KU   = reshape(pagemtimes(KE0,    reshape(Uall,ndof,1,numele)), ndof, numele);
dKU  = reshape(pagemtimes(dKE0th, reshape(Uall,ndof,1,numele)), ndof, numele);

xdens_all = xphy(1:numele);
dg_dx = zeros(numele,nmode); dg_dtheta = zeros(numele,nmode);
for k = 1:nmode
    lambda_k = lambda(:,k);
    Lall_k = lambda_k(edofMat)';        % ndof x numele

    quadK  = sum(Lall_k .* KU,  1)';    % numele x 1
    quaddK = sum(Lall_k .* dKU, 1)';    % numele x 1

    dg_dx(:,k)     = fac(:,k).*dHdx(:,k)  - (penal./xdens_all).*quadK;
    dg_dtheta(:,k) = fac(:,k).*dHdth(:,k) - quaddK;
end
end
