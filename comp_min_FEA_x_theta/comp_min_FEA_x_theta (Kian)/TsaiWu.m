
function [g_tw, dgtw_dx, dgtw_dtheta, TW, TW_gp, vonMises] = TsaiWu( ...
    U, K, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs)
% Tsai–Wu stress constraint with adjoint sensitivities
% Strength allowables 
Xt = strength.Xt; Xc = strength.Xc;
Yt = strength.Yt; Yc = strength.Yc;
S  = strength.S;
% Tsai-Wu
F1  = 1/Xt - 1/Xc; F2  = 1/Yt - 1/Yc;
F11 = 1/(Xt*Xc); F22 = 1/(Yt*Yc);
F66 = 1/S^2; F12 = -0.5*sqrt(F11*F22);
% Material stiffness (material coordinates) 
E1   = matprop.E1; E2   = matprop.E2; nu12 = matprop.nu12;
nu21 = matprop.nu21; G12  = matprop.G12;
% Elastic stiffness
C0 = [ E1/(1-nu12*nu21), nu21*E1/(1-nu12*nu21),   0;
       nu12*E2/(1-nu12*nu21), E2/(1-nu12*nu21),    0;
       0,                                    0, G12];

% Initialisation 
TW        = zeros(numele,1); dTWdx    = zeros(numele,1); 
dTWdth   = zeros(numele,1); vonMises = zeros(numele,1);
dKE0th = cell(numele,1); ndof = size(edofMat,2);
gp   = [-1 1]/sqrt(3);

% Store element adjoint RHS contributions
fadj_elem = cell(numele,1);
% Stress interpolation exponent
q = 0.8;

% Element loop 
for e = 1:numele

    xe = coords(1,conn(:,e))'; ye = coords(2,conn(:,e))';
    Ue = U(edofMat(e,:)); xdens = xphy(e); theta = xphy(numele + e);
    c = cos(theta); s = sin(theta);
    T_eps = [ c^2, s^2,  c*s;
              s^2, c^2, -c*s;
             -2*c*s, 2*c*s, c^2-s^2 ];
    TW_e   = 0.0; dTWdx_e = 0.0;
    fadj_e = zeros(ndof,1); dTW_dth_e = 0.0;
    dKE_dth   = zeros(8,8); gcount = (e-1)*4;    

    for i = 1:2
        for j = 1:2

            gcount = gcount + 1;
            xi = gp(i); eta = gp(j);
            wt  = gs(6,gcount);
            jac = gs(7,gcount);

            % Shape derivatives
            dNdxi = 0.25 * ...
                [-(1-eta),  (1-eta),  (1+eta), -(1+eta);
                 -(1-xi),  -(1+xi),   (1+xi),   (1-xi)];
            J = [dNdxi(1,:)*xe, dNdxi(1,:)*ye;
                 dNdxi(2,:)*xe, dNdxi(2,:)*ye];
            dNdx = J \ dNdxi;

            % B-matrix
            B = zeros(3,ndof);
            for a = 1:4
                B(:,2*a-1:2*a) = ...
                    [dNdx(1,a), 0;
                     0, dNdx(2,a);
                     dNdx(2,a), dNdx(1,a)];
            end

            Tinv = [c^2,  s^2, -2*c*s;       % inverse rotation (stress form)
                    s^2,  c^2,  2*c*s;
                    c*s, -c*s,  c^2-s^2];
            dTinv_dth = [-sin(2*theta),  sin(2*theta), -2*cos(2*theta);
                          sin(2*theta), -sin(2*theta),  2*cos(2*theta);
                          cos(2*theta), -cos(2*theta), -2*sin(2*theta)];
            dCxy_dth = dTinv_dth * C0 * Tinv' + Tinv * C0 * dTinv_dth';
            dKE_dth = dKE_dth + jac * wt * B' * dCxy_dth * B;

            % Strain/stress
            eps_l = T_eps * (B * Ue);
            sig_unscaled   = C0 * eps_l;
            sig = xdens^q * sig_unscaled;
            s1 = sig(1); s2 = sig(2); t12 = sig(3);
            sig_global = Tinv * sig;
            sx = sig_global(1); sy = sig_global(2); txy = sig_global(3);
            vm_gp = sqrt(sx^2 - sx*sy + sy^2 + 3*txy^2);
            vonMises(e) = vonMises(e) + vm_gp / 4;

            % Tsai–Wu at Gauss point
            TW_gp = F1*s1 + F2*s2 + ...
                    F11*s1^2 + F22*s2^2 + ...
                    F66*t12^2 + 2*F12*s1*s2;
            TW_e = TW_e + TW_gp / 4; % * wt * jac;

            % adjoint RHS (local) 
            psi = [ F1  + 2*F11*s1 + 2*F12*s2;
                    F2  + 2*F22*s2 + 2*F12*s1;
                    2*F66*t12 ];

            % DENSITY SENSITIVITY (Chain Rule: dTW/dsig * dsig/dx)
            dsig_dx = q * xdens^(q-1) * sig_unscaled;
            dTWdx_gp = psi' * dsig_dx;
            dTWdx_e = dTWdx_e + dTWdx_gp / 4;

            % derivative of T_eps with respect to theta
            dT_eps_dth = [-sin(2*theta),  sin(2*theta),   cos(2*theta);
                           sin(2*theta), -sin(2*theta),  -cos(2*theta);
                          -2*cos(2*theta), 2*cos(2*theta), -2*sin(2*theta)];          
            dTW_dth_gp = psi' * (xdens^q * C0 * dT_eps_dth * (B * Ue));
            dTW_dth_e  = dTW_dth_e + dTW_dth_gp / 4;            
            fadj_gp = B' * T_eps' * C0' * psi * xdens^q;
            fadj_e = fadj_e + fadj_gp / 4;
        end
    end

    TW(e)     = TW_e;
    dTWdx(e)  = dTWdx_e;    
    dTWdth(e) = dTW_dth_e;
    fadj_elem{e} = fadj_e;
    dKE0th{e} = dKE_dth;
end

% p‑norm aggregation
p = 8; % 8 16 32
TWp  = (sum(TW.^p))^(1/p);
g_tw = TWp - 1;
fac = (TW.^(p-1)) / (TWp^(p-1)); % numele *

% assemble adjoint RHS 
fadj = zeros(size(U));
for e = 1:numele
    fadj(edofMat(e,:)) = fadj(edofMat(e,:)) + fac(e) * fadj_elem{e};
end

% adjoint solve
lambda = zeros(size(U));
lambda(freedofs) = K(freedofs,freedofs) \ fadj(freedofs);

% final sensitivities
dgtw_dx = zeros(numele,1);
dgtw_dtheta = zeros(numele,1);
for e = 1:numele
    xdens = xphy(e);
    theta  = xphy(numele + e);
    Ue = U(edofMat(e,:));
    le = lambda(edofMat(e,:));
    Ke0 = KE0{e};
    dgtw_dx(e) = fac(e)*dTWdx(e) ...
               - (penal/xdens) * (le' * Ke0 * Ue);
    dgtw_dtheta(e) = fac(e) * dTWdth(e) ...
                   - xdens^penal * (le' * dKE0th{e} * Ue);
end
end
