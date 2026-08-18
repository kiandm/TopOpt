
function [g_tw, dgtw_dx, dgtw_dtheta, TW, vonMises] = Copy_of_TsaiWu( ...
    U, K, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength, freedofs)
% Tsai–Wu stress constraint with adjoint sensitivities (at centre of element)
% Strength allowables 
Xt = strength.Xt; Xc = strength.Xc;
Yt = strength.Yt; Yc = strength.Yc;
S  = strength.S;

% Tsai-Wu criteria constants
F1  = 1/Xt - 1/Xc; F2  = 1/Yt - 1/Yc;
F11 = 1/(Xt*Xc);  F22 = 1/(Yt*Yc);
F66 = 1/S^2;      F12 = -0.5*sqrt(F11*F22);

% Material stiffness
E1 = matprop.E1; E2 = matprop.E2; nu12 = matprop.nu12; 
nu21 = matprop.nu21; G12 = matprop.G12;
C0 = [ E1/(1-nu12*nu21), nu21*E1/(1-nu12*nu21), 0;
       nu12*E2/(1-nu12*nu21), E2/(1-nu12*nu21), 0;
       0, 0, G12];

% Initialisation
TW = zeros(numele,1); dTWdx = zeros(numele,1); dTWdth = zeros(numele,1);
fadj_elem = cell(numele,1); dKE0th = cell(numele,1); vonMises = zeros(numele,1);
q = 0.8; ndof = size(edofMat,2);

for e = 1:numele
    xe = coords(1,conn(:,e))'; ye = coords(2,conn(:,e))';
    Ue = U(edofMat(e,:)); xdens = xphy(e); theta = xphy(numele + e);
    
    % Evaluate at centre (xi=0, eta=0)
    dN = 0.25 * [-1 1 1 -1; -1 -1 1 1];
    J = [dN(1,:)*xe, dN(1,:)*ye; dN(2,:)*xe, dN(2,:)*ye];
    dNdx = J \ dN;
    B = zeros(3,ndof);
    B(1,1:2:end) = dNdx(1,:); B(2,2:2:end) = dNdx(2,:);
    B(3,1:2:end) = dNdx(2,:); B(3,2:2:end) = dNdx(1,:);
    
    c = cos(theta); s = sin(theta);
    T_eps = [c^2, s^2, c*s; s^2, c^2, -c*s; -2*c*s, 2*c*s, c^2-s^2];
    Tinv  = [c^2, s^2, -2*c*s; s^2, c^2, 2*c*s; c*s, -c*s, c^2-s^2];
    
    eps_l = T_eps * (B * Ue);
    sig = xdens^q * (C0 * eps_l);
    
    % Tsai-Wu Index
    TW(e) = F1*sig(1) + F2*sig(2) + F11*sig(1)^2 + F22*sig(2)^2 + F66*sig(3)^2 + 2*F12*sig(1)*sig(2);
    
    % Adjoint & Sensitivities
    psi = [F1 + 2*F11*sig(1) + 2*F12*sig(2); F2 + 2*F22*sig(2) + 2*F12*sig(1); 2*F66*sig(3)];
    dTWdx(e) = psi' * (q * xdens^(q-1) * (C0 * eps_l));
    
    dT_eps_dth = [-sin(2*theta), sin(2*theta), cos(2*theta); sin(2*theta), -sin(2*theta), -cos(2*theta); -2*cos(2*theta), 2*cos(2*theta), -2*sin(2*theta)];
    dTWdth(e)  = psi' * (xdens^q * C0 * dT_eps_dth * (B * Ue));
    fadj_elem{e} = B' * T_eps' * C0' * psi * xdens^q;
    
    % Stiffness derivative (scaled by area)
    dTinv_dth = [-sin(2*theta), sin(2*theta), -2*cos(2*theta); sin(2*theta), -sin(2*theta), 2*cos(2*theta); cos(2*theta), -cos(2*theta), -2*sin(2*theta)];
    dCxy_dth = dTinv_dth * C0 * Tinv' + Tinv * C0 * dTinv_dth';
    dKE0th{e} = det(J) * (B' * dCxy_dth * B);
end

% P-norm Aggregation
p = 8;
TWp = (sum(TW.^p))^(1/p);
g_tw = TWp - 1;
fac = (TW.^(p-1)) / (TWp^(p-1));

% Adjoint solve
fadj = zeros(size(U));
for e = 1:numele
    fadj(edofMat(e,:)) = fadj(edofMat(e,:)) + fac(e) * fadj_elem{e};
end
lambda = zeros(size(U));
lambda(freedofs) = K(freedofs,freedofs) \ fadj(freedofs);

% Final Sensitivities
dgtw_dx = zeros(numele,1); dgtw_dtheta = zeros(numele,1);
for e = 1:numele
    xdens = xphy(e); Ue = U(edofMat(e,:)); le = lambda(edofMat(e,:));
    dgtw_dx(e)     = fac(e)*dTWdx(e) - (penal/xdens) * (le' * KE0{e} * Ue);
    dgtw_dtheta(e) = fac(e)*dTWdth(e) - xdens^penal * (le' * dKE0th{e} * Ue);
end
end