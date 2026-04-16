
function [g_tw, dgtw_dx, dgtw_dtheta, TW, sigmax] = TsaiWu( ...
    U, K, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength)
%--------------------------------------------------------------------------
% Tsai–Wu stress constraint with adjoint sensitivities
% Fully consistent formulation (explicit + adjoint terms)
%
% Inputs:
%   U        : global displacement vector
%   K        : global stiffness matrix
%   KE0      : cell array of unscaled element stiffness matrices
%   xphy     : physical design variables [x; theta]
%--------------------------------------------------------------------------

% Strength allowables 
Xt = strength.Xt; Xc = strength.Xc;
Yt = strength.Yt; Yc = strength.Yc;
S  = strength.S;

F1  = 1/Xt - 1/Xc;
F2  = 1/Yt - 1/Yc;
F11 = 1/(Xt*Xc);
F22 = 1/(Yt*Yc);
F66 = 1/S^2;
F12 = -0.5*sqrt(F11*F22);

% Material stiffness (material coordinates) 
E1   = matprop.E1;
E2   = matprop.E2;
nu12 = matprop.nu12;
nu21 = matprop.nu21;
G12  = matprop.G12;

C0 = [ E1/(1-nu12*nu21),  nu21*E1/(1-nu12*nu21), 0;
       nu12*E2/(1-nu12*nu21), E2/(1-nu12*nu21),  0;
       0,                   0,                  G12 ];

% Initialisation 
TW        = zeros(numele,1);
dTWdx    = zeros(numele,1);
dTWdth   = zeros(numele,1);
sigmax   = zeros(numele,1);

ndof = size(edofMat,2);
gp   = [-1 1]/sqrt(3);

% Store element adjoint RHS contributions
fadj_elem = cell(numele,1);

% Stress interpolation exponent
q = 0.8;

% element loop (no aggregation yet) 
for e = 1:numele

    xe = coords(1,conn(:,e))';
    ye = coords(2,conn(:,e))';
    Ue = U(edofMat(e,:));

    xdens = xphy(e);
    theta = xphy(numele + e);

    c = cos(theta); s = sin(theta);
    T_eps = [ c^2, s^2,  c*s;
              s^2, c^2, -c*s;
             -2*c*s, 2*c*s, c^2-s^2 ];

    TW_e   = 0.0;
    fadj_e = zeros(ndof,1);
    gcount = (e-1)*4;

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

            % Strain/stress
            eps_l = T_eps * (B * Ue);
            sig   = C0 * eps_l;

            s1 = sig(1); s2 = sig(2); t12 = sig(3);
            sigmax(e) = max(sigmax(e), sqrt(s1^2 + s2^2 + 2*t12^2));

            % Tsai–Wu at Gauss point
            TW_gp = F1*s1 + F2*s2 + ...
                    F11*s1^2 + F22*s2^2 + ...
                    F66*t12^2 + 2*F12*s1*s2;

            TW_e = TW_e + TW_gp * wt * jac;

            % adjoint RHS (local) 
            psi = [ F1  + 2*F11*s1 + 2*F12*s2;
                    F2  + 2*F22*s2 + 2*F12*s1;
                    2*F66*t12 ];

            fadj_e = fadj_e + ...
                B' * T_eps' * C0' * psi * wt * jac;
        end
    end

    % Stress interpolation
    TW(e)     = xdens^q * TW_e;
    dTWdx(e)  = q * xdens^(q-1) * TW_e;

    fadj_elem{e} = xdens^q * fadj_e;
end

% p‑norm aggregation
p = 8;
TWp  = (mean(TW.^p))^(1/p);
g_tw = TWp - 1;

fac = (TW.^(p-1)) / (numele * TWp^(p-1));

% assemble adjoint RHS 
fadj = zeros(size(U));
for e = 1:numele
    fadj(edofMat(e,:)) = fadj(edofMat(e,:)) + fac(e) * fadj_elem{e};
end

% adjoint solve
lambda = K \ fadj;

% final sensitivities
dgtw_dx = zeros(numele,1);

for e = 1:numele
    xdens = xphy(e);
    Ue = U(edofMat(e,:));
    le = lambda(edofMat(e,:));
    Ke0 = KE0{e};

    dgtw_dx(e) = fac(e)*dTWdx(e) ...
               - penal * xdens^(penal-1) * (le' * Ke0 * Ue);
end

dgtw_dtheta = zeros(numele,1);  % θ‑adjoint not yet implemented

end
