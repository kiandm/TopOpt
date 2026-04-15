function [g_tw, dgtw_dx, dgtw_dtheta,TW,sigmax] = TsaiWu(U, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength)
% TSAI-WU criterion stress constraint by Kian Das
% Computes failure index using stresses in material coordinartes

% Strength allowables (X = longitudinal & Y = transverse, t = tensile & c = compressive)
Xt = strength.Xt; Xc = strength.Xc;
Yt = strength.Yt; Yc = strength.Yc;
S  = strength.S;

% Tsai-Wu equations
F1  = 1/Xt - 1/Xc;
F2  = 1/Yt - 1/Yc;
F11 = 1/(Xt*Xc);
F22 = 1/(Yt*Yc);
F66 = 1/S^2;
F12 = -0.5*sqrt(F11*F22);

% Material properties
E1   = matprop.E1;
E2   = matprop.E2;
nu12 = matprop.nu12;
nu21 = matprop.nu21;
G12  = matprop.G12;

% Orthotropic plane-stress stiffness matrix
C0 = [ E1/(1-nu12*nu21),  nu21*E1/(1-nu12*nu21), 0;
       nu12*E2/(1-nu12*nu21), E2/(1-nu12*nu21),  0;
       0,                   0,                  G12 ];

%% --- Initialize ---
TW     = zeros(numele,1); % elementwise Tsai-Wu index
dTWdx  = zeros(numele,1); % density sens
dTWdth = zeros(numele,1); % angle sens (INCOMPLETE)
sigmax = zeros(numele,1); % max stress magnitude per element for plotting

gp = [-1  1]/sqrt(3);
ndof = size(edofMat,2);
%% --- Element loop ---
for e = 1:numele
    % Extract element data
    xe = coords(1, conn(:,e))';
    ye = coords(2, conn(:,e))';

    Ue = U(edofMat(e,:)); % Element displacement vector
    % Desing variables
    theta = xphy(numele + e);
    xdens = xphy(e);

    TW_e = 0;
    dTWdtheta_e = 0; % NEED TO IMPLEMENT THIS NEXT
    gcount = 0; % (e-1)*4
    % Rotation matrices
    c = cos(theta); s = sin(theta);
    T_eps = [ c^2, s^2,  c*s;
              s^2, c^2, -c*s;
             -2*c*s, 2*c*s, c^2-s^2 ];

   
    for i = 1:2
        for j = 1:2

            gcount = gcount + 1;
            xi = gp(i); eta = gp(j);

            wt  = gs(6,gcount);
            jac = gs(7,gcount);

            % Shape functions and derivatives
            dNdxi = 0.25 * ...
              [- (1-eta),  (1-eta),  (1+eta), -(1+eta);
                 - (1-xi), -(1+xi),   (1+xi),   (1-xi)];

            J = [dNdxi(1,:)*xe, dNdxi(1,:)*ye;
                 dNdxi(2,:)*xe, dNdxi(2,:)*ye];

            dNdx = J \ dNdxi;

            % B-matrix
            B = zeros(3,ndof);
            for a = 1:4
                B(:,2*a-1:2*a) = ...
                  [ dNdx(1,a), 0;
                    0, dNdx(2,a);
                    dNdx(2,a), dNdx(1,a) ];
            end

            % Strain and stress
            eps_g = B * Ue;
            eps_l = T_eps * eps_g;
            sig   = C0 * eps_l;

            s1 = sig(1); s2 = sig(2); t12 = sig(3);
            
            sigma_mag = sqrt(s1^2 + s2^2 + 2*t12^2);
            sigmax(e) = max(sigmax(e), sigma_mag);
            % Tsai-Wu evaluation at Gauss point
            TW_gp = F1*s1 + F2*s2 + ...
                    F11*s1^2 + F22*s2^2 + ...
                    F66*t12^2 + 2*F12*s1*s2;

            TW_e = TW_e + TW_gp * wt * jac;
        end
    end

    % SIMP scaling
    q=0.5;
    TW(e) = xdens^q * TW_e;

    dTWdx(e)= q * xdens^(q-1) * TW_e;

    dTWdth(e) = 0;   % NOT DONE
end


% p-norm aggregation 
p = 8;
TWp = (mean(TW.^p))^(1/p);
g_tw = TWp - 1;

% Sensitivities of aggregation
fac = (TW.^(p-1)) / (numele * TWp^(p-1));
dgtw_dx     = fac .* dTWdx;
dgtw_dtheta = fac .* dTWdth;

end