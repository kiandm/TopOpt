function out = gpEval(U, xphy, numele, gs, edofMat, matprop, strength, opts, dphix_ref, dphiy_ref, xmin_solid)
% EXPERIMENT: GAUSS-POINT-level failure indices of a design, no aggregation.  Classical Hashin with HARD sign gating
% (ft/fc by the sign of s1, mt/mc by the sign of s2) using the same relaxed stress x^q*C*eps as the optimiser, plus the
% Dong-type index H_gp = smooth-max(FI_f, FI_m).  Max is over the Gauss points of elements with density > xmin_solid.
if nargin < 11, xmin_solid = 0.5; end
Xt = strength.Xt; Xc = strength.Xc; Yt = strength.Yt; Yc = strength.Yc; S = strength.S;
q = opts.q;  delta_f = opts.delta_frac*min(Xt,Xc);  delta_m = opts.delta_frac*min(Yt,Yc);  eps_max = opts.eps_max;
E1 = matprop.E1; E2 = matprop.E2; nu12 = matprop.nu12; nu21 = matprop.nu21; G12 = matprop.G12;
C0 = [E1/(1-nu12*nu21), nu21*E1/(1-nu12*nu21), 0; nu12*E2/(1-nu12*nu21), E2/(1-nu12*nu21), 0; 0 0 G12];
xd = xphy(1:numele);  theta = xphy(numele+1:2*numele);  c = cos(theta);  s = sin(theta);
Uall = U(edofMat)';
solid = xd > xmin_solid;
xq = xd.^q;
FI = zeros(numele, 4, 4);   H = zeros(numele, 4);       % element x mode x Gauss point
for ii = 1:4
    B = zeros(3,8);  dphix = dphix_ref(:,ii)'; dphiy = dphiy_ref(:,ii)';
    B(1,1:2:end) = dphix;  B(2,2:2:end) = dphiy;  B(3,1:2:end) = dphiy;  B(3,2:2:end) = dphix;
    eg = (B * Uall)';   e1 = eg(:,1);  e2 = eg(:,2);  e3 = eg(:,3);
    l1 = c.^2.*e1 + s.^2.*e2 + c.*s.*e3;
    l2 = s.^2.*e1 + c.^2.*e2 - c.*s.*e3;
    l3 = -2*c.*s.*e1 + 2*c.*s.*e2 + (c.^2 - s.^2).*e3;
    s1 = xq.*(C0(1,1)*l1 + C0(1,2)*l2);  s2 = xq.*(C0(2,1)*l1 + C0(2,2)*l2);  t12 = xq.*(C0(3,3)*l3);
    FI(:,1,ii) = (s1 > 0) .* sqrt((s1/Xt).^2 + (t12/S).^2);
    FI(:,2,ii) = (s1 < 0) .* sqrt((s1/Xc).^2 + (t12/S).^2);
    FI(:,3,ii) = (s2 > 0) .* sqrt((s2/Yt).^2 + (t12/S).^2);
    FI(:,4,ii) = (s2 < 0) .* sqrt((s2/Yc).^2 + (t12/S).^2);
    s1c = min(max(s1, -delta_f), delta_f);  Hf = 0.75*(Xt-Xc)*(s1c/delta_f - s1c.^3/(3*delta_f^3)) + (Xt+Xc)/2;
    s2c = min(max(s2, -delta_m), delta_m);  Hm = 0.75*(Yt-Yc)*(s2c/delta_m - s2c.^3/(3*delta_m^3)) + (Yt+Yc)/2;
    FI_f = sqrt((s1./Hf).^2 + (t12/S).^2 + opts.eps_sq);  FI_m = sqrt((s2./Hm).^2 + (t12/S).^2 + opts.eps_sq);
    H(:,ii) = 0.5*((FI_f+FI_m) + sqrt((FI_f-FI_m).^2 + eps_max) - sqrt(eps_max));
end
out.gpMax    = reshape(max(max(FI(solid,:,:), [], 3), [], 1), 1, 4);   % 1x4  worst Gauss point of each mode
out.gpMaxAll = max(out.gpMax);                                       % worst Gauss point of any mode
out.gpMaxH   = max(max(H(solid,:)));                                 % worst Gauss point of the Dong index
out.nGPover  = reshape(sum(sum(FI(solid,:,:) > 1, 3), 1), 1, 4);        % 1x4 number of Gauss points with index > 1
out.nGPoverH = sum(sum(H(solid,:) > 1));
elAvg = mean(FI, 3);                                                 % element average of the HARD-gated index (no sigmoid)
out.elMax    = reshape(max(elAvg(solid,:), [], 1), 1, 4);            % 1x4
out.dilution = out.gpMax ./ out.elMax;                               % >1: the element average hides a worse Gauss point
wm = squeeze(max(FI(solid,:,:), [], 2));  wm = sort(wm(:));            % worst-mode index of every solid Gauss point
out.gpP99 = wm(max(1, ceil(0.99*numel(wm))));  out.gpP95 = wm(max(1, ceil(0.95*numel(wm))));
out.nsolid   = nnz(solid);
end