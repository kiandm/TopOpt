function KE = assembleElementStiffnessWithFixedC(xe, ye, C_fixed, gp, w)
    KE = zeros(8,8);
    for ixi = 1:2
        for ieta = 1:2
            [~, dNdxi] = shapeFcnQ4(gp(ixi), gp(ieta));
            J    = [dNdxi(1,:)*xe, dNdxi(1,:)*ye;
                    dNdxi(2,:)*xe, dNdxi(2,:)*ye];
            dNdx = (J\eye(2)) * dNdxi;
            B = zeros(3,8);
            for a = 1:4
                B(:,2*a-1:2*a) = [dNdx(1,a), 0; 0, dNdx(2,a); dNdx(2,a), dNdx(1,a)];
            end
            KE = KE + B' * C_fixed * B * abs(det(J)) * w(ixi) * w(ieta);
        end
    end
end