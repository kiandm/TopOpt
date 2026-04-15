function KE = assembleElementStiffnessQ4(xe, ye, C, gp, w, theta)
    KE = zeros(8,8);
    for ixi = 1:2
        xi = gp(ixi); wx = w(ixi);
        for ieta = 1:2
            eta = gp(ieta); wy = w(ieta);
            [N, dNdxi] = shapeFcnQ4(xi, eta);
            % Jacobian & inverse
            J = [dNdxi(1,:)*xe, dNdxi(1,:)*ye;
                 dNdxi(2,:)*xe, dNdxi(2,:)*ye];
            detJ = det(J);          % determinant
            invJ = J \ eye(2);
            dNdx = invJ * dNdxi;    % 2x4
            % B matrix (3x8)
            B = zeros(3,8);
            for a = 1:4
                B(:,2*a-1:2*a) = [ dNdx(1,a), 0;
                                   0,         dNdx(2,a);
                                   dNdx(2,a), dNdx(1,a) ];
            end
            % Element elastic matrix
             c = cos(theta);  s = sin(theta);
             T = [   c.^2,   s.^2,   2*c.*s;
                        s.^2,   c.^2,  -2*c.*s;
                       -c.*s,   c.*s,  c.^2-s.^2   ];
            
             R = diag([1, 1, 2]);
             Tinv = inv(T);
             Rinv = inv(R);
             C_rot = Tinv * C * R * T * Rinv;

             KE = KE + (B' * C_rot * B) * abs(detJ) * (wx * wy);

        end
    end
end
