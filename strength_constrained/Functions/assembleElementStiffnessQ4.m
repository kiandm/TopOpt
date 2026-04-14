function KE = assembleElementStiffnessQ4(e,nel,xe, ye, C, gp, w, E1, E2, nu12, ...
                                                    nu21, G12, theta )
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
            De = [1/E1,   -nu12/E2,       0; 
                 -nu12/E2,    1/E2,       0; 
                  0,             0,  1/G12]; 
            %for e=1:nel
                % Rotation matrix
             T_theta = [cos(theta(e)).^2, sin(theta(e)).^2, -sin(2*theta(e));
                       sin(theta(e)).^2, cos(theta(e)).^2, sin(2*theta(e));
                       sin(theta(e)).*cos(theta(e)), -sin(theta(e)).*cos(theta(e)), cos(2*theta(e))];
            
             R = diag([1, 1, 2]);
             %De_rot = T_theta * De * T_theta';
                %De_rot = (T_theta / De) / (R .* T_theta / R );
             C_rot = (T_theta \ C) * R * T_theta / R;

             KE = KE + (B' * C_rot * B) * abs(detJ) * (wx * wy);
            %end
        end
    end
end
