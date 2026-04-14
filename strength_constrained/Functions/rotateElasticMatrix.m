function De_rot = rotateElasticMatrix(De, theta)
% Rotates plane-stress elastic matrix by angle theta (radians)
% Uses Voigt-notation stress transformation T such that:
%   sigma_global = T^{-1} * sigma_local
%   De_rot = T^{-1} * De * T^{-T}

    c = cos(theta);  s = sin(theta);
    c2 = c^2;  s2 = s^2;  cs = c*s;

    % Stress transformation matrix (Voigt, plane stress)
    T = [ c2,    s2,    2*cs;
          s2,    c2,   -2*cs;
         -cs,    cs,  c2-s2 ];

    % Reuter matrix (converts engineering shear strain to tensor shear)
    R = diag([1, 1, 2]);

    % Rotated stiffness: De_rot = T^{-1} * De * R * T * R^{-1}
    % (standard formula for rotating compliance/stiffness in Voigt notation)
    De_rot = (T \ De) / (R * T / R);   % equivalent to inv(T)*De*inv(T')
end