function dCrot = dRotateElasticMatrix(C, theta)
% Analytic derivative of the rotated elastic matrix w.r.t. theta

    c = cos(theta);  s = sin(theta);

    % dT/dtheta
    dT = [ -2*c*s,   2*c*s,  2*(c^2-s^2);
            2*c*s,  -2*c*s, -2*(c^2-s^2);
           s^2-c^2, c^2-s^2, -4*c*s      ];

    R = diag([1, 1, 2]);
    T = [ c^2,   s^2,   2*c*s;
          s^2,   c^2,  -2*c*s;
         -c*s,   c*s,  c^2-s^2];

    Tinv  = inv(T);
    Rinv = inv(R);
    RTRi  = R * T * Rinv;   % R*T*R^{-1}

    % Product rule: d/dtheta [ Tinv * De * inv(RTRi) ]
    dTinv   = -Tinv * dT * Tinv;
    dRTRi   = R * dT / R;
    %dDe_dtheta = dTinv * De / RTRi + Tinv * De * (-RTRi \ dRTRi / RTRi);
    dCrot = dTinv * C * RTRi + Tinv * C * dRTRi;
end