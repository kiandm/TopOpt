function [v_con, dv_dx, dv_theta] = volume_const(xphy, volfrac, numele, ve)

    %To impose equility constraint for voluem Ax-b=0;
    %we will need two constraints, i.e. Ax-b<=0 and Ax-b>=0 (or -Ax+b<=0)

    x = xphy(1:numele);
    %normalised volume constraint
    v_con=sum(x.*ve)/(volfrac*sum(ve)) - 1;   % Ax-b<=0

    % Sensitivity of the volume constraint with respect to density
    dv_dx = ve/(volfrac*sum(ve));
    % dv_dx = ve;

    % Sensitivity of the volume constraint with respect to fiber orientation
    % For volume constraint, the sensitivity with respect to fiber orientation is zero
    dv_theta = zeros(numele,1);
end
