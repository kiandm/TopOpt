function [v, dv_dx, dv_theta] = volume_constraint(xphy, volfrac, numele, ve)

    % Compute the total volume
    % v = sum(x)/numele - volfrac;
    x = xphy(1:numele);

    v=sum(x.*ve)/(volfrac*sum(ve))-1;  

    % Sensitivity of the volume constraint with respect to density
    dv_dx = ve/(volfrac*sum(ve));

    % Sensitivity of the volume constraint with respect to fiber orientation
    % For volume constraint, the sensitivity with respect to fiber orientation is zero
    dv_theta = zeros(numele,1);
end
