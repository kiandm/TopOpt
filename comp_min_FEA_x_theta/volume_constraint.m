function [v, dv_dx, dv_theta] = volume_constraint(x, volfrac, numele, ve, H, Hs)

    %desing variables
    xphy=x(1:numele); 
    xphy = H*(xphy./Hs);

    %normalised volume constraint
    v=sum(xphy(1:numele).*ve)/(volfrac*sum(ve)) - 1;  

    % Sensitivity of the volume constraint with respect to density
    dv_dx = ve;

    % Sensitivity of the volume constraint with respect to fiber orientation
    % For volume constraint, the sensitivity with respect to fiber orientation is zero
    dv_theta = zeros(numele,1);
end
