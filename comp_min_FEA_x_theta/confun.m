function [c, ceq, dc, dceq] = confun(x, volfrac, numele, ve, H, Hs)
    
    %volume constraint
    [v, dv_dx, dv_theta] = volume_constraint(x, volfrac, numele, ve, H, Hs); 

    dv_dx = H*(dv_dx./Hs);
    dv_theta = H*(dv_theta./Hs);
    
    % inquality constraint and its derivative
    c=[]; 
    dc= []; 

    % Equailty constraint and its derivative
    ceq=[v];
    dceq=[dv_dx; dv_theta]; 
end %confun
