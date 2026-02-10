function mms = mmsFunctions(E, nu)
    % Define all MMS functions in a structure
    mms.E = E;
    mms.nu = nu;
    mms.Estar = E/(1 - nu^2);
    
    % Displacements (Trigonometric)
    mms.ux_exact = @(x,y) sin(2*pi*x) .* sin(2*pi*y);
    mms.uy_exact = @(x,y) sin(2*pi*x) .* sin(2*pi*y);
    
    % Strains
    mms.strain_xx = @(x,y) 2*pi*cos(2*pi*x) .* sin(2*pi*y);
    mms.strain_yy = @(x,y) 2*pi*sin(2*pi*x) .* cos(2*pi*y);
    mms.strain_xy = @(x,y) 2*pi*(cos(2*pi*x).*sin(2*pi*y) + sin(2*pi*x).*cos(2*pi*y));
    
    % Stresses
    mms.sigma_xx = @(x,y) 2*pi*mms.Estar * (cos(2*pi*x).*sin(2*pi*y) + nu*sin(2*pi*x).*cos(2*pi*y));
    mms.sigma_yy = @(x,y) 2*pi*mms.Estar * (nu*cos(2*pi*x).*sin(2*pi*y) + sin(2*pi*x).*cos(2*pi*y));
    mms.sigma_xy = @(x,y) pi*mms.Estar*(1-nu) * (cos(2*pi*x).*sin(2*pi*y) + sin(2*pi*x).*cos(2*pi*y));
    
    % Stress derivatives for plane stress
    mms.dsigmaxx_dx = @(x,y)  4*pi^2*mms.Estar * (-sin(2*pi*x).*sin(2*pi*y) + nu*cos(2*pi*x).*cos(2*pi*y));
    mms.dsigmaxy_dy = @(x,y)  2*pi^2*mms.Estar * ((1-nu) * ((cos(2*pi*x).*cos(2*pi*y) - sin(2*pi*x).*sin(2*pi*y))));
    mms.dsigmayy_dy = @(x,y)  4*pi^2*mms.Estar * (nu*cos(2*pi*x).*cos(2*pi*y) - sin(2*pi*x).*sin(2*pi*y));
    mms.dsigmaxy_dx = @(x,y)  2*pi^2*mms.Estar * ((1-nu) * ((cos(2*pi*x).*cos(2*pi*y) - sin(2*pi*x).*sin(2*pi*y))));

    % Body forces
    mms.body_force_x = @(x,y) -4*pi^2*mms.Estar * (((nu+1)/2)*cos(2*pi*x).*cos(2*pi*y) - ((3-nu)/2)*sin(2*pi*x).*sin(2*pi*y));
    mms.body_force_y = @(x,y) -4*pi^2*mms.Estar * (((nu+1)/2)*cos(2*pi*x).*cos(2*pi*y) - ((3-nu)/2)*sin(2*pi*x).*sin(2*pi*y));
end