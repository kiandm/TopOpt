function [comp, dc_dx, dc_theta] = objective_function(U, xphy, penal, numele, gs, edofMat, coords, conn, matprop)
    
    % Update design variables
    x = xphy(1:numele);
    theta = xphy(numele+1:end);

    % Initialize objective and sensitivities
    comp = 0.0;
    dc_dx = zeros(numele,1);
    dc_theta = zeros(numele,1);

    % Material properties
    E1 = matprop.E1; % Young's modulus in fiber direction
    E2 = matprop.E2; % Young's modulus perpendicular to fiber direction
    nu12 = matprop.nu12; % Major Poisson's ratio
    nu21 = matprop.nu21; % Minor Poisson's ratio
    G12 = matprop.G12; % Shear modulus
    
    % Compute elastic stiffness matrix in local matrical coordinates 1,2
    C12 = [E1/(1 - nu12 * nu21), nu12 * E2 / (1 - nu12 * nu21), 0;
         nu21 * E1 / (1 - nu12 * nu21), E2 / (1 - nu12 * nu21), 0;
         0, 0, G12];

    gcount=0; 
    for ee=1:numele
        % Compute rotation matrix
        c = cos(theta(ee));
        s = sin(theta(ee));

        Tinv = [c^2,        s^2,     -2*c*s;
                s^2,        c^2,     2*c*s;
                c*s,        -c*s,    c^2-s^2];
     
       % Rotate the elastic stiffness matrix C12 from material to global
       % coordiantes x,y
       Cxy= Tinv * C12 * Tinv';
        
       % Derivative of Tinv with respect to theta
        dTinv_dtheta = [-sin(2*theta(ee)),  sin(2*theta(ee)), -2*cos(2*theta(ee));
                         sin(2*theta(ee)), -sin(2*theta(ee)),  2*cos(2*theta(ee));
                         cos(2*theta(ee)), -cos(2*theta(ee)), -2*sin(2*theta(ee))];        
        dCxy_dtheta=dTinv_dtheta * C12 * Tinv' + Tinv * C12 * dTinv_dtheta';

        KE_s=zeros(8,8); %element stiffness for solid material
        dKE_dtheta = zeros(8, 8); % Initialize element stiffness matrix

        for ii=1:4
            gcount=gcount+1; gg=gs(:,gcount); 
            weight=gg(6); jac=gg(7);  
            [phi, dphix, dphiy]=SF_FE(gg,coords,conn);
            Bmat=zeros(3,8);
            Bmat(1,1:2:end)=dphix;  Bmat(2,2:2:end)=dphiy;
            Bmat(3,1:2:end)=dphiy;  Bmat(3,2:2:end)=dphix;
            KE_s=KE_s+jac*weight*Bmat'*Cxy*Bmat;
            dKE_dtheta=dKE_dtheta + jac*weight*Bmat'*dCxy_dtheta*Bmat;
        end 

        % Element displacement vector
        Ue=U(edofMat(ee,:));

        % Compute element compliance for solid element
        ce=Ue'*KE_s*Ue;

        % Total compliance with using E_eff
        comp=comp+x(ee)^penal * ce;

        % Sensitivity of compliance with respect to density
        dc_dx(ee) = -penal * x(ee)^(penal - 1) * ce;

        % Sensitivity of compliance with respect to theta
        dc_theta(ee) = -x(ee)^penal * Ue'*dKE_dtheta*Ue;
    end 
end