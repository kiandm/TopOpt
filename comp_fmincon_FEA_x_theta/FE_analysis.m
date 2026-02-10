function [U] = FE_analysis(x, penal, numnode, numele, gs, edofMat, coords, conn, freedofs, F, matprop, H, Hs)
    % Material properties


    %desing variables
    xphy=x(1:numele); 
    theta=x(numele+1:end);
    
    xphy = H*(xphy./Hs);
    theta = H*(theta./Hs);

    % Material properties for the composite
    E1 = matprop.E1; % Young's modulus in fiber direction
    E2 = matprop.E2; % Young's modulus perpendicular to fiber direction
    nu12 = matprop.nu12; % Major Poisson's ratio
    nu21 = matprop.nu21; % Minor Poisson's ratio
    G12 = matprop.G12; % Shear modulus

    % Compute the reduced stiffness matrix Q
    C12 = [E1/(1 - nu12 * nu21), nu12 * E2 / (1 - nu12 * nu21), 0;
         nu21 * E1 / (1 - nu12 * nu21), E2 / (1 - nu12 * nu21), 0;
         0, 0, G12];

    % Total degrees of freedom
    ndof = 2 * numnode;

    % Initialize global stiffness matrix and force vector
    K = sparse(ndof, ndof);
    % F = sparse(ndof, 1);

    % Apply loads (example: load at the middle of the top edge)
    % F(2,1) = -1;

    gcount=0; 
    for ee=1:numele
        KE=zeros(8,8); 

        % Compute rotation matrix
        c = cos(theta(ee));
        s = sin(theta(ee));
        Tinv = [c^2,        s^2,     -2*c*s;
                s^2,        c^2,     2*c*s;
                c*s,        -c*s,    c^2-s^2];
    
       % Rotate the elastic stiffness matrix C12 from material to global
       % coordiantes x,y
       Cxy= Tinv * C12 * Tinv';
       Cxy = xphy(ee)^penal*Cxy;

        for ii=1:4
            gcount=gcount+1; gg=gs(:,gcount); 
            weight=gg(6); jac=gg(7);  
            [phi, dphix, dphiy]=SF_FE(gg,coords,conn);
            Bmat=zeros(3,8);
            Bmat(1,1:2:end)=dphix;  Bmat(2,2:2:end)=dphiy;
            Bmat(3,1:2:end)=dphiy;  Bmat(3,2:2:end)=dphix;
            KE=KE+jac*weight*Bmat'*Cxy*Bmat;
        end 
        K(edofMat(ee,:),edofMat(ee,:))=K(edofMat(ee,:),edofMat(ee,:))+KE; 
    end 

    % Solve the system
    U = zeros(ndof, 1);
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);

    % figure
    % scale=1e-2;
    % coords_def(1,:)=coords(1,:)+U(1:2:end)'*scale;
    % coords_def(2,:)=coords(2,:)+U(2:2:end)'*scale;
    % 
    % patch('Faces',conn','Vertices',coords_def','FaceVertexCData',U(2:2:end),...
    % 'FaceColor','flat','EdgeColor',[0 0 0]); colorbar; %axis equal off;
    % axis equal;
    
end