function [principalStress, stressDirection] = calculatePrincipalStress(U, numele, gs, edofMat, coords, conn, matprop)
    % Initialise outputs
    principalStress = zeros(numele, 1);
    stressDirection = zeros(numele, 1);
    
    % Isotropic Constitutive Matrix (D)
    E = matprop.E1; nu = matprop.nu12;
    D = (E/(1-nu^2)) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    
    % Loop over elements
    for e = 1:numele
        % Get element nodes and displacements
        nodes = conn(:,e);
        ue = U(edofMat(e,:));
        
        % Get coordinates of the 4 nodes
        xe = coords(1, nodes); ye = coords(2, nodes);
        
        % Calculate stress at the element center (local xi=0, eta=0)
        % Standard B-matrix for a 4-node quad
        xi = 0; eta = 0;
        dN = 1/4 * [-(1-eta)  (1-eta) (1+eta) -(1+eta);
                    -(1-xi)  -(1+xi)  (1+xi)   (1-xi)];
        J = dN * [xe' ye'];
        % invJ = inv(J);
        % dN_global = invJ * dN;
        dN_global = J \ dN;
        
        B = zeros(3,8);
        B(1,1:2:7) = dN_global(1,:);
        B(2,2:2:8) = dN_global(2,:);
        B(3,1:2:7) = dN_global(2,:);
        B(3,2:2:8) = dN_global(1,:);
        
        % Stress tensor: [sigma_x; sigma_y; tau_xy]
        sigma = D * B * ue;
        
        sx = sigma(1); sy = sigma(2); sxy = sigma(3);
        
        % Principal Stress Angle
        % theta_p = 0.5 * atan( 2*tau_xy / (sigma_x - sigma_y) )
        % use atan2 to handle the quadrants safely
        stressDirection(e) = 0.5 * atan2(2*sxy, sx - sy);
        
        % Major Principal Stress Magnitude (for reference)
        principalStress(e) = (sx + sy)/2 + sqrt(((sx - sy)/2)^2 + sxy^2);
    end
end