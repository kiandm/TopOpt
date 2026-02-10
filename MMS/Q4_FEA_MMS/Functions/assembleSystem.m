function [K, F] = assembleSystem(nodes, enodes, edof, De, mms)
    ndof = 2 * size(nodes, 1);
    K = sparse(ndof, ndof);
    F = sparse(ndof, 1);
    
    % Gauss quadrature
    gp = [-1, 1]/sqrt(3);
    w = [1, 1];
    
    % Loop over elements
    for e = 1:size(enodes, 1)
        [KE, FE] = assembleElement(e, nodes, enodes, edof, De, mms, gp, w, @shapeFcn);
        
        % Global assembly
        dofs = edof(e,:);
        K(dofs, dofs) = K(dofs, dofs) + KE;
        F(dofs) = F(dofs) + FE;
    end
end

function [KE, FE] = assembleElement(e, nodes, enodes, edof, De, mms, gp, w, shapeFcn)
    dofs = edof(e,:);
    n    = enodes(e,:);
    xe   = nodes(n,2);
    ye   = nodes(n,3);
    % Initialise element stiffness matrix
    KE = zeros(8,8);
    % Initialise element force vector
    FE = zeros(8,1);
    % For each of the 2 gauss points and weights define the natural coords
    for ixi = 1:2
        xi = gp(ixi); wx = w(ixi);
        for ieta = 1:2
            eta = gp(ieta); wy = w(ieta);
            
            % Define the derivatives of the 4 bilinear shape functions
            [N,dNdxi] = shapeFcn(xi,eta);

            % Jacobian (2x2)
            J = [dNdxi(1,:)*xe, dNdxi(1,:)*ye;
                 dNdxi(2,:)*xe, dNdxi(2,:)*ye];
            detJ = abs(det(J));

            % dN/dx, dN/dy
            dNdx = J \ dNdxi; % 2x4

            % B matrix (3x8)
            % B=[B1 B2 B3 B4] where each Bi is 3x2 -> 3x8
            B = zeros(3,8); % Initialise
            for i = 1:4    % Loop over the 4 nodes of the element
                % Fill columns for node i: columns (2*i-1) and (2*i) 
                % correspond to Ux_i and Uy_i so for all rows, each 2 columns 
                B(:,2*i-1:2*i) = [ dNdx(1,i),            0;
                                   0,            dNdx(2,i);
                                   dNdx(2,i),  dNdx(1,i) ];
            end

            KE = KE + (B' * De * B) * detJ * (wx * wy);

            % Adding the body forces for MMS
            % Real coordinates of gauss points
            x_gp = N * xe;
            y_gp = N * ye;

            bx = mms.body_force_x(x_gp, y_gp);
            by = mms.body_force_y(x_gp, y_gp);

            % Distributing forces to the nodes
            for i = 1:4
                FE(2*i-1) = FE(2*i-1) + N(i) * bx * detJ * wx * wy;
                FE(2*i)   = FE(2*i)   + N(i) * by * detJ * wx * wy;
            end
        end
    end
end