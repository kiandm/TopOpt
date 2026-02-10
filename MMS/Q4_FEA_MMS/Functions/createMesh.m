function [nodes, enodes, edof, nelx, nely, nn, ndof] = createMesh(test)
    % Mesh parameters
    nelx = test; 
    nely = test;
    nnx = nelx + 1; 
    nny = nely + 1;
    nn = nnx * nny;
    nel = nelx * nely;
    ndof = 2 * nn;
    
    % Node coordinates and DOFs
    nodes = zeros(nn, 5);
    node = 0;
    for y = 0:nely
        for x = 0:nelx
            node = node + 1;
            nodes(node,1) = node;
            nodes(node,2) = x/nelx;          % x coordinate
            nodes(node,3) = y/nely;          % y coordinate
            nodes(node,4) = 2*node - 1;      % Ux DOF
            nodes(node,5) = 2*node;          % Uy DOF
        end
    end
    
    % Element connectivity
    node_id = @(x,y) (y*(nelx+1) + x + 1);
    enodes = zeros(nel, 4);
    elem = 0;
    for ex = 1:nelx
        for ey = 1:nely
            elem = elem + 1;
            x0 = ex-1; y0 = ey-1;
            n1 = node_id(x0, y0);
            n2 = node_id(x0+1, y0);
            n3 = node_id(x0+1, y0+1);
            n4 = node_id(x0, y0+1);
            enodes(elem,:) = [n1, n2, n3, n4];
        end
    end
    
    % DOF connectivity
    edof = zeros(nel, 8);
    for e = 1:nel
        n = enodes(e,:);
        edof(e,:) = [2*n(1)-1, 2*n(1), 2*n(2)-1, 2*n(2), 2*n(3)-1, 2*n(3), 2*n(4)-1, 2*n(4)];
    end
end