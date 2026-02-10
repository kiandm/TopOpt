function [nodes, enodes, edof, nelx, nely, nn, ndof] = createMesh(test)
    % Mesh parameters
    nelx = test; % Number of elements in each direction
    nely = test;
    Nx = 2*nelx; % Number of half steps
    Ny = 2*nely;
    
    % Initialise a matrix of -1's at every possible node location
    idmap = -ones(Nx+1, Ny+1);

    % Generate all node coordinates and DOFs
    node = 0;
    for j = 0:Ny
        for i = 0:Nx
            % Skip middle nodes for Q8
            if mod(i,2)==1 && mod(j,2)==1                
                continue 
            end 
                node = node + 1;        % Keep track of number of nodes
                idmap(i+1,j+1) = node;  % Assign node number to idmap position
                nodes(node,1) = node;
                nodes(node,2) = i/Nx;          % x coordinate
                nodes(node,3) = j/Ny;          % y coordinate
                nodes(node,4) = 2*node - 1;      % Ux DOF
                nodes(node,5) = 2*node;          % Uy DOF
            
        end
    end
    nn = node;    % Number of nodes
    ndof = 2*nn;  % Number of DoFs

    % Element connectivity
    nel = nelx*nely;          % Number of elements
    enodes = zeros(nel, 8);   % Initialise table of elements with associated nodes
    elem = 0;
    for ey = 0:nely-1
        for ex = 0:nelx-1
            elem = elem + 1;
            i = 2*ex; 
            j = 2*ey;
            n1 = idmap(i+1, j+1);
            n2 = idmap(i+2+1, j+1);
            n3 = idmap(i+2+1, j+2+1);
            n4 = idmap(i+1, j+2+1);
            n5 = idmap(i+1+1,j+1);
            n6 = idmap(i+2+1, j+1+1);
            n7 = idmap(i+1+1,j+2+1);
            n8 = idmap(i+1,j+1+1);
            enodes(elem,:) = [n1, n2, n3, n4, n5, n6, n7, n8];
        end
    end

    % DOF connectivity
    edof = zeros(nel, 16);
    for e = 1:nel
        n = enodes(e,:);
        edof(e,:) = [2*n(1)-1, 2*n(1), 2*n(2)-1, 2*n(2), 2*n(3)-1, 2*n(3), 2*n(4)-1, 2*n(4), ...
                     2*n(5)-1, 2*n(5), 2*n(6)-1, 2*n(6), 2*n(7)-1, 2*n(7), 2*n(8)-1, 2*n(8)];
    end
end