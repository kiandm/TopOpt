function [nn,nodes,nel,enodes,ndof,edof] = createMesh()
% Imports list of node coordinates, as well as element construction from
% Cubit, then constructs connectivity.
addpath('C:\Users\sfxz32\OneDrive - Durham University\Desktop\Cubit')
nodes_file    = readmatrix("lbeam_nodes.dat");      % [id, x, y]
elements_file = readmatrix("lbeam_elements.dat");   % [eid, n1, n2, n3, n4]
% Node coordinates & DoFs
nn = size(nodes_file, 1);           % Number of nodes
nodes = zeros(nn, 5);               % [id, x, y, dofUx, dofUy]
nodes(:,1) = nodes_file(:,1);       % id
nodes(:,2) = nodes_file(:,2);       % x
nodes(:,3) = nodes_file(:,3);       % y
nodes(:,4) = 2*nodes(:,1) - 1;      % Ux DOF = 2*id - 1
nodes(:,5) = 2*nodes(:,1);          % Uy DOF = 2*id
% Element connectivity
nel = size(elements_file, 1);
enodes = elements_file(:, 2:5);     % take columns n1 n2 n3 n4
edof = zeros(nel, 8);
for e = 1:nel
    n = enodes(e,:);
    edof(e,:) = [2*n(1)-1, 2*n(1), ...
                 2*n(2)-1, 2*n(2), ...
                 2*n(3)-1, 2*n(3), ...
                 2*n(4)-1, 2*n(4)];
end
% Total DOFs
ndof = 2 * nn;
end