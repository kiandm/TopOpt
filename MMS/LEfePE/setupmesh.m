function [coord,etpl,fext,bc,ngp,E,v]=setupmesh

%Mesh and material information
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   30/09/2020
% Description:
% Function to provide the finite element mesh information (such as nodal
% positions and mesh connectivity), the boundary conditions on the mesh and
% the properties of the analysed material. 
%
%--------------------------------------------------------------------------
% [coord,etpl,fext,bc,E,v,ngp,dens,sf] = SETUPMESH
%--------------------------------------------------------------------------
% Input(s):
% 
%--------------------------------------------------------------------------
% Ouput(s);
% coord   - nodal coordiantes (nodes,nD)
% etpl    - element connectivity (nels,nen)
% fext    - external force vector (nodes*nD,1)
% bc      - boundary condition matrix (*,2)
% ngp     - number of Gauss points per element
% E       - Young's modulus 
% v       - Poission's ratio
%--------------------------------------------------------------------------
% See also:
% FORMCOORD4NQ  - 4-noded element mesh generation
% FORMCOORD8NQ  - 8-noded element mesh generation
%--------------------------------------------------------------------------


lx    = 1;                                                                  % mesh length in the x direction
ly    = 1;                                                                  % mesh length in the y direction
nelsx = 2^9;                                                                % number of elements in the x direction
nelsy = nelsx;                                                              % number of elements in the y direction


E       = 1;                                                                % Young's modulus
v       = 0.2;                                                              % Poisson's ratio
ngp     = 9;                                                                % number of Gauss points

% [etpl,coord] = formCoord3NT(nelsx,nelsy,lx,ly);                             % mesh generation (LINEAR TRIANGLE)
% [etpl,coord] = formCoord4NQ(nelsx,nelsy,lx,ly);                             % mesh generation (LINEAR QUAD)
[etpl,coord] = formCoord8NQ(nelsx,nelsy,lx,ly);                             % mesh generation (QUADRATIC QUAD)

[nodes,nD] = size(coord);                                                   % number of nodes and dimensions

tol=1e-6; bc=zeros(nodes*nD,2);                                             % boundary conditions
for i=1:nodes
  if coord(i,1)<tol || abs(coord(i,1)-lx)<tol || coord(i,2)<tol || abs(coord(i,2)-ly)<tol 
      u = dispSolution(coord(i,:));
      bc(i*2-1,:)=[i*2-1 u(1)];
      bc(i*2  ,:)=[i*2   u(2)];
  end 
end
bc=sortrows(bc,1);                                                     
bc=bc(bc(:,1)>0,:);                                                         % remove excess rows from bc

fext=zeros(nodes*nD,1);                                                     % external force vector

end