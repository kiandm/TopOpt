function [etpl,coord]=formCoord3NT(nelsx,nelsy,lx,ly)

%Two dimensional finite element mesh generation
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   06/05/2015
% Description:
% Function to generate a 2D finite element grid of linear triangular
% elements.
%
%--------------------------------------------------------------------------
% [etpl,coord] = FORMCOORD3NT(nelsx,nelsy,lx,ly)
%--------------------------------------------------------------------------
% Input(s):
% nelsx - number of elements in the x direction
% nelsy - number of elements in the y direction
% lx    - length in the x direction
% ly    - length in the y direction
%--------------------------------------------------------------------------
% Ouput(s);
% etpl  - element topology
% coord - nodal coordinates
%--------------------------------------------------------------------------
% See also:
%
%--------------------------------------------------------------------------


nels=nelsx*nelsy;                                                           % half the number of elements

%% node generation
nodes=0;
coord=zeros(nels*4,2); 
for j=0:nelsy
  y=ly*j/nelsy;  
  for i=0:nelsx
    nodes=nodes+1;
    x=lx*i/nelsx;
    coord(nodes,:)=[x y];
  end
end
coord(nodes+1:end,:)=[];                                                    % remove unused nodes                                              

%% element generation
i=0; 
etpl=zeros(2*nels,3);
for nely=1:nelsy
  for nelx=1:nelsx
    i=i+1;
    etpl(i,1)=(nely-1)*(nelsx+1)+nelx;
    etpl(i,2)=nely*(nelsx+1)+nelx;  
    etpl(i,3)=etpl(i,1)+1; 
    
    i=i+1;
    etpl(i,1)=etpl(i-1,2)+1;
    etpl(i,2)=etpl(i-1,2);  
    etpl(i,3)=etpl(i-1,3); 
  end
end