function [etopol,coord]=formCoord8NQ(nelsx,nelsy,lx,ly)

%Two dimensional finite element grid generation
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   06/05/2015
% Description:
% Function to generate a 2D finite element grid of quadratic quadrilateral 
% elements.
%
%--------------------------------------------------------------------------
% [etpl,coord] = FORMCOORD8NQ(nelsx,nelsy,lx,ly)
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

nels=nelsx*nelsy;
nnodesxO = (2*nelsx)+1;
nnodesxE =   nelsx +1;
nodes    = nnodesxO*(nelsy+1) + nnodesxE*nelsy;
node=0; elem=0;

coord = zeros(nodes,2);
for ylayr=1:(2*nelsy+1)
  y=(ylayr-1)*ly/(2*nelsy);
  for xlayr=1:(2*nelsx+1)
    x=(xlayr-1)*lx/(2*nelsx);
    if(oddeve(ylayr)==1)
      node=node+1;
      coord(node,1:2)=[x y];
    elseif((((oddeve(ylayr)==0))||((oddeve(ylayr)==1)))&&(oddeve(xlayr)==1))
      node=node+1;
      coord(node,1:2)=[x y];
    end
  end
end
etopol=zeros(nels,8);
for elemy=1:nelsy
  for elemx=1:nelsx
      elem=elem+1;
      n1=((elemy-1)*(nnodesxO+nnodesxE))+(2*(elemx-1))+1;
      n2=(nnodesxO)*elemy+(nnodesxE)*(elemy-1)+elemx;
      n3=(nnodesxO+nnodesxE)*elemy+(elemx-1)*2+1;
      n4=n3+1;
      n5=n4+1;
      n6=n2+1;
      n7=n1+2;
      n8=n1+1;
      etopol(elem,1:8)=[n1 n2 n3 n4 n5 n6 n7 n8];
  end
end
end

function [oddeven]=oddeve(number)
if ((number/2)~=(floor(number/2)))
  oddeven=1; % odd
else
  oddeven=0; % even
end
end