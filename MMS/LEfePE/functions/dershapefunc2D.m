function [wp,dNr]=dershapefunc2D(ngp,nen)

%2D Basis function derivatives and Gauss point weights
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   30/09/2020
% Description:
% Function to determine derivatives of the basis functions (shape
% functions) with respect to the local coordinates and the weights of the 
% Gauss points that are used to integrate the finite elements. 
% 
% Elements included:
%   -  3 noded linear triangle
%   -  4 noded linear quadrilateral
%   -  8 noded quadratic quadrilateral
%   -  9 noded quadratic quadrilateral
%
%--------------------------------------------------------------------------
% [wp,dNr,gp] = DERSHAPEFUNC2D(ngp,nen)
%--------------------------------------------------------------------------
% Input(s):
% ngp   - number of Gauss points used to integrate the element
% nen   - number of nodes for each element
%--------------------------------------------------------------------------
% Ouput(s);
% wp    - Gauss point weights (ngp,1)
% dNr   - derivatives of the basis functions with respect to local position
%         evaluated at all Gauss point locations (nD*ngp,nen) 
%--------------------------------------------------------------------------
% See also:
% 
%--------------------------------------------------------------------------

%% Gauss point information
if ngp==1 && nen==3
  gp = [1 1]/4;
  wp = 1/2;
elseif ngp==1
  gp = [0 0];
  wp = 8;
elseif ngp==4
  g2=1/sqrt(3);
  gp(:,1)=[-1 -1 1 1].'*g2;
  gp(:,2)=[-1 1 1 -1].'*g2;
  wp=ones(4,1);
else
  g2=sqrt(3/5);
  gp(:,1)=[-1 -1 -1  0  0  0  1  1  1].'*g2;
  gp(:,2)=[-1  1  0 -1  1  0 -1  1  0].'*g2;
  wp1=[5/9 5/9 5/9 8/9 8/9 8/9 5/9 5/9 5/9].';
  wp2=[5/9 5/9 8/9 5/9 5/9 8/9 5/9 5/9 8/9].';
  wp=wp1.*wp2;  
end

%% Basis/shape function derivatives
xsi=gp(:,1); eta=gp(:,2); r2=ngp*2;
if nen==8
  dNr(1:2:r2,1)=1/4*(2*xsi+eta).*(1-eta);
  dNr(1:2:r2,2)=-1/2*(1-eta.^2);
  dNr(1:2:r2,3)=1/4*(2*xsi-eta).*(1+eta); 
  dNr(1:2:r2,4)=-xsi.*(1+eta);
  dNr(1:2:r2,5)=1/4*(1+eta).*(2*xsi+eta);
  dNr(1:2:r2,6)=1/2*(1-eta.^2);
  dNr(1:2:r2,7)=1/4*(1-eta).*(2*xsi-eta); 
  dNr(1:2:r2,8)=-xsi.*(1-eta);
  dNr(2:2:r2+1,1)=1/4*(1-xsi).*(xsi+2*eta);
  dNr(2:2:r2+1,2)=-eta.*(1-xsi);
  dNr(2:2:r2+1,3)=1/4*(1-xsi).*(2*eta-xsi); 
  dNr(2:2:r2+1,4)=1/2*(1-xsi.^2);
  dNr(2:2:r2+1,5)=1/4*(1+xsi).*(2*eta+xsi);
  dNr(2:2:r2+1,6)=-eta.*(1+xsi);
  dNr(2:2:r2+1,7)=1/4*(1+xsi).*(2*eta-xsi);
  dNr(2:2:r2+1,8)=-1/2*(1-xsi.^2);
elseif nen==9
  dNr(1:2:r2,1)= 1/4*eta.*(eta-1).*(2*xsi-1);
  dNr(1:2:r2,2)=-1/2*(eta+1).*(eta-1).*(2*xsi-1);
  dNr(1:2:r2,3)= 1/4*eta.*(eta+1).*(2*xsi-1);
  dNr(1:2:r2,4)=-1*eta.*(eta+1).*xsi;
  dNr(1:2:r2,5)= 1/4*eta.*(eta+1).*(2*xsi+1);
  dNr(1:2:r2,6)=-1/2*(eta+1).*(eta-1).*(2*xsi+1);
  dNr(1:2:r2,7)= 1/4*eta.*(eta-1).*(2*xsi+1);
  dNr(1:2:r2,8)=-1*eta.*(eta-1).*xsi;
  dNr(1:2:r2,9)= 2*(eta+1).*(eta-1).*xsi;
  
  dNr(2:2:r2+1,1)= 1/4*xsi.*(xsi-1).*(2*eta-1);
  dNr(2:2:r2+1,2)= -1*xsi.*(xsi-1).*eta;
  dNr(2:2:r2+1,3)= 1/4*xsi.*(xsi-1).*(2*eta+1);
  dNr(2:2:r2+1,4)=-1/2*(xsi+1).*(xsi-1).*(2*eta+1);
  dNr(2:2:r2+1,5)= 1/4*xsi.*(xsi+1).*(2*eta+1);
  dNr(2:2:r2+1,6)=-1*xsi.*(xsi+1).*eta;
  dNr(2:2:r2+1,7)= 1/4*xsi.*(xsi+1).*(2*eta-1);
  dNr(2:2:r2+1,8)=-1/2*(xsi+1).*(xsi-1).*(2*eta-1); 
  dNr(2:2:r2+1,9)= 2*(xsi+1).*(xsi-1).*eta;
elseif nen==4
  dNr(1:2:r2,1)  =-1/4*(1-eta);
  dNr(1:2:r2,2)  =-1/4*(1+eta);
  dNr(1:2:r2,3)  = 1/4*(1+eta);
  dNr(1:2:r2,4)  = 1/4*(1-eta);
  dNr(2:2:r2+1,1)=-1/4*(1-xsi);
  dNr(2:2:r2+1,2)= 1/4*(1-xsi);
  dNr(2:2:r2+1,3)= 1/4*(1+xsi);
  dNr(2:2:r2+1,4)=-1/4*(1+xsi);
elseif nen==3
  dNr(1:2:r2,1)  =-1;
  dNr(1:2:r2,2)  = 1;
  dNr(1:2:r2,3)  = 0;
  dNr(2:2:r2+1,1)=-1;
  dNr(2:2:r2+1,2)= 0;
  dNr(2:2:r2+1,3)= 1;
end
end