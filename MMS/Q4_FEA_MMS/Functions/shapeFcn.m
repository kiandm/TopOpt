function [N,dNdxi] = shapeFcn(xi,eta)
%SHAPE FUNCTIONS for 4-node quad element
%N = [ (1-xi)*(1-eta)/4, (1+xi)*(1-eta)/4, (1+xi)*(1+eta)/4, (1-xi)*(1+eta)/4 ];
N1 = (1-xi)*(1-eta)/4;
N2 = (1+xi)*(1-eta)/4;
N3 = (1+xi)*(1+eta)/4;
N4 = (1-xi)*(1+eta)/4;
N = [N1, N2, N3, N4];

% Derivatives of shape functions
% dN/dxi and dN/deta (rows: [xi; eta], cols: N1..N4)
% dNdxi = [ (eta-1)/4,  (1-eta)/4,  (eta+1)/4, -(eta+1)/4; (xi -1)/4, -(xi+1)/4,  (xi+1)/4,   (1 -xi)/4 ];
dN1dxi = (eta-1)/4;
dN2dxi = (1-eta)/4;
dN3dxi = (eta+1)/4;
dN4dxi = -(eta+1)/4;
dN1deta = (xi -1)/4;
dN2deta = -(xi+1)/4;
dN3deta = (xi+1)/4;
dN4deta = (1 -xi)/4;

dNdxi = [dN1dxi, dN2dxi, dN3dxi, dN4dxi;
         dN1deta, dN2deta, dN3deta, dN4deta];

end

