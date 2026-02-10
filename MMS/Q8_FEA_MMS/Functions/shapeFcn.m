function [N,dNdxi] = shapeFcn(xi,eta)
%SHAPE FUNCTIONS for 4-node quad element
%N = [ (1-xi)*(1-eta)/4, (1+xi)*(1-eta)/4, (1+xi)*(1+eta)/4, (1-xi)*(1+eta)/4 ];
N1 = (1/4) * (1 - xi) * (eta - 1) * (xi + eta + 1);
N2 = (1/4) * (1 + xi) * (eta - 1) * (-xi + eta + 1);
N3 = (1/4) * (1 + xi) * (eta + 1) * (xi + eta - 1);
N4 = (1/4) * (xi - 1) * (eta + 1) * (xi - eta + 1);
N5 = (1/2) * (1 - eta) * (1 - xi^2);
N6 = (1/2) * (1 + xi) * (1 - eta^2);
N7 = (1/2) * (1 + eta) * (1 - xi^2);
N8 = (1/2) * (1 - xi) * (1 - eta^2);
N = [N1, N2, N3, N4, N5, N6, N7, N8];

% Derivatives of shape functions
% dN/dxi and dN/deta (rows: [xi; eta], cols: N1..N4)
% dNdxi = [ (eta-1)/4,  (1-eta)/4,  (eta+1)/4, -(eta+1)/4; (xi -1)/4, -(xi+1)/4,  (xi+1)/4,   (1 -xi)/4 ];
dN1dxi = (1/4) * (eta - eta^2 + 2*xi - 2*xi*eta);
dN2dxi = (1/4) * (-eta + eta^2 + 2*xi - 2*xi*eta);
dN3dxi = (1/4) * (eta + eta^2 + 2*xi + 2*xi*eta);
dN4dxi = (1/4) * (-eta - eta^2 + 2*xi + 2*xi*eta);
dN5dxi = xi * (eta - 1);
dN6dxi = (1/2) * (1 - eta^2);
dN7dxi = -xi * (1 + eta);
dN8dxi = -(1/2) * (1 - eta^2);

dN1deta = (1/4) * (xi - xi^2 + 2*eta - 2*xi*eta);
dN2deta = (1/4) * (-xi - xi^2 + 2*eta + 2*xi*eta);
dN3deta = (1/4) * (xi + xi^2 + 2*eta + 2*xi*eta);
dN4deta = (1/4) * (-xi + xi^2 + 2*eta - 2*xi*eta);
dN5deta = (1/2) * (xi^2 - 1);
dN6deta = -eta * (1 + xi);
dN7deta = (1/2) * (1 - xi^2);
dN8deta = eta * (xi - 1);

dNdxi = [dN1dxi, dN2dxi, dN3dxi, dN4dxi, dN5dxi, dN6dxi, dN7dxi, dN8dxi;
         dN1deta, dN2deta, dN3deta, dN4deta, dN5deta, dN6deta, dN7deta, dN8deta];

end

