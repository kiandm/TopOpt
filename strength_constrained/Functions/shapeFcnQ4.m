function [N, dNdxi] = shapeFcnQ4(xi, eta)
    % Bilinear Q4 shape functions and derivatives wrt (xi, eta)
    N = 0.25 * [ (1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta) ];
    dNdxi = 0.25 * [  -(1-eta),   (1-eta),   (1+eta),  -(1+eta);
                      -(1-xi),   -(1+xi),    (1+xi),    (1-xi) ];
end