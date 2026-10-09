function [x_phy,dxphy] = heavisideProjection(x,beta,eta)
%HEAVISIDE PROJECTION Summary of this function goes here
%   Detailed explanation goes here
x_phy = (tanh(beta * eta) + tanh(beta * (x - eta))) / ...
    (tanh(beta * eta) + tanh(beta * (1 - eta)));
dxphy = beta * (1 - tanh(beta * (x - eta)).^2) / ...
        (tanh(beta * eta) + tanh(beta * (1 - eta)));
end
