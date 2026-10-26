function [x_phy,dxphy] = heavisideProjection(x,beta,eta)
% Smoothed Heaviside projection of the filtered density x (threshold eta, sharpness beta).
% beta = 0 means NO projection: the limit of the formula for beta -> 0 is x_phy = x and dxphy = 1,
% while the formula itself is 0/0 at beta = 0.
if beta == 0
    x_phy = x;
    dxphy = ones(size(x));
    return
end
x_phy = (tanh(beta * eta) + tanh(beta * (x - eta))) / ...
    (tanh(beta * eta) + tanh(beta * (1 - eta)));
dxphy = beta * (1 - tanh(beta * (x - eta)).^2) / ...
        (tanh(beta * eta) + tanh(beta * (1 - eta)));
end