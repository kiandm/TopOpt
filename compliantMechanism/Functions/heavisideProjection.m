function [x,beta,eta] = heavisideProjection(x_phy)
%HEAVISIDEPROJECTION Summary of this function goes here
%   Detailed explanation goes here
x_phy = (tanh(beta * eta) + tanh(beta * (x - eta))) / (tanh(beta * eta) + tanh(beta * (1 - eta)));

end