function [sig] = stressSolution(x,E,v)

Es = E/((1-2*v)*(1+v));

sig = zeros(3,1);

sig(1) = 2*pi*Es*((1-v)*cos(2*pi*x(1))*sin(2*pi*x(2))+v*sin(2*pi*x(1))*cos(2*pi*x(2)));
sig(2) = 2*pi*Es*(v*cos(2*pi*x(1))*sin(2*pi*x(2))+(1-v)*sin(2*pi*x(1))*cos(2*pi*x(2)));
sig(3) = pi*Es*(1-2*v)*(sin(2*pi*x(1))*cos(2*pi*x(2))+sin(2*pi*x(2))*cos(2*pi*x(1)));