function [u] = dispSolution(x)

u = [sin(2*pi*x(1))*sin(2*pi*x(2)) ;
     sin(2*pi*x(1))*sin(2*pi*x(2))];