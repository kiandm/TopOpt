function [fb] = bodyForce(x,E,v)

Es = E/((1-2*v)*(1+v));

dsigxx_dx = 4*pi^2*Es*( -(1-v)*sin(2*pi*x(1))*sin(2*pi*x(2)) + v*cos(2*pi*x(1))*cos(2*pi*x(2)) );
dsigyy_dy = 4*pi^2*Es*( -(1-v)*sin(2*pi*x(1))*sin(2*pi*x(2)) + v*cos(2*pi*x(1))*cos(2*pi*x(2)) );
dsigxy_dy = 2*pi^2*(1-2*v)*Es*( -sin(2*pi*x(1))*sin(2*pi*x(2)) + cos(2*pi*x(1))*cos(2*pi*x(2)) );
dsigyx_dx = 2*pi^2*(1-2*v)*Es*( -sin(2*pi*x(1))*sin(2*pi*x(2)) + cos(2*pi*x(1))*cos(2*pi*x(2)) );

fb = -[dsigxx_dx + dsigxy_dy ;
       dsigyx_dx + dsigyy_dy];