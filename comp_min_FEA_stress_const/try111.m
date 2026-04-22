clc; clear; 

% theta=2000
% c=cosd(theta); 
% s=sind(theta); 
% 
% Tinv = [c^2,        s^2,     -2*c*s;
%         s^2,        c^2,     2*c*s;
%         c*s,        -c*s,    c^2-s^2];
% 
% T     = [c^2,        s^2,     2*c*s;
%         s^2,        c^2,     -2*c*s;
%         -c*s,       c*s,    c^2-s^2];
% T*Tinv

p=100; 
n=20; 
ge = rand(n,1);


gpn = ((1/n)*sum(ge.^p))^(1/p)
gmax=max(ge)