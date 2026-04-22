function [S,Shat]=stress(x,theta,s_penl,coords,conn,numele,U,edofMat,gs1,matprop) 
            
% Material properties for the composite
E1 = matprop.E1; % Young's modulus in fiber direction
E2 = matprop.E2; % Young's modulus perpendicular to fiber direction
nu12 = matprop.nu12; % Major Poisson's ratio
nu21 = matprop.nu21; % Minor Poisson's ratio
G12 = matprop.G12; % Shear modulus

% Compute the reduced stiffness matrix Q
C12 = [E1/(1 - nu12 * nu21), nu12 * E2 / (1 - nu12 * nu21), 0;
     nu21 * E1 / (1 - nu12 * nu21), E2 / (1 - nu12 * nu21), 0;
     0, 0, G12];


S=zeros(numele,3); Shat=S; 
for ee=1:numele
    c = cos(theta(ee));
    s = sin(theta(ee));
    
    Tinv = [c^2,        s^2,     -2*c*s;
            s^2,        c^2,     2*c*s;
            c*s,        -c*s,    c^2-s^2];

    T     = [c^2,        s^2,     2*c*s;
            s^2,        c^2,     -2*c*s;
            -c*s,       c*s,    c^2-s^2];

    % Rotate the elastic stiffness matrix C12 from material to global
    % coordiantes x,y
    Cxy= Tinv * C12 * Tinv';
        

    gg=gs1(:,ee);
    [phi, dphix, dphiy]=SF_FE(gg,coords,conn);    
    Bmat=zeros(3,8);
    Bmat(1,1:2:end)=dphix;  Bmat(2,2:2:end)=dphiy;
    Bmat(3,1:2:end)=dphiy;  Bmat(3,2:2:end)=dphix;
    
    u=U(edofMat(ee,:)); 
    
    S(ee,:)=T*Cxy*Bmat*u; %S is original stress for solid matearil in material coordinates
    Shat(ee,:)=x(ee)^s_penl*S(ee,:); %Shat is relaxed stress
end
