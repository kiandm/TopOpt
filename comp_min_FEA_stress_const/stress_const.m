function [gpn, gmax_con, ge, dgmax_x, dgmax_theta, cI] = stress_const(xphy, coords, conn, ...
          matprop, pl, p, s_pl, gs, gs1, numele, U, edofMat, K, numnode,freedofs,cI_1,alpha_I,gpn_I_1,ge_I_1,H,Hs)
    % Zahur's code
    % Design variables
    x = xphy(1:numele);
    theta = xphy(numele+1:end);

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

    % Material properties for Tsai-Wu criterion
    F1 = matprop.F1;
    F2 = matprop.F2;
    F11 = matprop.F11;
    F22 = matprop.F22;
    F12 = matprop.F12;
    F66 = matprop.F66;

    % Calculaiton of strress (S) and relaxed stress (S_hat) in material
    % coordiantes 12
    [s,s_hat]=stress(x,theta,s_pl,coords,conn,numele,U,edofMat,gs1,matprop);      

    % Initialize constraint values and derivatives in local coordinates
    ge = zeros(numele, 1);
    dge_s12 = zeros(3, numele);

    % Calculate local Tsai-Wu criterion and its derivatives
    for e = 1:numele
        s1 = s_hat(e, 1); s2 = s_hat(e, 2); s6 = s_hat(e, 3);
        
        % Local Tsai-Wu failure criterion
        ge(e) = F1*s1 + F2*s2 + F11*s1^2 + F22*s2^2 + 2*F12*s1*s2 + F66*s6^2;
        
        % Derivatives of g_e with respect to sigma1, sigma2, sigma6
        dge_s12(1,e) = F1 + 2*F11*s1 + 2*F12*s2; %∂g/∂σ1
        dge_s12(2,e) = F2 + 2*F22*s2 + 2*F12*s1;  
        dge_s12(3,e) = 2*F66*s6;
    end
    
    %filter ge before using
    % ge = H*(ge./Hs);

    % Global constraint (p-norm aggregation) ge<=1
    % gpn = (sum(ge.^p)/numele)^(1/p); 
    gpn = (sum(ge.^p))^(1/p); 

    % calculate gmax from gpn using CI
    cI = alpha_I*(max(ge_I_1)/gpn_I_1) + (1-alpha_I)*cI_1;
    gmax=cI*gpn;

    gmax_con = gmax-1; %ge<=1

    % if iter<=1
    %     CI= 1; 
    % else 
    %     CI = alpha_I*(max(ge_I_1)/gpn_I_1) + (1-alpha_I)*CI_1;
    % end 
    
    % Derivatives of the global constraint
    % dgpn_dge = ((sum(ge.^p)/numele)^(1/p-1)) * (1/numele)*(ge.^(p-1));
    dgpn_dge = (sum(ge.^p))^(1/p-1) * (ge.^(p-1));


    % Calculation of adjoint variable
    % Initialize right side of adjoint equation K * lambda = gama
    gama=zeros(2*numnode,1);
    for ee = 1:numele
        % Compute rotation matrix
        c = cos(theta(ee));
        s = sin(theta(ee));
        
        Tinv = [c^2,        s^2,     -2*c*s;
                s^2,        c^2,     2*c*s;
                c*s,        -c*s,    c^2-s^2];
        
        gg=gs1(:,ee);
        [phi, dphix, dphiy]=SF_FE(gg,coords,conn);

        Bmat=zeros(3,8);
        Bmat(1,1:2:end)=dphix;  Bmat(2,2:2:end)=dphiy;
        Bmat(3,1:2:end)=dphiy;  Bmat(3,2:2:end)=dphix;

        % Assemble the term inside the summation
        gama(edofMat(ee,:))=gama(edofMat(ee,:)) + ...
            dgpn_dge(ee) * x(ee)^s_pl * Bmat' * C12' * Tinv * dge_s12(:,ee);
    end
    LAMBDA=zeros(2*numnode,1);
    LAMBDA(freedofs,:)=K(freedofs,freedofs)\gama(freedofs,:);

    % Initialize sensitivities
    dgmax_x = zeros(numele, 1);
    dgmax_theta = zeros(numele, 1);
    
    % Loop over elements to calculate sensitivities
    gcount=0; 
    for ee=1:numele

        % Compute rotation matrix
        c = cos(theta(ee));
        s = sin(theta(ee));

        Tinv = [c^2,        s^2,     -2*c*s;
                s^2,        c^2,     2*c*s;
                c*s,        -c*s,    c^2-s^2];
     
       % Rotate the elastic stiffness matrix C12 from material to global
       % coordiantes x,y
       Cxy= Tinv * C12 * Tinv';
        
       % Derivative of Tinv with respect to theta
        dTinv_dtheta = [-sin(2*theta(ee)),  sin(2*theta(ee)), -2*cos(2*theta(ee));
                         sin(2*theta(ee)), -sin(2*theta(ee)),  2*cos(2*theta(ee));
                         cos(2*theta(ee)), -cos(2*theta(ee)), -2*sin(2*theta(ee))];        
        dCxy_dtheta=dTinv_dtheta * C12 * Tinv' + Tinv * C12 * dTinv_dtheta';

        KE_s=zeros(8,8); %element stiffness for solid material
        dKE_dtheta = zeros(8, 8); % Initialize element stiffness matrix

        for ii=1:4
            gcount=gcount+1; gg=gs(:,gcount); 
            weight=gg(6); jac=gg(7);  
            [phi, dphix, dphiy]=SF_FE(gg,coords,conn);
            Bmat=zeros(3,8);
            Bmat(1,1:2:end)=dphix;  Bmat(2,2:2:end)=dphiy;
            Bmat(3,1:2:end)=dphiy;  Bmat(3,2:2:end)=dphix;
            KE_s=KE_s+jac*weight*Bmat'*Cxy*Bmat;
            dKE_dtheta=dKE_dtheta + jac*weight*Bmat'*dCxy_dtheta*Bmat;
        end 

        ue=U(edofMat(ee,:));           
        lambda =LAMBDA(edofMat(ee,:)); 


        %calculate Bmat again at the centre of element to use in the Eqs
        gg=gs1(:,ee);
        [phi, dphix, dphiy]=SF_FE(gg,coords,conn);

        Bmat=zeros(3,8);
        Bmat(1,1:2:end)=dphix;  Bmat(2,2:2:end)=dphiy;
        Bmat(3,1:2:end)=dphiy;  Bmat(3,2:2:end)=dphix;


        % Sensitivity of g with respect to density x
        dgmax_x(ee) = cI*dgpn_dge(ee)*dge_s12(:,ee)'*s_pl*x(ee)^(s_pl-1)*C12*Tinv'*Bmat*ue ...
                     - cI*lambda'*pl*x(ee)^(pl-1)*KE_s*ue;       
        
        % Sensitivity of g with respect to fiber orientation theta
        dgmax_theta(ee) = cI*dgpn_dge(ee)*dge_s12(:,ee)'*x(ee)^s_pl*C12*dTinv_dtheta'*Bmat*ue ...
                       - cI*lambda'*x(ee)^pl*dKE_dtheta*ue;       

    end

end