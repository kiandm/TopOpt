% compliance minimisation x and theta as design variables and volume and
% p-norm stress constraints 
clear; clc; close all; warning off

% Parameters
volfrac = 0.4; % Volume fraction
pl = 3.0; % Penalization factor
rmin = 1.5; % Filter radius

% Material properties composites
matprop.E1=10; % Young's modulus in fiber direction
matprop.E2=1; % Young's modulus perpendicular to fiber direction
matprop.nu12=0.3; % Major Poisson's ratio
matprop.nu21=matprop.nu12*matprop.E2/matprop.E1; % Minor Poisson's ratio
matprop.G12=matprop.E1/(2*(1+matprop.nu12)); % Shear modulus


[coords, conn, edofMat, numnode, numele, freedofs, F, H]= problem_setup_Lbrac(rmin);
% [coords, conn, edofMat, numnode, numele, freedofs, F, H]= problem_setup_Lbrac(rmin);
Hs = sum(H,2);

U = zeros(2*numnode,1);
gs=gauss_domain(coords,numele,conn,2);
% gs1=gauss_domain(coords,numele,conn,1);  %1 gauss pint at the centre


% volume of each element
ve=zeros(numele,1);
gcount=0;
for i=1:numele 
    for gg=1:2*2 %for 2x2 gauss points
        gcount=gcount+1;
        weight=gs(6,gcount); jac=gs(7,gcount);
        ve(i)=ve(i) + jac*weight;
    end
end

% Store in global variables to be used in myOutputFcn
global coords_global conn_global ve_global;
coords_global = coords;
conn_global = conn;
ve_global = ve; 


% Initialize design variables
dens = volfrac * ones(numele,1); % Density variables
theta = pi/2 * ones(numele,1); % Fiber direction variables

% Combine design variables
x0 = [dens; theta];

% Bounds for densities and fiber directions
xmin = 1e-4 * ones(numele,1); % Lower bound for densities
xmax = 1 * ones(numele,1); % Upper bound for densities

xmin_theta =  -pi * ones(numele,1); % Lower bound for fiber directions
xmax_theta =   pi * ones(numele,1); % Upper bound for fiber directions

% Combine bounds
LB = [xmin; xmin_theta];
UB = [xmax; xmax_theta];

% Set options for fmincon                           
options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
    'MaxIter',100,...
    'TolX',1e-12, ...
    'TolFun', 1e-12, ... 
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'OutputFcn',@(x,optimValues,state) myOutputFcn(x,optimValues,state));

% [xopt, fopt]=fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
[x,fval]=fmincon(@(x) comp_dc(x, pl, numele, numnode, gs, edofMat, coords, conn, matprop, freedofs, F, H, Hs), ...
                 x0,[],[],[],[],LB,UB, ...
                 @(x) confun(x, volfrac, numele, ve, H, Hs), options); 

%plot angles
figure
patch('Faces',conn','Vertices',coords','FaceVertexCData',x(1:numele),...
      'FaceColor','flat','EdgeColor',[0 0 0]); colorbar; %axis equal off;
axis equal; hold on 
len=0.05; 
theta=x(numele+1:end); 
theta = H*(theta./Hs);

u_len=len*cos(theta); 
v_len=len*sin(theta); 

x_cen=(coords(1,conn(1,:))+coords(1,conn(4,:)))/2;    
y_cen=(coords(2,conn(1,:))+coords(2,conn(2,:)))/2;    

ind=find(x(1:numele)>0.5);
quiver(x_cen(ind)',y_cen(ind)',u_len(ind),v_len(ind),'Color','r','LineWidth',2); 




% Define the custom output function
function stop = myOutputFcn(x,optimValues,state)
    % this will make sure to only record when it is actual iteration
    % otherwise it will get values multiple times for each iteration    
    if strcmp(state, 'iter') 
        %PRINT RESULTS 
        loop=optimValues.iteration;  c=optimValues.fval; 
        change=optimValues.stepsize;
        
        global coords_global conn_global ve_global;
        coords = coords_global;
        conn = conn_global;
        ve = ve_global;

        numele = length(conn(1,:));
        vol_frac=sum(x(1:numele).*ve)/(sum(ve)); 

        disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
        ' Vol.: ' sprintf('%6.3f',vol_frac) ...
        ' ch.: ' sprintf('%6.3f',change)]) 

        % PLOT DENSITIES  
        % Store U and K in global variables

        colormap(jet);
        patch('Faces',conn','Vertices',coords','FaceVertexCData',x(1:numele),...
              'FaceColor','flat','EdgeColor',[0,0,0]); axis equal off; colorbar;  
        drawnow; clim([0 1]); 
    
    end % switch

    % Return false to continue optimization
    stop = false;
end % myOutputFcn
