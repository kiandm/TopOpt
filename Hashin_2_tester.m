% Hashin Fibre Angle Comparison
clearvars;
clc;
close all;

% Inputs

L = 200;              % Plate length (mm)
H = 100;              % Plate height (mm)
thickness = 1;        % Plate thickness (mm)

nx = 20;              % Elements along length
ny = 10;              % Elements along height

force = -1000;        % Downward force (N)
penal = 3;            % SIMP penalty
initial_density = 1;  % Full material

% Test angles 
angles = -90:5:90;

% Elastic properties (MPa)
matprop.E1 = 135000;
matprop.E2 = 10000;
matprop.nu12 = 0.30;
matprop.nu21 = matprop.nu12*matprop.E2/matprop.E1;
matprop.G12 = 5000;

% Strength properties (MPa)
strength.Xt = 1500;
strength.Xc = 1200;
strength.Yt = 40;
strength.Yc = 200;
strength.S = 70;

%Create Mesh

numele = nx*ny;

[X,Y] = meshgrid(linspace(0,L,nx+1), ...
                 linspace(0,H,ny+1));

coords = [X(:),Y(:)];
nnode = size(coords,1);
ndofs = 2*nnode;

conn = zeros(numele,4);
edofMat = zeros(numele,8);
centres = zeros(numele,2);

e = 0;

for ix = 1:nx
    for iy = 1:ny

        e = e+1;

        n1 = (ix-1)*(ny+1)+iy;
        n2 = ix*(ny+1)+iy;
        n3 = n2+1;
        n4 = n1+1;

        nodes = [n1 n2 n3 n4];

        conn(e,:) = nodes;

        edofMat(e,:) = reshape( ...
            [2*nodes-1;2*nodes],1,[]);

        centres(e,:) = mean(coords(nodes,:),1);

    end
end

% Define Supports and Load

% Fix the left edge
fixednodes = find(abs(coords(:,1)) < 1e-9);

fixeddofs = sort([2*fixednodes-1; ...
                  2*fixednodes]);

fixeddofs = fixeddofs(:);

freedofs = setdiff((1:ndofs)',fixeddofs);

% Load at top-right corner
[~,loadnode] = min( ...
    (coords(:,1)-L).^2 + ...
    (coords(:,2)-H).^2);

F = zeros(ndofs,1);
F(2*loadnode) = force;

% Initialise Results

results = zeros(length(angles),6);

modes = {'Fibre tension', ...
         'Fibre compression', ...
         'Matrix tension', ...
         'Matrix compression'};

% Loop through Fibre Angles

for a = 1:length(angles)

    fibre_angle = angles(a);

    fprintf('\nTesting fibre angle: %g degrees\n', ...
            fibre_angle);

    %Set Density and Fibre Orientation

    density = initial_density*ones(numele,1);

    theta = deg2rad(fibre_angle)*ones(numele,1);

    xphy = [density;theta];

    %Material Stiffness Matrix

    E1 = matprop.E1;
    E2 = matprop.E2;
    nu12 = matprop.nu12;
    nu21 = matprop.nu21;
    G12 = matprop.G12;

    C0 = [E1, nu21*E1, 0;
          nu12*E2, E2, 0;
          0, 0, G12] / (1-nu12*nu21);

    C0(3,3) = G12;

  

    gpts = [-1 1]/sqrt(3);

    gs = zeros(7,4*numele);
    dphix_ref = zeros(4,4);
    dphiy_ref = zeros(4,4);

    KE0 = zeros(8,8,numele);

    rows = zeros(64*numele,1);
    cols = zeros(64*numele,1);
    vals = zeros(64*numele,1);

    counter = 0;

    for e = 1:numele

        nodes = conn(e,:);
        xy = coords(nodes,:);

        theta_e = theta(e);
        c = cos(theta_e);
        s = sin(theta_e);

        Tinv = [c^2,s^2,-2*c*s;
                s^2,c^2,2*c*s;
                c*s,-c*s,c^2-s^2];

        Cxy = Tinv*C0*Tinv';

        Ke = zeros(8,8);
        gp = 0;

        for i = 1:2
            for j = 1:2

                gp = gp+1;

                xi = gpts(i);
                eta = gpts(j);

                dN = 0.25 * ...
                    [-(1-eta), -(1-xi);
                       (1-eta), -(1+xi);
                       (1+eta),  (1+xi);
                      -(1+eta),  (1-xi)];

                J = xy'*dN;
                detJ = det(J);

                dNxy = dN/J;

                dphix_ref(:,gp) = dNxy(:,1);
                dphiy_ref(:,gp) = dNxy(:,2);

                B = zeros(3,8);

                B(1,1:2:8) = dNxy(:,1)';
                B(2,2:2:8) = dNxy(:,2)';
                B(3,1:2:8) = dNxy(:,2)';
                B(3,2:2:8) = dNxy(:,1)';

                idx = (e-1)*4+gp;

                gs(6,idx) = 1;
                gs(7,idx) = detJ*thickness;

                Ke = Ke+B'*Cxy*B*detJ*thickness;

            end
        end

        Ke = density(e)^penal*Ke;
        KE0(:,:,e) = Ke;

        dofs = edofMat(e,:);

        [rr,cc] = ndgrid(dofs,dofs);
        ids = counter+(1:64);

        rows(ids) = rr(:);
        cols(ids) = cc(:);
        vals(ids) = Ke(:);

        counter = counter+64;

    end

    K = sparse(rows,cols,vals,ndofs,ndofs);

    % Solve Displacements

    U = zeros(ndofs,1);

    dK = K(freedofs,freedofs);

    U(freedofs) = dK\F(freedofs);

    %Use Hashin Function

    [g, dg_dx, dg_dtheta, FailIdx, vonMises] = ...
        Hashin_2(U,dK,KE0,xphy,penal,numele, ...
        gs,edofMat,coords,conn, ...
        matprop,strength,freedofs, ...
        dphix_ref,dphiy_ref);

    % Calculate Total Displacement

    Uxy = reshape(U,2,[])';

    displacement_magnitude = sqrt(sum(Uxy.^2,2));

    maxDisp = max(displacement_magnitude);

    % Store Results

    results(a,:) = [fibre_angle, ...
                    g(:)', ...
                    maxDisp];

    % Individual Failure Figures

    % figure('Name',sprintf('Fibre angle %g',fibre_angle), ...
    %        'NumberTitle','off');
    %
    % for k = 1:4
    %
    %     subplot(2,2,k);
    %
    %     scatter(centres(:,1),centres(:,2), ...
    %             70,FailIdx(:,k),'filled','s');
    %
    %     axis equal tight;
    %     colorbar;
    %     clim([0 1]);
    %
    %     title(modes{k});
    %     xlabel('X (mm)');
    %     ylabel('Y (mm)');
    %
    % end
    %
    % sgtitle(sprintf('Fibre angle = %g degrees', ...
    %                 fibre_angle));
    %
    % drawnow;

end

    %Create Results Table

T = array2table(results, ...
    'VariableNames', ...
    {'Angle','FibreTension','FibreCompression', ...
     'MatrixTension','MatrixCompression', ...
     'MaxDisplacement'});

disp(' ');
disp('--- FIBRE ANGLE COMPARISON ---');
disp(T);

% Hashin Failure Constrains

figure('Name','Hashin failure comparison', ...
       'NumberTitle','off');

plot(angles,results(:,2:5),'o-', ...
     'LineWidth',1.5, ...
     'MarkerSize',4);

hold on;

yline(0,'r--','Failure threshold');

grid on;

xlabel('Fibre angle (degrees)');
ylabel('Hashin constraint g');

legend(modes,'Location','best');

title('Hashin failure constraints vs fibre orientation');

% Maximum Displacement

figure('Name','Displacement comparison', ...
       'NumberTitle','off');

plot(angles,results(:,6),'o-', ...
     'LineWidth',1.5, ...
     'MarkerSize',4);

grid on;

xlabel('Fibre angle (degrees)');
ylabel('Maximum displacement magnitude (mm)');

title('Displacement vs fibre orientation');

% Max Hashin Constraint

worst_failure = max(results(:,2:5),[],2);

figure('Name','Worst Hashin failure constraint', ...
       'NumberTitle','off');

plot(angles,worst_failure,'k-o', ...
     'LineWidth',1.5, ...
     'MarkerSize',4);

hold on;

yline(0,'r--','Failure threshold');

grid on;

xlabel('Fibre angle (degrees)');
ylabel('Maximum Hashin constraint g');

title('Critical Hashin failure constraint vs fibre orientation');

% Identify Best Angles

[~,idx_disp] = min(results(:,6));

[~,idx_failure] = min(worst_failure);

fprintf('\n======================================\n');
fprintf('SUMMARY OF RESULTS\n');
fprintf('======================================\n');

fprintf('Angle giving minimum displacement: %.1f degrees\n', ...
        angles(idx_disp));

fprintf('Minimum displacement: %.4f mm\n', ...
        results(idx_disp,6));

fprintf('\nAngle giving lowest worst Hashin constraint: %.1f degrees\n', ...
        angles(idx_failure));

fprintf('Lowest worst Hashin constraint: %.4f\n', ...
        worst_failure(idx_failure));

fprintf('======================================\n');
