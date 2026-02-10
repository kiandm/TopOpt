%LEfe2dpe: 2D plane strain Linear Elastic finite element solver
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   30/09/2020
% Description:
% Linear Elastic finite element (LEfe) solver for two dimensional analysis
% with a plane strain assumption in the third direction. 
%
%--------------------------------------------------------------------------
% LEFE2DPE
%--------------------------------------------------------------------------
% Input(s):
% 
%--------------------------------------------------------------------------
% Ouput(s);
%  
%--------------------------------------------------------------------------
% See also:
% SETUPMESH         - analysis mesh and material information
% DERSHAPEFUNC2D    - 2D basis function derivatives
%--------------------------------------------------------------------------

clear; addpath('functions'); tic;
[coord,etpl,fext,bc,ngp,E,v]=setupmesh;                                     % problem set up information
[nodes,nD]=size(coord);                                                     % number of nodes and dimensions
[nels,nen]=size(etpl);                                                      % number of ekements and nodes per element
nDoF=nodes*nD;                                                              % number of degrees of freedom
neDoF=(nen*nD)^2;                                                           % number of entries in the element stiffness matrix
krow=zeros(neDoF*nels,1); kcol=krow; kval=krow;                             % zero global stiffness storage
B=zeros(3,nD*nen);                                                          % zero the strain-displacement matrix
[wp,dNr]=dershapefunc2D(ngp,nen);                                           % Gauss point weights and shape function derivatives
bm1=[1 1 0].'; I6=[1 0 0; 0 1 0; 0 0 0.5]; 
D=E/((1+v)*(1-2*v))*((1-2*v)*I6+v*(bm1*bm1.'));                             % elastic stiffness matrix
for nel=1:nels                                                              % START of element loop
  ed=ones(nD,1)*(etpl(nel,:)-1)*nD+(1:nD).'*ones(1,nen);                    
  ed=reshape(ed,1,nen*nD);                                                  % element degrees of freedom 
  JT=dNr*coord(etpl(nel,:),:);                                              % Jacobian matrix (all Gauss points)
  ke=zeros(nen*nD);                                                         % zero element stiffness matrix
  fe = zeros(nen*nD,1);
  for gp=1:ngp                                                              % START of Gauss point loop
    indx=nD*gp-(nD-1:-1:0);                                                 % matrix location of information assocuated with current gp
    detJ=det(JT(indx,:));                                                   % determinant of the Jacobian
    dNx=(JT(indx,:))\dNr(indx,:);                                           % global shape function derivatives
    B([1 3],1:2:end)=dNx(:,:);                                              % strain-displacement matrix 
    B([3 2],2:2:end)=dNx(:,:); 
    ke=ke+(B.'*D*B*detJ*wp(gp));                                            % element stiffness matrix
    
    N = shapefunc(gp,ngp,nen);
    
    NN = zeros(nD,nen*nD);
    NN(1,1:nD:end) = N;
    NN(2,2:nD:end) = N;
    
    eC = coord(etpl(nel,:),:);
    x  = N.'*eC;
    
    fb = bodyForce(x,E,v);
    
    fe=fe+NN.'*fb*detJ*wp(gp);
    
  end                                                                       % END of Gauss point loop
  krow((nel-1)*neDoF+1:nel*neDoF)=reshape(ed.'*ones(1,nen*nD),neDoF,1);     % stiffness storage: row locations
  kcol((nel-1)*neDoF+1:nel*neDoF)=reshape(ones(nen*nD,1)*ed  ,neDoF,1);     % stiffness storage: column locations
  kval((nel-1)*neDoF+1:nel*neDoF)=reshape(ke,neDoF,1);                      % stiffness storgae: stiffness values
  fext(ed) = fext(ed) + fe;
end                                                                         % END of element loop
K=sparse(krow,kcol,kval,nDoF,nDoF);                                         % form global stiffness matrix
c=toc; fprintf('formation time = %8.4f s\n',c); tic;                        % output assembly time

d=zeros(nDoF,1);                                                            % zero displacements 
d(bc(:,1))=bc(:,2);                                                         % imposed displacement BCs
fd=(1:nDoF); fd(bc(:,1))=[];                                                % free degrees of freedom
d(fd)=K(fd,fd)\(fext(fd)-K(fd,bc(:,1))*bc(:,2));                            % solve f = K*d for d 
c=toc; fprintf('solution time  = %8.4f s\n',c); tic;                        % output solver time

sigErr = 0;
disErr = 0;
sigma=zeros(3,ngp,nels);                                                    % initialisation of the sigma array
for nel=1:nels                                                              % START of element loop
  ed=ones(nD,1)*(etpl(nel,:)-1)*nD+(1:nD).'*ones(1,nen);                    
  ed=reshape(ed,1,nen*nD);                                                  % element degrees of freedom 
  JT=dNr*coord(etpl(nel,:),:);                                              % Jacobian matrix (all Gauss points)                       
  for gp=1:ngp                                                              % START of Gauss point loop
    indx=nD*(gp-1)+(1:nD);                                                  % index for current Gp information
    dNx=JT(indx,:)\dNr(indx,:);                                             % global shape func. derivatives                            
    B([1 3],1:2:end)=dNx(:,:);                                              % strain-displacement matrix 
    B([3 2],2:2:end)=dNx(:,:);
    strain=B*d(ed);                                                         % Strain                                                    
    sigma(:,gp,nel)=D*strain;                                               % Cauchy stress    
    
    N  = shapefunc(gp,ngp,nen);
    eC = coord(etpl(nel,:),:);
    x  = N.'*eC;
    sig = stressSolution(x,E,v);
    
    detJ = det(JT(indx,:));
    
    sigErr = sigErr+norm(sigma(:,gp,nel)-sig)/norm(sig)*detJ*wp(gp);
    
    NN = zeros(nD,nen*nD);
    NN(1,1:nD:end) = N;
    NN(2,2:nD:end) = N;
    u = dispSolution(x);
    disErr = disErr + norm(NN*d(ed)-u)*detJ*wp(gp);
    
  end                                                                       % END of Gauss point loop
end                                                                         % END of element loop

uErr = 0;
for node = 1:nodes
    u = dispSolution(coord(node,:));
    uErr = uErr + norm(d([node*nD-1 node*nD].')-u)/nodes;
end

c=toc; fprintf('post pro time  = %8.4f s\n',c);                             % output post-processing time

% makeVtk(coord,etpl,d,fext,'mesh.vtk')