function [coords, conn, edofMat, numnode, numele, freedofs, ...
                F, W]= problem_setup_portalframe(rmin) 

addpath('meshes/Portal_frame')
load elements.dat
conn=elements(:,2:5)';

load nodes.dat
coords=nodes(:,2:3)'; 

numnode=length(coords(1,:)); 
numele=length(conn(1,:));
edofMat=zeros(numele,8); 
edofMat(:,1:2:end)=2*conn'-1'; 
edofMat(:,2:2:end)=2*conn'; 

figure
patch('Faces',conn','Vertices',coords','LineWidth',1,...
      'FaceColor',[0,153/255,153/255],'EdgeColor','k'); 
hold on 
plot(coords(1,:),coords(2,:),'ok','markersize',4,'markerfacecolor',[1,0,0])
xlabel('x','fontsize',14); ylabel('y','fontsize',14);
axis equal 

%Find out the dirichlet dofs fully fixed 
bc_data;
plot(coords(1,bc_xy),coords(2,bc_xy),'ok','markersize',6,'markerfacecolor',[0,0,1])
plot(coords(1,bc_y),coords(2,bc_y),'ok','markersize',6,'markerfacecolor',[0,0,1])

%Find out the dirichlet dofs 
fix_dofs1=2*bc_xy-1; 
fix_dofs2=2*bc_xy; 
fix_dofs3=2*bc_y; 

fixeddofs=[fix_dofs1 fix_dofs2 fix_dofs3]; 

alldofs = 1:2*numnode;
freedofs = setdiff(alldofs,fixeddofs);

%Traction boundary (point load)
coords_trac = coords(:,nodes_trac);
plot(coords_trac(1,:),coords_trac(2,:),'ok','markersize',8,'markerfacecolor',[0,0.8,0])
xlabel('x','fontsize',14); ylabel('y','fontsize',14);


%Traction boundary (point load)
force = -300/6;   %2*10/nodes
F=zeros(2*numnode,1); 
F(2*nodes_trac,1)=force; 

% area=abs(coords(1,nodes_trac(2))-coords(1,nodes_trac(1))); 
% for el=1:length(nodes_trac)-1
%     F(2*nodes_trac(el))  =F(2*nodes_trac(el))  +pres*area/2;
%     F(2*nodes_trac(el+1))=F(2*nodes_trac(el+1))+pres*area/2;
% 
% end

%filter for density
%coordinate of background cells centres
gs=gauss_domain(coords,numele,conn,1); 
coords_cells(1,:) = gs(4,:); 
coords_cells(2,:) = gs(5,:); 

%to find out domain of influence of the node
xy=coords(:,conn(:,1)); 
r1=sqrt( (xy(1,1)-xy(1,3))^2  + (xy(2,1)-xy(2,3))^2  ); 
r2=sqrt( (xy(1,2)-xy(1,4))^2  + (xy(2,2)-xy(2,4))^2  ); 
r=max(r1, r2); 
xspac=r;
yspac=r;

dm_cells(1,1:numele)=rmin*(xspac*ones(1,numele));
dm_cells(2,1:numele)=rmin*(yspac*ones(1,numele));

W=zeros(numele,numele); 
for cc=1:numele
    gpos=coords_cells(:,cc); 
    [v,L]=nodes_in_support(numele, coords_cells, gpos, dm_cells);
    xi=coords_cells(:,v); 
    difx=abs((gpos(1,1)-xi(1,:))); 
    dify=abs((gpos(2,1)-xi(2,:))); 
    dif=sqrt(difx.^2 + dify.^2); 
    rij=dif./sqrt(dm_cells(1,v).^2 + dm_cells(2,v).^2);
    wij=(rmin-rij)./rmin; 
    W(cc,v)=wij; 
end
W=W./sum(W,2); 
W=sparse(W); 