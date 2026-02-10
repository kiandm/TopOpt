function [coords, conn, edofMat, numnode, numele, freedofs, F, W]= problem_setup_mbb(rmin)
nely=20;  nelx=60; 
Hb=nely;  Lb=nelx; 

[x2d,y2d]=meshgrid(0:Lb/nelx:Lb, Hb:-Hb/nely:0);

size2d=size(x2d);
count=1;
numnode=size2d(1)*size2d(2);
numele=nelx*nely;

for i=1:size2d(2)
    for j=1:size2d(1)
        coords(1,count)= x2d(j,i);
        coords(2,count)= y2d(j,i);
        node_num(j,i)=count;
        count=count+1;
    end
end

%to find the nodes for each cell
count1=1;
count2=1;

for i=1:nelx
    for j=1:nely
        conn(1,count1)=node_num(count2);
        conn(2,count1)=node_num(count2+1);
        conn(3,count1)=node_num(count2+1+(nely+1));
        conn(4,count1)=node_num(count2+(nely+1));
        count1=count1+1;
        count2=count2+1;
    end
    count2=count2+1;
end

edofMat=zeros(numele,8); 
edofMat(:,1:2:end)=2*conn'-1'; 
edofMat(:,2:2:end)=2*conn'; 


%Find out the dirichlet dofs (left edge and 3 nodes for roller)
%left
fixeddofs1 = 2*(find(coords(1,:)==0))-1; 

%right roller
fixeddofs2 = find(coords(2,:)==0);
len = length(fixeddofs2);
fixeddofs2 = 2*fixeddofs2(len-0:len); 
fixeddofs = [fixeddofs1 fixeddofs2];  
         
alldofs = 1:2*numnode;
freedofs = setdiff(alldofs,fixeddofs);


%to find out domain of influence of the node
xspac=Lb/nelx;
yspac=Hb/nely;


% %Traction boundary (point load)
% nodes_trac=find(coords(2,:)==40); 
% nodes_trac=nodes_trac(1:1);
% coords_trac = coords(:,nodes_trac); 
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
%filter for density
%coordinate of background cells centres
gs=gauss_domain(coords,numele,conn,1); 
coords_cells(1,:) = gs(4,:); 
coords_cells(2,:) = gs(5,:); 

%to find out domain of influence of the node
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

% figure
% patch('Faces',conn','Vertices',coords','LineWidth',1,...
%       'FaceColor',[0,153/255,153/255],'EdgeColor','k'); 
% hold on 
% plot(coords(1,:),coords(2,:),'ok','markersize',4,'markerfacecolor',[1,0,0])
% xlabel('x','fontsize',14); ylabel('y','fontsize',14);
% axis equal 
% 
