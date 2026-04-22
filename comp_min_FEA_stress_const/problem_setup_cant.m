function [coords, conn, edofMat, numnode, numele, freedofs, F, W]= ...
                                    problem_setup_cant(rmin)

Hb=50; Lb=2*Hb; 
nely=50;  nelx=2*nely; 

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

numele=length(conn(1,:)); 

edofMat=zeros(numele,8); 
edofMat(:,1:2:end)=2*conn'-1'; 
edofMat(:,2:2:end)=2*conn'; 

%Find out the dirichlet dofs
%left
fixed_nodes=find(coords(1,:)==0); 
fixeddofs=zeros(1, 2*length(fixed_nodes));
fixeddofs(1:2:end)=2*fixed_nodes-1; 
fixeddofs(2:2:end)=2*fixed_nodes; 

alldofs = 1:2*numnode;
freedofs = setdiff(alldofs,fixeddofs);


%to find out domain of influence of the node
xspac=Lb/nelx;
yspac=Hb/nely;

%to find out domain of influence of the node
dm_cells(1,1:numele)=rmin*(xspac*ones(1,numele));
dm_cells(2,1:numele)=rmin*(yspac*ones(1,numele));


%Traction boundary (point load)
nodes_trac=find(coords(1,:)==Lb); 
mid_node_index= ceil(length(nodes_trac)/2); 
nodes_trac=nodes_trac(mid_node_index-2:1:mid_node_index+2);
coords_trac = coords(:,nodes_trac); 


pres=-100; 
F=zeros(2*numnode,1); 

%Traction boundary (point load)
area=abs(coords(2,nodes_trac(2))-coords(2,nodes_trac(1))); 
for el=1:length(nodes_trac)-1
    F(2*nodes_trac(el))  =F(2*nodes_trac(el))  +pres*area/2;
    F(2*nodes_trac(el+1))=F(2*nodes_trac(el+1))+pres*area/2;

end

%filter for density
%coordinate of background cells centres
gs=gauss_domain(coords,numele,conn,1); 
coords_cells(1,:) = gs(1,:); 
coords_cells(2,:) = gs(2,:); 

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


figure
patch('Faces',conn','Vertices',coords','LineWidth',1,...
      'FaceColor',[0,153/255,153/255],'EdgeColor','k'); 
hold on 
plot(coords(1,:),coords(2,:),'ok','markersize',4,'markerfacecolor',[1,0,0])
plot(coords(1,fixed_nodes),coords(2,fixed_nodes),'ok','markersize',4,...
    'markerfacecolor',[0,0,0])
plot(coords_trac(1,:),coords_trac(2,:),'ok','markersize',8,'markerfacecolor',[0,0.8,0])

xlabel('x','fontsize',14); ylabel('y','fontsize',14);
axis equal 
