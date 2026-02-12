function [coords, conn, edofMat, numnode, numele, freedofs, F, H]= ...
    problem_setup_Lbrac(rmin)                     
                           
ndiv=50;  %150, 100, 50, 30 use multiple of 5                           
size_cell=100/ndiv;

[x2d1 y2d1]=meshgrid(0:size_cell:50,   100:-size_cell:0); 
[x2d2 y2d2]=meshgrid(50+size_cell:size_cell:100, 50:-size_cell:0); 

size1=size(x2d1); numnode1=size1(1)*size1(2); 
size2=size(x2d2); numnode2=size2(1)*size2(2); 
numnode=numnode1+numnode2; 

x1=size(2,numnode1); x2=size(2,numnode2); 
for i=1:numnode1
    x1(1,i)=x2d1(i); 
    x1(2,i)=y2d1(i); 
end 

for i=1:numnode2
    x2(1,i)=x2d2(i); 
    x2(2,i)=y2d2(i); 
end 
x=[x1 x2]; coords=x; 

%to find the nodes for each cell
count1=1;
count2=1; 
for i=1:50/size_cell %length wise 
    for j=1:100/size_cell %width wise
        conn1(1,count1)=count2+1;
        conn1(2,count1)=conn1(1,count1)+ndiv+1;
        conn1(3,count1)=conn1(2,count1)-1;
        conn1(4,count1)=conn1(1,count1)-1;
        count1=count1+1;
        count2=count2+1;
    end
    count2=count2+1;
end
                            
count1=1;
count2=ndiv*((50/size_cell)+1)+1; 
for i=1:(100-50)/size_cell
    for j=1:50/size_cell
        conn2(1,count1)=count2+1;
        conn2(2,count1)=conn2(1,count1)+(50/size_cell)+1;
        conn2(3,count1)=conn2(2,count1)-1;
        conn2(4,count1)=conn2(1,count1)-1;
        count1=count1+1;
        count2=count2+1;
    end
    count2=count2+1;
end
conn=[conn1 conn2]; 
numele=length(conn(1,:)); 

edofMat=zeros(numele,8); 
edofMat(:,1:2:end)=2*conn'-1'; 
edofMat(:,2:2:end)=2*conn'; 

%Find out the dirichlet dofs (left edge and 3 nodes for roller)
%left
fixed_nodes = find(coords(2,:)==100); 
fixeddofs=zeros(1,2*length(fixed_nodes)); 
fixeddofs(1, 1:2:end)=2*fixed_nodes-1; 
fixeddofs(1, 2:2:end)=2*fixed_nodes; 

alldofs = 1:2*numnode;
freedofs = setdiff(alldofs,fixeddofs);

% %to find out domain of influence of the node
xspac=size_cell;
yspac=size_cell;

% dm(1,1:numnode)=dmax*(xspac*ones(1,numnode));
% dm(2,1:numnode)=dmax*(yspac*ones(1,numnode));

%Traction boundary (point load)
nodes_trac=find(coords(2,:)==50 & coords(1,:)>=95); 
% nodes_trac=find(coords(1,:)==100 & coords(2,:)>=30); 
coords_trac=coords(:,nodes_trac); 

pres=-1; 
F=zeros(2*numnode,1); 

%Traction boundary (point load)
area=abs(coords(1,nodes_trac(2))-coords(1,nodes_trac(1))); 
for el=1:length(nodes_trac)-1
    F(2*nodes_trac(el))  =F(2*nodes_trac(el))  +pres*area/2;
    F(2*nodes_trac(el+1))=F(2*nodes_trac(el+1))+pres*area/2;

end


%filter for density
%coordinate of background cells centres
gs=gauss_domain(coords,numele,conn,1); 
coords_cells(1,:) = gs(4,:); 
coords_cells(2,:) = gs(5,:); 

%to find out domain of influence of the node
dm_cells(1,1:numele)=rmin*(xspac*ones(1,numele));
dm_cells(2,1:numele)=rmin*(yspac*ones(1,numele));

H=zeros(numele,numele); 
for cc=1:numele
    gpos=coords_cells(:,cc); 
    [v,L]=nodes_in_support(numele, coords_cells, gpos, dm_cells);
    xi=coords_cells(:,v); 
    difx=abs((gpos(1,1)-xi(1,:))); 
    dify=abs((gpos(2,1)-xi(2,:))); 
    dif=sqrt(difx.^2 + dify.^2); 
    rij=dif./sqrt(dm_cells(1,v).^2 + dm_cells(2,v).^2);
    wij=(rmin-rij)./rmin; 
    H(cc,v)=wij; 
end
H=H./sum(H,2); 
H=sparse(H); 

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
                        