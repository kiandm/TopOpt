function [coords, conn, edofMat, numnode, numele, freedofs, F, W] = problem_setup_cantilever(rmin_phys)
% Cantilever beam (Fig. 7): length Lb, height Lb/2, fully clamped at
% x=0 (both DOFs), point load at mid-height of the free (right) edge,
% pointing in -y (0.25*Lb above it, 0.25*Lb below it).

ndiv = 100;                    % elements along the length
Lb = 100;                      % physical length
Hb = Lb/2;                     % physical height (cantilever aspect ratio from Fig. 7)
size_cell = Lb/ndiv;
nelx = ndiv;
nely = round(Hb/size_cell);
rmin = rmin_phys / size_cell;  % cell-count radius, used only for the candidate search below

[x2d,y2d]=meshgrid(0:size_cell:Lb, Hb:-size_cell:0);

size2d=size(x2d);
numnode=size2d(1)*size2d(2);
numele=nelx*nely;

count=1;
for i=1:size2d(2)
    for j=1:size2d(1)
        coords(1,count)= x2d(j,i);
        coords(2,count)= y2d(j,i);
        node_num(j,i)=count;
        count=count+1;
    end
end

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

% Fully clamped left edge (cantilever wall) - BOTH dofs, unlike the MBB symmetry BC
fixed_nodes = find(coords(1,:)==0);
fixeddofs=zeros(1,2*length(fixed_nodes));
fixeddofs(1,1:2:end)=2*fixed_nodes-1;
fixeddofs(1,2:2:end)=2*fixed_nodes;

alldofs = 1:2*numnode;
freedofs = setdiff(alldofs,fixeddofs);

xspac=size_cell;
yspac=size_cell;

% % Point load at mid-height of the free (right) edge, per Fig. 7
% mid_y = Hb/2;
% node_load = find(coords(1,:)==Lb & coords(2,:)==mid_y);
% P_total = -200;   % physical load magnitude - tune so stresses approach the strength allowables
% F = sparse(2*node_load, 1, P_total, 2*numnode, 1);
% Load distributed over 5 nodes, centred on the mid-height of the free (right) edge
mid_y = Hb/2;
n_node_span = 5;
half_span = floor(n_node_span/2)*size_cell;      % = 2*size_cell either side of mid_y
y_lo = mid_y - half_span;
y_hi = mid_y + half_span;

nodes_trac = find(coords(1,:)==Lb & coords(2,:)>=y_lo & coords(2,:)<=y_hi);
[~, order] = sort(coords(2,nodes_trac));
nodes_trac = nodes_trac(order);                  % must be sorted along the edge for the segment loop below
coords_trac = coords(:,nodes_trac);

P_total = -1000;                                   % total applied load - tune so stresses approach the strength allowables
L_trac = coords_trac(2,end) - coords_trac(2,1);   % physical length spanned by the traction nodes
pres = P_total / L_trac;

F = zeros(2*numnode,1);
for el = 1:length(nodes_trac)-1
    area_el = abs(coords(2,nodes_trac(el+1)) - coords(2,nodes_trac(el)));
    F(2*nodes_trac(el))   = F(2*nodes_trac(el))   + pres*area_el/2;
    F(2*nodes_trac(el+1)) = F(2*nodes_trac(el+1)) + pres*area_el/2;
end
F = sparse(F);

%filter for density
gs=gauss_domain(coords,numele,conn,1);
coords_cells(1,:) = gs(4,:);
coords_cells(2,:) = gs(5,:);

dm_cells(1,1:numele)=rmin*(xspac*ones(1,numele));
dm_cells(2,1:numele)=rmin*(yspac*ones(1,numele));

W=zeros(numele,numele);
for cc=1:numele
    gpos=coords_cells(:,cc);
    [v,~]=nodes_in_support(numele, coords_cells, gpos, dm_cells);
    xi=coords_cells(:,v);
    difx=abs((gpos(1,1)-xi(1,:)));
    dify=abs((gpos(2,1)-xi(2,:)));
    dif=sqrt(difx.^2 + dify.^2);
    wij = max(rmin_phys - dif, 0) / rmin_phys;
    W(cc,v)=wij;
end
W=W./sum(W,2);
W=sparse(W);
end