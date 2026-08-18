function gs=gauss_domain(coords,numcell,conn,nQ)%Domains 

if nQ==1
    gauss(1,1) = 0;
    gauss(2,1) = 2;
elseif nQ==2
    gauss(1,1) =-0.5773502691896257645091488;
    gauss(1,2) =+0.5773502691896257645091488;
    
    gauss(2,1) = 1;
    gauss(2,2) = 1;
end

count=1; 
for e=1:numcell
    for i=1:4
        je=conn(i,e); xe(i)=coords(1,je); ye(i)=coords(2,je);
    end
    one=ones(1,4);
    psiJ=[-1 +1 +1 -1]; etaJ=[-1 -1  +1  +1];
    for i=1:nQ %4 gauss points
        for j=1:nQ %4 gauss points
            eta=gauss(1,i); psi=gauss(1,j);
            N=0.25*(one+psi*psiJ).*(one+eta*etaJ);
            xq=N*xe'; yq=N*ye';% Position of the gauss point
            
            NJpsi=.25*psiJ.*(one+eta*etaJ);
            NJeta=.25*etaJ.*(one+psi*psiJ);
            xpsi=NJpsi*xe';ypsi=NJpsi*ye';xeta=NJeta*xe';yeta=NJeta*ye';
            jcob=xpsi*yeta-xeta*ypsi;
            
            gs(1,count)=e; 
            gs(2,count)=psi;
            gs(3,count)=eta;
            gs(4,count)=xq;
            gs(5,count)=yq;
            gs(6,count)=gauss(2,i)*gauss(2,j);
            gs(7,count)=jcob;
            count=count+1;
        end
    end
end


% figure
% patch('Faces',conn','Vertices',coords','LineWidth',1,...
%       'FaceColor',[0,153/255,153/255],'EdgeColor','k'); 
% hold on 
% plot(coords(1,:),coords(2,:),'ok','markersize',8,'markerfacecolor',[1,0,0])
% xlabel('x','fontsize',14); ylabel('y','fontsize',14);
% axis equal 
% plot(gs(4,:),gs(5,:),'xg','markersize',8,'linewidth',2)
