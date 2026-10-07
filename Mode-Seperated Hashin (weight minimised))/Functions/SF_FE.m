function [phi, dphix, dphiy]=SF_FE(gg, x, conn)

for i=1:4
    je=conn(i,gg(1)); xe(i)=x(1,je); ye(i)=x(2,je);
end

one=ones(1,4);
psiJ=[-1 +1 +1 -1]; etaJ=[-1 -1  +1  +1];
psi=gg(2); eta=gg(3); %position in psi, eta coordinates i.e. (-1,+1 and -1,+1)
phi=0.25*(one+psi*psiJ).*(one+eta*etaJ);

NJpsi=.25*psiJ.*(one+eta*etaJ);
NJeta=.25*etaJ.*(one+psi*psiJ);
xpsi=NJpsi*xe';ypsi=NJpsi*ye';xeta=NJeta*xe';yeta=NJeta*ye';
jcob_mat=[xpsi  ypsi; xeta  yeta];
Bmat1=jcob_mat\[NJpsi;NJeta];
dphix=Bmat1(1,:);
dphiy=Bmat1(2,:);
