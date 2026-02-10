function [N]=shapefunc(gp,ngp,nen)
if ngp==9
  g2=sqrt(3/5);
  xsi=[-1 -1 -1  0  0  0  1  1  1].'*g2;
  eta=[-1  1  0 -1  1  0 -1  1  0].'*g2;
elseif ngp==4
  g2=1/sqrt(3);
  xsi=[-1 -1 1 1].'*g2;
  eta=[-1 1 1 -1].'*g2;  
end
N=zeros(nen,1);
if nen==8
  N(1)=1/4*(1-xsi(gp))*(1-eta(gp))*(-xsi(gp)-eta(gp)-1);
  N(2)=1/2*(1-xsi(gp))*(1-eta(gp)^2);
  N(3)=1/4*(1-xsi(gp))*(1+eta(gp))*(-xsi(gp)+eta(gp)-1);
  N(4)=1/2*(1-xsi(gp)^2)*(1+eta(gp));
  N(5)=1/4*(1+xsi(gp))*(1+eta(gp))*(xsi(gp)+eta(gp)-1);
  N(6)=1/2*(1+xsi(gp))*(1-eta(gp)^2);
  N(7)=1/4*(1+xsi(gp))*(1-eta(gp))*(xsi(gp)-eta(gp)-1);
  N(8)=1/2*(1-xsi(gp)^2)*(1-eta(gp));
elseif nen==9
  N(1)= 1/4*xsi(gp)*(xsi(gp)-1)*eta(gp)*(eta(gp)-1);
  N(2)=-1/2*xsi(gp)*(xsi(gp)-1)*(eta(gp)+1)*(eta(gp)-1);
  N(3)= 1/4*xsi(gp)*(xsi(gp)-1)*eta(gp)*(eta(gp)+1);
  N(4)=-1/2*(xsi(gp)+1)*(xsi(gp)-1)*eta(gp)*(eta(gp)+1);
  N(5)= 1/4*xsi(gp)*(xsi(gp)+1)*eta(gp)*(eta(gp)+1);
  N(6)=-1/2*xsi(gp)*(xsi(gp)+1)*(eta(gp)+1)*(eta(gp)-1);
  N(7)= 1/4*xsi(gp)*(xsi(gp)+1)*eta(gp)*(eta(gp)-1);
  N(8)=-1/2*(xsi(gp)+1)*(xsi(gp)-1)*eta(gp)*(eta(gp)-1);
  N(9)=(xsi(gp)+1)*(xsi(gp)-1)*(eta(gp)+1)*(eta(gp)-1);
elseif nen==4
  N(1)=(1-xsi(gp))*(1-eta(gp))/4;
  N(2)=(1-xsi(gp))*(1+eta(gp))/4;
  N(3)=(1+xsi(gp))*(1+eta(gp))/4;
  N(4)=(1+xsi(gp))*(1-eta(gp))/4;
end
