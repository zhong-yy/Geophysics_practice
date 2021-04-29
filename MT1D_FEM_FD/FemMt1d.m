function [rs,phase]=FemMt1d(rho,DZ,f)
mu=4e-7*pi;
ep=8.85e-12;
ne=length(rho);%element number
np=ne+1;%node number
if ne~=length(DZ)
    error(['rho: ',num2str(ne),'£¬DZ: ',num2str(length(DZ))]);
end
for f1=1:length(f)
    K1=zeros(np);
    K2=zeros(np);
    K3=zeros(np);
    %%  K1  %%    
    for h=1:ne
        k1e=[1,-1;
            -1,1]/DZ(h);
        K1([h,h+1],[h,h+1])=K1([h,h+1],[h,h+1])+k1e;
    end
    %%  K2  %%
    for h=1:ne
        ksquare=-1i*2*pi*f(f1)*mu*(1/rho(h)-1i*2*pi*f(f1)*ep);
        k2e=ksquare*DZ(h)*...
            [1/3,1/6;
            1/6,1/3];
        K2([h,h+1],[h,h+1])=K2([h,h+1],[h,h+1])+k2e;
    end
    K3(np,np)=sqrt(-1i*2*pi*f(f1)*mu*(1/rho(ne)-1i*2*pi*f(f1)*ep));
    %% 
    K=K1+K2+K3;
    
    %%  u|AB=1 %%
%     K(1,1)=K(1,1)*10^12;
%     v(1,1)=v(1,1)*10^12;
%     p(1)=K(1,1)*1;
    p=zeros(np-1,1);%Earth-air interface: u1=1£¬move u1 to the right hand side of equations
    p(1)=-1*K(2,1);
    K(1,:)=[];
    K(:,1)=[];
    v=sparse(K);
   %% solve equations
   u(:,f1)=ThomasAlgorithm(K,p);
end
u=[ones(1,length(f));u];


rs=zeros(length(f),1);
phase=zeros(length(f),1);
z=zeros(length(f),1);
uz=zeros(length(f),1);
d=sum(DZ(1:3));
for j=1:length(f)
    w=2*pi*f(j);
    uz(j)=(-11*u(1,j)+18*u(2,j)-9*u(3,j)+2*u(4,j))/(2*d);
    z(j)=1i*w*mu*u(1,j)/uz(j);
    rs(j)=(abs(z(j))^2)/(w*mu);
    phase(j)=angle(z(j));
    phase(j)=-phase(j);%Because of time dependence: -iwt, the phase should be negated
    if(phase(j)<0)
        phase(j)=phase(j)+2*pi;
    end
     phase(j)=phase(j)*180/pi;%make the phase range from 0 to 360
end