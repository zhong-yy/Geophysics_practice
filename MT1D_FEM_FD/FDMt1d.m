function [rs,phase]=FDMt1d(rho,dL,f)

mu=4e-7*pi;
ep=8.85e-12;
ne=length(rho);%element number
np=ne+1;%node number
u1=1;
for f1=1:length(f)
    K=zeros(np-1,np-1);
    p=zeros(np-1,1);
    sig=(1/rho(1)+1/rho(2))/2;
    w=2*pi*f(f1);
    ksquare=-i*w*mu*(sig-i*w*ep);
    K(1,1:2)=[-ksquare-2/(dL^2),1/(dL^2)];
    p(1)=-u1/(dL^2);
    for h=2:ne-1
        sig=(1/rho(h+1)+1/rho(h))/2;
        K(h,h-1:h+1)=[1/(dL^2),-2/(dL^2)+i*w*mu*(sig-i*w*ep),1/(dL^2)];
    end
    K(end,[np-2,np-1])=[-1/dL,sqrt(-i*w*mu*(1/rho(end)-i*w*ep))+1/dL];
    u(:,f1)=ThomasAlgorithm(K,p);
end
% figure
% spy(K)
u=[ones(1,length(f));u];


rs=zeros(length(f),1);
phase=zeros(length(f),1);
z=zeros(length(f),1);
uz=zeros(length(f),1);
d=3*dL;
for j=1:length(f)
    w=2*pi*f(j);
    uz(j)=(-11*u(1,j)+18*u(2,j)-9*u(3,j)+2*u(4,j))/(2*d);
    z(j)=1i*w*mu*u(1,j)/uz(j);
    rs(j)=(abs(z(j))^2)/(w*mu);
    phase(j)=angle(z(j));
    phase(j)=-phase(j);%因为时间因子取-iwt，所以计算结果的相位要取反才能得到真实的相位
    if(phase(j)<0)
        phase(j)=phase(j)+2*pi;
    end
     phase(j)=phase(j)*180/pi;%范围为0到360
end