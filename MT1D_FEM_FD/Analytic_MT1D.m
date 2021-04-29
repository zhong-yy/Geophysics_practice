function [rhos,phase]=Analytic_MT1D(rho,h,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1ά��״���ʽ�����
%Author: ����Դ
%���ϴ�ѧ�������ѧ����Ϣ����ѧԺ
%E-mail: zhongyiyuan@csu.edu.cn
%Date: 2017 10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1=length(rho);
n2=length(h);
n3=length(f);
omega=2*pi*f;
if n1~=n2+1
    error('��������ʸ����������ƥ��')
end
mu=4*pi*10^(-7);
k=zeros(n1,n3);


for u=1:n1
    k(u,:)=sqrt(-i*omega*mu/rho(u));
end



R=ones(1,n3);
for u=n2:-1:1
    A=1+exp(-2*k(u,:)*h(u)); 
    B=1-exp(-2*k(u,:)*h(u));
    R=(A*sqrt(rho(u+1)/rho(u)).*R+B)./(B*sqrt(rho(u+1)/rho(u)).*R+A);
end
rhos=rho(1)*abs(R).^2;
phase=-(angle(R)-pi/4);%ȡ���õ�ʵ����λ
for k=1:length(phase)
    if(phase(k)<0)
       phase(k)=phase(k)+2*pi;
    end
end
phase=phase*180/pi;%��ΧΪ0��360��

