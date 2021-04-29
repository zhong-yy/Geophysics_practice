function [rs,phase]=MT2D_FEM_BILINEAR_TM(rho,DX,DY,f)
%rho 单元电阻率
%DX 行向量，横向的各个单元边长
%DY 行向量，纵向的各个单元边长
%f 频点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%钟乙源
%中南大学，2016,8,16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=4e-7*pi;
% ep=8.85e-12;
nx=size(rho,2);
ny=size(rho,1);
if (nx~=length(DX))||(ny~=length(DY))
    error(['输入电阻率单元矩阵rho与单元边长向量DX或DY不匹配，rho为',num2str(nx),'by',num2str(nx),'，DX=',num2str(length(DX)),'DY=',num2str(length(DY))]);
end
np=(nx+1)*(ny+1);%总节点数
ne=nx*ny;%单元总数

%节点编号
indexE=zeros(4,ne);
for ix=1:nx
    for iy=1:ny
        n=(ix-1)*ny+iy;
        ibegin=(ix-1)*(ny+1)+iy;   
        indexE(1,n)=ibegin;
        indexE(2,n)=ibegin+1;
        indexE(3,n)=indexE(2,n)+ny+1;
        indexE(4,n)=indexE(1,n)+ny+1;
    end    
end

for f1=1:size(f,2)
    K1=sparse(np,np);
    K2=sparse(np,np);
    K3=sparse(np,np);
    p=sparse(np,1);
    %%  形成K1  %%
    dx=ones(ny,1)*DX(:)';
    dy=DY(:)*ones(1,nx);
    dx=dx(:);
    dy=dy(:);
    a=dy./(6*dx).*rho(:);b=dx./(6*dy).*rho(:);
    ka=[2,;1;-1;-2
        1; 2;-2;-1
       -1;-2; 2; 1
       -2;-1; 1; 2];
    kb=[2;-2;-1; 1
       -2; 2; 1;-1
       -1; 1; 2;-2
        1;-1;-2; 2];
    ka=ka*ones(1,ne);
    kb=kb*ones(1,ne);
    a=ones(16,1)*a.';
    b=ones(16,1)*b.';
    ke=ka.*a+kb.*b;
    ir=(1:4)'*ones(1,4);
    ic=ir';
    ic=indexE(ic(:),:);
    ir=indexE(ir(:),:);
    K1=K1+sparse(ir(:),ic(:),ke(:) ,np,np);        

    %% 形成k2  %%%%%%%
    dx=ones(ny,1)*DX(:)';
    dy=DY(:)*ones(1,nx);
    dx=dx(:);
    dy=dy(:);
    a=1i*2*pi*f(f1)*mu*dx.*dy/36;
    ke=[4;2;1;2;
        2;4;2;1;
        1;2;4;2;
        2;1;2;4];
    a=ones(16,1)*a.';
    ke=ke*ones(1,ne).*a;
    ir=(1:4)'*ones(1,4);
    ic=ir';
    ic=indexE(ic(:),:);
    ir=indexE(ir(:),:);
    K2=K2+sparse(ir(:),ic(:),ke(:) ,np,np);    
        

    %%    形成k3   %%%%%%%
    h=ny:ny:ne;
    dx=DX(ceil(h/ny));
    mk=dx(:)'.*sqrt(-1i*2*pi*f(f1)*mu.*rho(h))/6; 
    mk=ones(4,1)*mk;
    ke=[2;1
        1;2];
    ke=ke*ones(1,nx).*mk;
    ir=[2;3;2;3];
    ic=[2;2;3;3];
    ir=indexE(ir(:),h);
    ic=indexE(ic(:),h);
    K3=K3+sparse(ir(:),ic(:),ke(:),np,np);    

    %% 
    K=K1-K2+K3;
    %%  u|AB=1 %%
    m=1:(ny+1):(np-ny);
    K(sub2ind([np,np],m,m))=K(sub2ind([np,np],m,m))*10^10;
    p=p+sparse(m,ones(1,nx+1),K(sub2ind([np,np],m,m))*1,np,1);    

    %% 解方程    
    tol=eps;
    maxit=1000;
    setup.type='ilutp';
    setup.milu='row';
    setup.droptol=1e-10;
    [L,U]=ilu(K,setup);
    x(:,f1)=bicgstab(K,p,tol,maxit,L,U);
%     x(:,f1)=K\p;
end
%% 求视电阻率和相位
surface=1:(ny+1):(np-ny);

d=sum(DY(1:3));
w=2*pi*f;
w=ones(length(surface),1)*w(:)';
u1=x(surface,:);
u2=x(surface+1,:);
u3=x(surface+2,:);
u4=x(surface+3,:);
uz=(-11*u1+18*u2-9*u3+2*u4)/(2*d);%这个求导公式的前提是4个点间的间隔相等
z=-rho(1)*uz./u1;%对于TE模式，应取Ey/(-Hx)
phase=angle(z);
phase=-phase;%因为时间因子取-iwt，所以计算结果的相位要取反才能得到真实的相位,
phase=phase+(phase<0)*2*pi;
phase=phase*180/pi;%范围为0到360
rs=(rho(1)*rho(1)*abs(uz./u1).^2)./(w*mu);
% for j=1:length(f)
%     for m=1:length(surface)
%         w=2*pi*f(j);
%         u1(m,j)=x(surface(m),j);
%         u2(m,j)=x(surface(m)+1,j);
%         u3(m,j)=x(surface(m)+2,j);
%         u4(m,j)=x(surface(m)+3,j);
%         uz(m,j)=(-11*u1(m,j)+18*u2(m,j)-9*u3(m,j)+2*u4(m,j))/(2*d);%这个求导公式的前提是4个点间的间隔相等
%         z(m,j)=-rho(1)*uz(m,j)/x(surface(m),j);
%         
%         phase(m,j)=angle(z(m,j));
%         phase(m,j)=-phase(m,j);
%         if(phase(m,j)<0)
%             phase(m,j)=phase(m,j)+2*pi;
%         end
%         phase(m,j)=phase(m,j)*180/pi;%范围为0到360
% %         if real(z(m,j))>0
% %             phase(m,j)=-atan(imag(z(m,j))/real(z(m,j)))*180/pi;%因为时间因子取-iwt，所以计算结果的相位要取反才能得到真实的相位
% %         elseif real(z(m,j))<0
% %             if imag(z(m,j))>0
% %                 phase(m,j)=-atan(imag(z(m,j))/real(z(m,j)))*180/pi+180;
% %             else
% %                 phase(m,j)=-atan(imag(z(m,j))/real(z(m,j)))*180/pi-180;
% %             end
% %         else 
% %             if imag(z(m,j))>0
% %                 phase(m,j)=90;
% %             else
% %                 phase(m,j)=-90;
% %             end
% %         end
%         rs(m,j)=(rho(1)*rho(1)*abs(uz(m,j)/x(surface(m),j))^2)/(w*mu);
%     end   
% end
% 
% semilogx(f,rs(5,:));
% figure
% imagesc([0,cumsum(DX)],log10(1./f),rs');colorbar
% figure;
% contourf([0,cumsum(DX)],log10(1./f),rs');colorbar