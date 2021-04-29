function [rs,phase]=MT2D_FEM_BILINEAR_TE(rho,air,DYair,DX,DY,f)
%rho 单元电阻率
%air 整数，空气层数
%DX 行向量，横向的各个单元边长
%DY 行向量，纵向的各个单元边长
%DYair 空气层各单元的纵向长度
%f 频点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%钟乙源
%中南大学，2016,8,16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=4e-7*pi;
ep=8.85e-12;
nx=size(rho,2);
ny=size(rho,1);
if (nx~=length(DX))||(ny~=length(DY))
    error(['输入电阻率单元矩阵rho与单元边长向量DX或DY不匹配，rho为',num2str(ny),'by',num2str(nx),'，DX=',num2str(length(DX)),',DY=',num2str(length(DY))]);
end
rho=[ones(air,size(rho,2))*1e10;rho];
DY=[DYair,DY];

nx=size(rho,2);
ny=size(rho,1);

np=(nx+1)*(ny+1);%总节点数
ne=nx*ny;%单元总数

%节点编号
indexE=ones(4,ne);
ibegin=(ones(ny,1)*(1:nx)-1)*(ny+1)+(1:ny)'*ones(1,nx);
indexE(1,:)=ibegin(:);
indexE(2,:)=ibegin(:)+1;
indexE(3,:)=indexE(2,:)+ny+1;
indexE(4,:)=indexE(1,:)+ny+1;
% for ix=1:nx
%     for iy=1:ny
%         n=(ix-1)*ny+iy;
%         ibegin=(ix-1)*(ny+1)+iy;
%         indexE(1,n)=ibegin;
%         indexE(2,n)=ibegin+1;
%         indexE(3,n)=indexE(2,n)+ny+1;
%         indexE(4,n)=indexE(1,n)+ny+1;
%     end    
% end

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
    a=dy./(12*1i*pi*f(f1)*mu*dx);b=dx./(12*1i*pi*f(f1)*mu*dy);
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
%     for h=1:ne
%         dx=DX(ceil(h/ny));
%         
%         if rem(h,ny)==0
%             dy=DY(end);
%         else
%             dy=DY(rem(h,ny));
%         end
%         a=dy/(6*dx);b=dx/(6*dy);        
%         ke=[2*(a+b)  a-2*b    -a-b   -2*a+b
%              a-2*b  2*(a+b) -2*a+b     -a-b
%               -a-b  -2*a+b  2*(a+b)   a-2*b
%             -2*a+b    -a-b    a-2*b  2*(a+b)]/(1i*2*pi*f(f1)*mu);
%         ir=indexE(:,h)*ones(1,4);
%         ic=ir';
%         K1=K1+sparse(ir(:),ic(:),ke(:),np,np);
%     end    
    %% 形成k2  %%%%%%%
    dx=ones(ny,1)*DX(:)';
    dy=DY(:)*ones(1,nx);
    dx=dx(:);
    dy=dy(:);
    a=(1./rho(:)-1i*2*pi*f(f1)*ep).*dx.*dy/36;
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
    
%     for h=1:ne
%         dx=DX(ceil(h/ny));
%         
%         if rem(h,ny)==0
%             dy=DY(end);
%         else
%             dy=DY(rem(h,ny));
%         end        
%         ke=(1/rho(h)-i*2*pi*f(f1)*ep)*dx*dy*[4,2,1,2;2,4,2,1;1,2,4,2;2,1,2,4]/36;
%         for ir=1:4
%            r1=indexE(ir,h);
%            for ic=1:4
%                c1=indexE(ic,h);
%                K2(r1,c1)=K2(r1,c1)+ke(ir,ic);
%            end
%         end        
%     end        
    %%    形成k3   %%%%%%%
    h=ny:ny:ne;
    dx=DX(ceil(h/ny));
    mk=dx(:)'.*sqrt(-1i*2*pi*f(f1)*mu./rho(h))/(12*1i*pi*f(f1)*mu);
    mk=ones(4,1)*mk;
    ke=[2;1
        1;2];
    ke=ke*ones(1,nx).*mk;
    ir=[2;3;2;3];
    ic=[2;2;3;3];
    ir=indexE(ir(:),h);
    ic=indexE(ic(:),h);
    K3=K3+sparse(ir(:),ic(:),ke(:),np,np);
%     for h=ny:ny:ne
%         dx=DX(ceil(h/ny));
%         mk=dx*sqrt(-i*2*pi*f(f1)*mu/rho(h))/6/(i*2*pi*f(f1)*mu);        
%         ke=[2,1
%             1,2]*mk;
%         ir=indexE([2,3],h)*ones(1,2);
%         ic=ir';
%         K3=K3+sparse(ir(:),ic(:),ke(:),np,np);
%     end
    %% 
    K=K1-K2+K3;
    %%  u|AB=1 %%
    m=1:(ny+1):(np-ny);
    K(sub2ind([np,np],m,m))=K(sub2ind([np,np],m,m))*10^10;
    p=p+sparse(m,ones(1,nx+1),K(sub2ind([np,np],m,m))*1,np,1);
%     for j=1:nx+1
%         m=(j-1)*(ny+1)+1;
%         K(sub2ind([np,np],m,m))=K(sub2ind([np,np],m,m))*10^10;
%         p=p+sparse(m,1,K(m,m)*1,np,1);
%     end
    %% solve the equations
    tol=eps;
    maxit=1000;
    setup.type='ilutp';
    setup.milu='row';
    setup.droptol=1e-10;
    [L,U]=ilu(K,setup);
    x(:,f1)=bicgstab(K,p,tol,maxit,L,U);
%     x(:,f1)=K\p;
    %setup.udiag=1;
    
    %x(:,f1)=bicgstab(v,p,tol,maxit,L,U); 
end
%% calculate the apparent resistivity and phase
surface=air+1:(ny+1):np;

d=sum(DY(air+1:air+3));
w=2*pi*f;
w=ones(length(surface),1)*w(:)';
u1=x(surface,:);
u2=x(surface+1,:);
u3=x(surface+2,:);
u4=x(surface+3,:);
uz=(-11*u1+18*u2-9*u3+2*u4)/(2*d);%numerical differentiation. This formula requires that the 4 points near the Earth-air interface are equally space
z=1i*mu*w.*u1./uz;%对于TE模式，应取Ey/(-Hx)
phase=angle(z);
phase=-phase;%Because the time harmonic dependence is exp(-iwt), the resulting phases should be negated
phase=phase+(phase<0)*2*pi;
phase=phase*180/pi;%范围为0到360

rs=mu*w.*(abs(u1./uz).^2); 
% T=-ux./uz;
