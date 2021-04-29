function [rs,phase]=MT2D_FEM_BIQUADRATIC_TE(rho,air,DYair,DX,DY,f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rho: a matrix, representing the values of resistivity in the mesh
%
%air: an integer number, representing the number of air layers
%
%DX: a row vector，representing lengths of element edges in the lateral direction
%
%DY: a row vector, representing lengths of elements edges in the subsurface in the vertical direction
%
%DYair: lengths of element edges in the air layers in the vertical direction
%f: a list of frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%钟乙源 Yiyuan Zhong
%中南大学 Central South University，2016,8,16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=4e-7*pi;
ep=8.85e-12;
nx=size(rho,2);
ny=size(rho,1);
if (nx~=length(DX))||(ny~=length(DY))
    error(['The matrix rho is not consistent with edge lengths DX and DY. rho is ',num2str(ny),'by',num2str(nx),'，DX=',num2str(length(DX)),',DY=',num2str(length(DY))]);
end
rho=[ones(air,size(rho,2))*1e10;rho];
DY=[DYair,DY];


ny=size(rho,1);

np=(2*ny+1)*(nx+1)+(ny+1)*nx;%Total number of nodes
ne=nx*ny;%number of elements

%indices of node
indexE=zeros(8,ne);
iy=(1:ny)'*ones(1,nx);
ibegin=(ones(ny,1)*(1:nx)-1)*(3*ny+2)+2*(1:ny)'*ones(1,nx)-1;
indexE(1,:)=ibegin(:);
indexE(2,:)=indexE(1,:)+2;
indexE(3,:)=indexE(2,:)+3*ny+2;
indexE(4,:)=indexE(1,:)+3*ny+2;
indexE(5,:)=indexE(1,:)+1;
indexE(6,:)=indexE(2,:)+2*ny+1-iy(:)';
indexE(7,:)=indexE(4,:)+1;
indexE(8,:)=indexE(6,:)-1;

for f1=1:size(f,2)
    K1=sparse(np,np);
    K2=sparse(np,np);
    K3=sparse(np,np);
    p=sparse(np,1);
    %%  construct K1  %%
    dx=ones(ny,1)*DX(:)';
    dy=DY(:)*ones(1,nx);
    dx=dx(:);
    dy=dy(:);
    a=dy./(90*dx)/(1i*2*pi*f(f1)*mu);
    b=dx./(90*dy)/(1i*2*pi*f(f1)*mu);
    ka=[52;  17;  23;  28;   6; -40;  -6; -80;
        17;  52;  28;  23;   6; -80;  -6; -40;
        23;  28;  52;  17;  -6; -80;   6; -40;
        28;  23;  17;  52;  -6; -40;   6; -80;
        6;   6;  -6;  -6;  48;   0; -48;   0;
        -40; -80; -80; -40;   0; 160;   0;  80;
        -6;  -6;   6;   6; -48;   0;  48;   0;
        -80; -40; -40; -80;   0;  80;   0; 160];
    kb=[  52;  28;  23;  17; -80;  -6; -40;   6;
        28;  52;  17;  23; -80;   6; -40;  -6;
        23;  17;  52;  28; -40;   6; -80;  -6;
        17;  23;  28;  52; -40;  -6; -80;   6;
        -80; -80; -40; -40; 160;   0;  80;   0;
        -6;   6;   6;  -6;   0;  48;   0; -48;
        -40; -40; -80; -80;  80;   0; 160;   0;
        6;  -6;  -6;   6;   0; -48;   0;  48];
    ka=ka*ones(1,ne);
    kb=kb*ones(1,ne);
    a=ones(64,1)*a.';
    b=ones(64,1)*b.';
    ke=ka.*a+kb.*b;
    ir=(1:8)'*ones(1,8);
    ic=ir';
    ic=indexE(ic(:),:);
    ir=indexE(ir(:),:);
    K1=K1+sparse(ir(:),ic(:),ke(:) ,np,np);
    clear dx dy a b ir r1 ic c1;
    %% construct k2  %%%%%%%
    ke=[6,2,3,2,-6,-8,-8,-6;
        2,6,2,3,-6,-6,-8,-8;
        3,2,6,2,-8,-6,-6,-8;
        2,3,2,6,-8,-8,-6,-6;
        -6,-6,-8,-8,32,20,16,20;
        -8,-6,-6,-8,20,32,20,16;
        -8,-8,-6,-6,16,20,32,20;
        -6,-8,-8,-6,20,16,20,32];
    ke=ke(:)*ones(1,ne);
    dx=ones(ny,1)*DX(:)';
    dy=DY(:)*ones(1,nx);
    dx=dx(:);
    dy=dy(:);
    a=(1./rho(:)-1i*2*pi*f(f1)*ep).*dx.*dy/180;
    a=ones(64,1)*a';
    ke=ke.*a;
    ir=(1:8)'*ones(1,8);
    ic=ir';
    ic=indexE(ic(:),:);
    ir=indexE(ir(:),:);
    K2=K2+sparse(ir(:),ic(:),ke(:) ,np,np);
    clear dx dy ir r1 ic c1;
    %%    construct k3   %%%%%%%
    h=ny:ny:ne;
    dx=DX(ceil(h/ny));
    mk=dx(:)'.*sqrt(-1i*2*pi*f(f1)*mu./rho(h))/(1i*2*pi*f(f1)*mu)/30;
    mk=ones(9,1)*mk;
    ke=[4;-1;2;
        -1;4;2;
        2;2;16]*ones(1,nx).*mk;
    ir=[2;3;6]*ones(1,3);
    ic=ir';
    ir=indexE(ir(:),h);
    ic=indexE(ic(:),h);
    K3=K3+sparse(ir(:),ic(:),ke(:),np,np);
    clear h dx a b c mk;
    %%
    K=K1-K2+K3;
    %%  u|AB=1 %%
    m=[1:(3*ny+2):(np-2*ny),2*ny+2:(3*ny+2):(np-3*ny-1)];
    K(sub2ind([np,np],m,m))=K(sub2ind([np,np],m,m))*10^10;
    p=p+sparse(m,ones(1,2*nx+1),K(sub2ind([np,np],m,m))*1,np,1);       
    clear j m;
    %% solve the equation
%     figure
%     spy(K)
    tol=eps;
    maxit=1000;
    setup.type='ilutp';
    setup.milu='row';
    setup.droptol=1e-10;
    [L,U]=ilu(K,struct('type','ilutp','droptol',1e-9));
    x(:,f1)=bicgstab(K,p,tol,maxit,L,U);
    %           x(:,f1)=v\p;
    %setup.udiag=1;
    
    %x(:,f1)=bicgstab(v,p,tol,maxit,L,U);
    %[L,U]=ilu(v)
end
%% calculate the apparent resistivity and phase
surface=(2*air+1):(3*ny+2):np;

d=DY(air+1)+1/2*DY(air+2);
w=2*pi*f;
w=ones(length(surface),1)*w(:)';
u1=x(surface,:);
u2=x(surface+1,:);
u3=x(surface+2,:);
u4=x(surface+3,:);
uz=(-11*u1+18*u2-9*u3+2*u4)/(2*d);%which requires that the 4 points near the Earth-air interface are equally space
z=1i*mu*w.*u1./uz;%For TE mode,Ey/(-Hx)
phase=angle(z);
phase=-phase;%Because the time harmonic dependence is exp(-iwt), the resulting phases should be negated
phase=phase+(phase<0)*2*pi;
phase=phase*180/pi;% make the phase range from 0 to 360

rs=mu*w.*(abs(u1./uz).^2);