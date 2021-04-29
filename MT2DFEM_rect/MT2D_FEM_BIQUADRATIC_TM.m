function [rs,phase]=MT2D_FEM_BIQUADRATIC_TM(rho,DX,DY,f)
%rho: a matrix, representing the values of resistivity in the mesh
%DX: a row vector，representing lengths of element edges in the lateral direction
%DY: a row vector, representing lengths of elements edges in the subsurface in the vertical direction
%f: a list of frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%钟乙源Yiyuan Zhong
%中南大学Central South University，2016,8,16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu=4e-7*pi;
% ep=8.85e-12;
nx=size(rho,2);
ny=size(rho,1);
if (nx~=length(DX))||(ny~=length(DY))
    error(['输入电阻率单元矩阵rho与单元边长向量DX或DY不匹配，rho为',num2str(nx),'by',num2str(nx),'，DX=',num2str(length(DX)),'DY=',num2str(length(DY))]);
end
np=(2*ny+1)*(nx+1)+(ny+1)*nx;%Total number of nodes
ne=nx*ny;%number of elements

%indexing
indexE=ones(8,ne);
% ix=(1:nx)'*ones(1,ny);
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
% for ix=1:nx
%     for iy=1:ny
%         n=(ix-1)*ny+iy;
%         ibegin=(ix-1)*(3*ny+2)+2*iy-1;
%         indexE(1,n)=ibegin;
%         indexE(2,n)=ibegin+2;
%         indexE(3,n)=indexE(2,n)+3*ny+2;
%         indexE(4,n)=indexE(1,n)+3*ny+2;
%         indexE(5,n)=ibegin+1;
%         indexE(6,n)=indexE(2,n)+2*ny+1-iy;
%         indexE(7,n)=indexE(4,n)+1;
%         indexE(8,n)=indexE(6,n)-1;
%     end    
% end

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
    a=dy./(90*dx).*rho(:);b=dx./(90*dy).*rho(:); 
    ka=[52;  17;  23;  28;   6; -40;  -6; -80;
        17;  52;  28;  23;   6; -80;  -6; -40;
        23;  28;  52;  17;  -6; -80;   6; -40;
        28;  23;  17;  52;  -6; -40;   6; -80;
         6;   6;  -6;  -6;  48;   0; -48;   0;
       -40; -80; -80; -40;   0; 160;   0;  80;
        -6;  -6;   6;   6; -48;   0;  48;   0;
       -80; -40; -40; -80;   0;  80;   0; 160];
    kb=[ 52;  28;  23;  17; -80;  -6; -40;   6;
         28;  52;  17;  23; -80;   6; -40;  -6;
         23;  17;  52;  28; -40;   6; -80;  -6;
         17;  23;  28;  52; -40;  -6; -80;   6;
        -80; -80; -40; -40; 160;   0;  80;   0;
         -6;   6;   6;  -6;   0;  48;   0; -48;
        -40; -40; -80; -80;  80;   0; 160;   0;
          6;  -6;  -6;   6;   0; -48;   0;  48];
    ka=ka*ones(1,ne);
    kb=kb*ones(1,ne);
    a=ones(64,1)*a';
    b=ones(64,1)*b';
    ke=ka.*a+kb.*b;
    ir=(1:8)'*ones(1,8);
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
%         a=dy/(90*dx)*rho(h);b=dx/(90*dy)*rho(h);        
%         ke=[52*a+52*b,17*a+28*b,23*a+23*b,28*a+17*b,6*a-80*b,-40*a-6*b,-6*a-40*b,-80*a+6*b;
%             17*a+28*b,52*a+52*b,28*a+17*b,23*a+23*b,6*a-80*b,-80*a+6*b,-6*a-40*b,-40*a-6*b;
%             23*a+23*b,28*a+17*b,52*a+52*b,17*a+28*b,-6*a-40*b,-80*a+6*b,6*a-80*b,-40*a-6*b;
%             28*a+17*b,23*a+23*b,17*a+28*b,52*a+52*b,-6*a-40*b,-40*a-6*b,6*a-80*b,-80*a+6*b
%             6*a-80*b,6*a-80*b,-6*a-40*b,-6*a-40*b,48*a+160*b,0,-48*a+80*b,0;
%             -40*a-6*b,-80*a+6*b,-80*a+6*b,-40*a-6*b,0,160*a+48*b,0,80*a-48*b;
%             -6*a-40*b,-6*a-40*b,6*a-80*b,6*a-80*b,-48*a+80*b,0,48*a+160*b,0
%             -80*a+6*b,-40*a-6*b,-40*a-6*b,-80*a+6*b,0,80*a-48*b,0,160*a+48*b];
%         ir=indexE(:,h)*ones(1,8);
%         ic=ir';
%         K1=K1+sparse(ir(:),ic(:),ke(:),np,np);
%         
% 
%     end    
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
    a=(1i*2*pi*f(f1)*mu)*dx.*dy/180;
    a=ones(64,1)*(a.');
    ke=ke.*a;
    ir=(1:8)'*ones(1,8);
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
%         ke=[6,2,3,2,-6,-8,-8,-6;
%             2,6,2,3,-6,-6,-8,-8;
%             3,2,6,2,-8,-6,-6,-8;
%             2,3,2,6,-8,-8,-6,-6;
%             -6,-6,-8,-8,32,20,16,20;
%             -8,-6,-6,-8,20,32,20,16;
%             -8,-8,-6,-6,16,20,32,20;
%             -6,-8,-8,-6,20,16,20,32];
%         ke=(1i*2*pi*f(f1)*mu)*dx*dy/180*ke;
%         ir=indexE(:,h)*ones(1,8);
%         ic=ir';
%         K2=K2+sparse(ir(:),ic(:),ke(:) ,np,np);
%     end        
    %%    construct k3   %%%%%%%
    h=ny:ny:ne;
    dx=DX(ceil(h/ny));
    mk=dx(:)'.*sqrt(-1i*2*pi*f(f1)*mu*rho(h))/30;
    mk=ones(9,1)*mk;
    ke=[4;-1;2;
        -1;4;2;
        2;2;16]*ones(1,nx).*mk;
    ir=[2;3;6]*ones(1,3);
    ic=ir';
    ir=indexE(ir(:),h);
    ic=indexE(ic(:),h);
    K3=K3+sparse(ir(:),ic(:),ke(:),np,np);
%     for h=ny:ny:ne
%         dx=DX(ceil(h/ny));
%         mk=dx*sqrt(-1i*2*pi*f(f1)*mu*rho(h))/30;
%         ke=[4,-1, 2
%            -1, 4, 2
%             2, 2,16]*mk;
%         ir=indexE([2,3,6],h)*ones(1,3);
%         ic=ir';
%         K3=K3+sparse(ir(:),ic(:),ke(:),np,np);
%     end
    %% 
    K=K1-K2+K3;
%     v=sparse(K);
    %%  u|AB=1 %%
    m=[1:(3*ny+2):(np-2*ny),2*ny+2:(3*ny+2):(np-3*ny-1)];
    K(sub2ind([np,np],m,m))=K(sub2ind([np,np],m,m))*10^10;
    p=p+sparse(m,ones(1,2*nx+1),K(sub2ind([np,np],m,m))*1,np,1);    
%     for j=1:nx+1
%         m=(j-1)*(3*ny+2)+1;
%         K(sub2ind([np,np],m,m))=K(sub2ind([np,np],m,m))*10^10;
% %         K(m,m)=K(m,m)*10^10;
%         p=p+sparse(m,1,K(m,m)*1,np,1);
%     end
%     for j=1:nx
%         m=j*(3*ny+2)-ny;
%         K(sub2ind([np,np],m,m))=K(m,m)*10^10;
% %         K(m,m)=K(m,m)*10^10;
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

end
%% calculate the resistivity and phase
surface=1:(3*ny+2):(np-ny);


d=DY(1)+1/2*DY(2);
w=2*pi*f;
w=ones(length(surface),1)*w(:)';
u1=x(surface,:);
u2=x(surface+1,:);
u3=x(surface+2,:);
u4=x(surface+3,:);
uz=(-11*u1+18*u2-9*u3+2*u4)/(2*d);%which requires that the 4 points near the Earth-air interface are equally space
z=-rho(1)*uz./u1;
phase=angle(z);
phase=-phase;%Because the time harmonic dependence is exp(-iwt), the resulting phases should be negated
phase=phase+(phase<0)*2*pi;
phase=phase*180/pi;% make the phase range from 0 to 360
rs=(rho(1)*rho(1)*abs(uz./u1).^2)./(w*mu);