close all
clear

DX=[8000,6000,4000,2000,ones(1,320)*500,2000,4000,6000,8000];
DY=[ones(1,8)*250,ones(1,8)*500,ones(1,20)*200,500,500,1000,1000,1000,2000,2000,5000,5000,10000,12000,10000,5000,5000,3000,1000,1000,1000,1000,3000,5000,5000,10000,20000,40000,40000];
% DY=[ones(1,500)*300,ones(1,50)*1000];
air=10;
DYair=[30000,20000,10000,8000,6000,4000,2000,1000,800,500];
nx=length(DX);
ny=length(DY);
rho=ones(ny,nx)*10;
rho(1:53,:)=100;
rho(22:31,121:156)=0.1;
rho(22:31,185:208)=0.1;
f=0.1;

tic
% [rs_te,phase_te]=MT2D_FEM_BIQUADRATIC_TE(rho,air,DYair,DX,DY,f);
[rs_te,phase_te]=MT2D_FEM_BILINEAR_TE(rho,air,DYair,DX,DY,f);
toc

tic
% [rs_tm,phase_tm]=MT2D_FEM_BIQUADRATIC_TM(rho,DX,DY,f);
[rs_tm,phase_tm]=MT2D_FEM_BILINEAR_TM(rho,DX,DY,f);
toc


rs1=rs_tm(5:2:325,:);
rs2=rs_te(5:2:325,:);
phase1=phase_tm(5:2:325,:);
phase2=phase_te(5:2:325,:);
x=(-80000:1000:80000)/1000;

H = load('2dmt_station_results_from_RenZY');
y = H(:,2)./1000;
app_xy = H(:,4);
app_yx = H(:,5);
pha_xy = -H(:,6);
pha_yx = -H(:,7);
pha_xy=pha_xy+(pha_xy<0)*2*pi;
pha_yx=pha_yx+(pha_yx<0)*2*pi;

msize=3;
lsize=1.5;
fsize=14;

figure('Position',[0,0,800,500])
plot(x,app_yx,'Color',[0.98,0.2,0.1],'LineWidth',lsize);hold on
plot(x,app_xy,'Color',[0,0.4,0.85],'LineWidth',lsize);hold on;
plot(x,rs1,'ro','MarkerSize',4,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0.98,0.2,0.1]);
plot(x,rs2,'ko','MarkerSize',4,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0,0.4,0.85]);

h=legend('Ren \rho_{TM}','Ren \rho_{TE}','My program \rho_{TM}','My program \rho_{TE}');
set(h,'Box','off','Location','SouthEast','FontSize',12);
ylabel('apparent resistivity\rho_s(\Omega¡¤m)','FontSize',14);
xlabel('stations(km)','FontSize',14);
ylim([20 120]);
xlim([min(x),max(x)])
set(gca,'FontSize',fsize);
title('Apparent Resistivity','FontSize',16);

figure('Position',[0,0,800,500]);
plot(x,pha_yx,'-','Color',[0.98,0.2,0.1],'LineWidth',lsize);hold on;
plot(x,pha_xy,'-','Color',[0,0.4,0.85],'LineWidth',lsize);

hold on;
plot(x,phase1,'o','MarkerSize',4,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0.98,0.2,0.1]);
plot(x,phase2,'o','MarkerSize',4,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0,0.4,0.85]);

h=legend('Ren \phi_{TM}','Ren \phi_{TE}','My program \phi_{TM}','My program \phi_{TE}');
set(h,'Box','off','Location','NorthEast','FontSize',12);
ylabel('Phase(degrees)','FontSize',14);
xlabel('stations(km)','FontSize',14);
ylim([40 90]);
title('Phase','FontSize',16);
set(gca,'FontSize',fsize);

figure;
h=plot(x,100*(rs1-app_yx)./app_yx,'-',x,100*(rs2-app_xy)./app_xy,'-');
set(h,'MarkerSize',6,'LineWidth',1.5);
ylabel('relative error(%)','FontSize',14);
xlabel('stations(km)','FontSize',14);
title('Error of Resistivity','FontSize',14);
legend('TM','TE')
set(gca,'FontSize',fsize);
xlim([min(x),max(x)])



figure;
h=plot(x,100*(phase1-pha_yx)./pha_yx,'-',x,100*(phase2-pha_xy)./(pha_xy),'-');
set(h,'MarkerSize',6,'LineWidth',1.5);
ylabel('relative error(%)','FontSize',14);
xlabel('stations(km)','FontSize',14);
title('Error of Phase','FontSize',14);
legend('TM','TE')
set(gca,'FontSize',fsize);
xlim([min(x),max(x)])


%calculate mt solutions for a series of frequencies
tic
f2=logspace(-3,3,13);
[rs_te2,phase_te2]=MT2D_FEM_BILINEAR_TE(rho,air,DYair,DX,DY,f2);
% [rs_te2,phase_te2]=MT2D_FEM_BILINEAR_TE(rho,air,DYair,DX,DY,f2);
toc

tic
[rs_tm2,phase_tm2]=MT2D_FEM_BILINEAR_TM(rho,DX,DY,f2);
% [rs_tm2,phase_tm2]=MT2D_FEM_BILINEAR_TM(rho,DX,DY,f2);
rs_tm2=rs_tm2(5:2:325,:);
rs_te2=rs_te2(5:2:325,:);
phase_tm2=phase_tm2(5:2:325,:);
phase_te2=phase_te2(5:2:325,:);
toc

rs_tm2=rs_tm2.';
rs_te2=rs_te2.';
phase_tm2=phase_tm2';
phase_te2=phase_te2';


figure
pcolor(x,log10(f2),rs_tm2);hold on;
shading interp
contour(x,log10(f2),rs_tm2,'w');
colormap(jet)
hcb=colorbar;
set(get(hcb,'Title'),'string','\Omega\cdot m');
set(hcb,'FontSize',14);
caxis([1,160])
xlabel('x(km)');
ylabel('log_{10}f(Hz)');

figure
pcolor(x,log10(f2),rs_te2);hold on;
shading interp
contour(x,log10(f2),rs_te2,'w');
colormap(jet)
hcb=colorbar;
set(get(hcb,'Title'),'string','\Omega\cdot m');
set(hcb,'FontSize',14);
caxis([1,160])