close all

f=logspace(-4,4,41);
[rhos_analytic,phase_analytic]=Analytic_MT1D([300,100,1000],[1000,500],f);

DZ=10;
rho=ones(1,500)*1000;
rho(1,1:100)=300;
rho(1,101:150)=100;


tic
[rs,phase]=FDMt1d(rho,DZ,f);
toc

fsize=16;
figure('Position',[300 100 850 450]);
loglog(f,rhos_analytic,'ko','MarkerSize',4,'LineWidth',1.5);
hold on
loglog(f,rs,'r-','MarkerSize',4,'LineWidth',1.5);
set(gca,'XDir','reverse');
xlabel('Frequency (Hz)');
ylabel('Apparent resistivity (\Omega \cdot m)')
legend('Analytical solutions','FD solutions','location','best')
set(gca,'fontsize',fsize);
figure('Position',[300 100 850 450]);
plot(f,(rs-rhos_analytic')./rhos_analytic'*100,'Linewidth',2);
set(gca,'XDir','reverse');
xlabel('Frequency (Hz)');
ylabel('Relative errors (%)')
set(gca,'fontsize',fsize);
