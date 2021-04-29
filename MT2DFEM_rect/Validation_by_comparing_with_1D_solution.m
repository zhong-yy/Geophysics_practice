clear
tic 

%frequencies
f=10.^(-3:0.25:3);

%resistivity
rho=zeros(30,40);

%number of air layers
air=10;

DYair=ones(1,air);
DYair(1:5)=2.^(7:-1:3)*100;
DYair(6:end)=[600,400,200,50,50];

%set up of a three layer model
rho(1:15,:)=1500;
rho(16:23,:)=100;
rho(24:30,:)=800;

DY=zeros(1,size(rho,1));
DX=zeros(1,size(rho,2));

DX(16:25)=500;
DX([11:15,26:30])=1000;
DX(1:5)=(5:-1:1)*2000;
DX(36:40)=(1:5)*2000;
DX([6:10,31:35])=1200;

DY(1:5)=50;
DY(6:10)=[200,500,1100,500,200];
DY(11:15)=50;
DY(16:23)=[50,150,250,350,350,250,150,50];
DY(24:30)=[50,150,250,500,1000,1500,2000];

%calculate the analytical solutions
[rs0,phase0]=Analytic_MT1D([1500,100,800],[3000,1600],f);

%start computation using FEM
[rs1,phase1]=MT2D_FEM_BILINEAR_TE(rho,air,DYair,DX,DY,f);
[rs2,phase2]=MT2D_FEM_BILINEAR_TM(rho,DX,DY,f);

[rs3,phase3]=MT2D_FEM_BIQUADRATIC_TE(rho,air,DYair,DX,DY,f);
[rs4,phase4]=MT2D_FEM_BIQUADRATIC_TM(rho,DX,DY,f);

%errors
error1=[norm(rs1(5,:)-rs0)/norm(rs0),norm(rs2(5,:)-rs0)/norm(rs0)];
error2=[norm(rs3(5,:)-rs0)/norm(rs0),norm(rs4(5,:)-rs0)/norm(rs0)];

%compare the analytical solutions and the numerical solutions for bilinear interpolation
figure%plot apparent resistivity
semilogx(1./f,rs0,'rs','MarkerSize',6.5,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0.1,0.2,0.1]);
hold on
semilogx(1./f,rs1(5,:),'bo','MarkerSize',6.5,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0.98,0.2,0.1]);
semilogx(1./f,rs2(5,:),'g+','MarkerSize',6.5,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0,0.4,0.85],'MarkerFaceColor',[0,0.4,0.85]);

legend('\fontsize{12}Analytic solution','\fontsize{12}FEM-TE mode','\fontsize{12}FEM-TM mode');
xlabel('Time(s)','FontSize',14);
ylabel('Apparent resistivity \rho_s(\Omega°§m)','FontSize',14);
title('Apparent Resistivity\rho_s(bilinear interpolation)','FontSize',16)
figure%œ‡Œª
semilogx(1./f,phase0,'rs','MarkerSize',6.5,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0.1,0.2,0.1]);
hold on
semilogx(1./f,phase1(5,:),'bo','MarkerSize',6.5,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0.98,0.2,0.1]);
semilogx(1./f,phase2(5,:),'g+','MarkerSize',6.5,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0,0.4,0.85],'MarkerFaceColor',[0,0.4,0.85]);
legend('\fontsize{12}Analytic solution','\fontsize{12}FEM-TE mode','\fontsize{12}FEM-TM mode');
xlabel('Time(s)','FontSize',14);
ylabel('Phase \circ','FontSize',14);
title('Phase(bilinear interpolation)','FontSize',16)

%compare the analytical solutions and the numerical solutions for biquadratic interpolation
figure
semilogx(1./f,rs0,'rs','MarkerSize',6.5,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0.1,0.2,0.1]);
hold on
semilogx(1./f,rs3(5,:),'bo','MarkerSize',6.5,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0.98,0.2,0.1]);
semilogx(1./f,rs4(5,:),'g+','MarkerSize',6.5,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0,0.4,0.85],'MarkerFaceColor',[0,0.4,0.85]);
legend('\fontsize{12}Analytic solution','\fontsize{12}FEM-TE mode','\fontsize{12}FEM-TM mode');
xlabel('Time(s)','FontSize',14);
ylabel('Apparent resistivity \rho_s(\Omega°§m)','FontSize',14);
title('Apparent Resistivity\rho_s(biquadratic interpolation)','FontSize',16);

figure
semilogx(1./f,phase0,'rs','MarkerSize',6.5,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0.1,0.2,0.1]);
hold on
semilogx(1./f,phase3(5,:),'bo','MarkerSize',6.5,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0.98,0.2,0.1]);
semilogx(1./f,phase4(5,:),'g+','MarkerSize',6.5,...
    'LineWidth',1.2,...
    'MarkerEdgeColor',[0,0.4,0.85],'MarkerFaceColor',[0,0.4,0.85]);
legend('\fontsize{12}Analytic solution','\fontsize{12}FEM-TE mode','\fontsize{12}FEM-TM mode');
xlabel('Time(s)','FontSize',14);
ylabel('Phase \circ','FontSize',14);
title('Phase(biquadratic interpolation)','FontSize',16);
toc


