clear
close all
fun=@(para,freq)(Analytic_MT1D(para(1:ceil(length(para)/2)),para(ceil(length(para)/2)+1:end),freq));
f=logspace(-4,4,40);
f=f';
nf=length(f);
rho={500;[500;200];[500;200;800]};
h={[];600;[600;1000]};
rhos=cell(1,4);
rho_inv=cell(1,4);
misfit=cell(1,4);
h_inv=cell(1,4);

rho_inv2=cell(3,4);
h_inv2=cell(3,4);
misfit2=cell(3,4);
maxit=500;


n_layer=4;%反演用的层数
fsize=13.1;
for j=1:3
    [rhos,phase_z]=fun([rho{j};h{j}],f);  
    %% 自动根据观测数据计算初值
    sample_space=nf/n_layer;%nlayer为层数，nf为频点数
    sample=floor(nf-sample_space/2:-sample_space:1);
    [h_ini,rho_ini]=Bostick(f(sample)',rhos(sample),phase_z(sample)*pi/180);
    %bostick反演要将相位的单位转换成弧度单位制
    h_ini=diff([0,h_ini(1:end-1)]);
    h_ini=h_ini';
    rho_ini=rho_ini';
    m_guess=[rho_ini;h_ini];

    %% 反演
    [minv,misfit,ite_num]=inversion_LM(fun,rhos,f,m_guess,50*ones(2*n_layer-1,1),ones(2*n_layer-1,1)*1e5,maxit,1e-5,3,0.25,0.75);
    showmsg(ite_num,maxit,misfit);
    
    rho_inv{j}=minv(1:n_layer);
    h_inv{j}=minv(n_layer+1:end);
    
    %% 绘图
    figure('Position',[80 100 1250 450]);
    subplot(121)
    loglog(f,rhos,'k*','MarkerSize',4,'LineWidth',1.5);
    set(gca,'XDir','reverse');
    hold on;
    loglog(f,fun(minv,f),'r-','MarkerSize',4,'LineWidth',2);
    set(gca,'XDir','reverse');%将频率坐标设置为由大到小
    xlim([min(f),max(f)])
    h_le=legend('Synthetic data (without noise)','Predicted data','location','best');
    title('Data fitting');
    set(h_le,'fontsize',13.6);
    % set(h_le,'box','off')
    xlabel('Frequency (Hz)');
    ylabel('Apparent resistivity (\Omega \cdot m)')
    
    set(gca,'fontsize',fsize);
    axis normal
    % 绘地层实际电阻率
    Depth_inv1=[1,cumsum(h_inv{j})'];
    Depth_inv2=[cumsum(h_inv{j})',100000];
    Depth_inv=[Depth_inv1;Depth_inv2];
    rho_plot_inv=[rho_inv{j}';rho_inv{j}'];
    
    Depth1=[1,cumsum(h{j}')];
    Depth2=[cumsum(h{j}'),100000];
    Depth=[Depth1;Depth2];
    rho_plot=[rho{j}';rho{j}'];
    if j==1
        ylim([200,1200]);
    elseif j==2
        ylim([150,700]);
    elseif j==3
        ylim([200,1120]);
    end
    
    subplot(122)
    loglog(Depth(:),rho_plot(:),'k','LineWidth',1.5);
    xlim([min(Depth(:)),max(Depth(:))]);
    hold on;
    loglog(Depth_inv(:),rho_plot_inv(:),'r--','LineWidth',1.5);
    if j==1
        ylim([200,1200]);
    elseif j==2
        ylim([100,1200]);
    elseif j==3
        ylim([150,1120]);
    end
    xlim([min(Depth_inv(:)),max(Depth_inv(:))]);
    h_le=legend('True model','Inversion result from LM method','location','best');
    title('Inversion result');
    set(h_le,'fontsize',14);
    % set(h_le,'box','off')
    xlabel('Depth (m)');
    ylabel('Resistivity (\Omega \cdot m)')
%    fsize=14;
    set(gca,'fontsize',fsize);
    axis normal
%     A=[rho_ini,rho_inv{j},rho{j}]
%     B=[h_ini,h_inv{j},h{j}]
end

err=[0.02,0.05,0.1];
colorline={[255,0,0]/255;[5,150,60]/255;[30,80,255]/255};
shape={'>';'s';'p'};
fsize=13.1;
for j=1:3
    [rhos0,phase_z0]=fun([rho{j};h{j}],f);
    figure('Position',[80 100 1250 450]);
    ax2=subplot(122);
    Depth1=[1,cumsum(h{j}')];
    Depth2=[cumsum(h{j}'),100000];
    Depth=[Depth1;Depth2];
    rho_plot=[rho{j}';rho{j}'];
    loglog(Depth(:),rho_plot(:),'k','LineWidth',1.5);
    xlim([min(Depth(:)),max(Depth(:))]);
    hold on;    
    for p=1:3
        rhos=rhos0+err(p)*rhos.*randn(size(rhos,1),size(rhos,2));
        %阻抗相位也加上噪音，因为下面会用阻抗相位做bostick反演以计算LM法的初值
        phase_z=phase_z0+err(p)*phase_z0.*randn(size(phase_z0,1),size(phase_z0,2));
        
    %     n_layer=length(rho{j});    
        %% 自动根据观测数据计算初值
        sample_space=nf/n_layer;
    %     sample=floor(sample_space/2:sample_space:nf);
    % %     sample=[1,sample_space,2*sample_space];
    %     sample=fliplr(sample);
        sample=floor(nf-sample_space/2:-sample_space:1);
        [h_ini,rho_ini]=Bostick(f(sample)',rhos(sample),phase_z(sample)*pi/180);
        %bostick反演要将相位的单位转换成弧度单位制
        h_ini=diff([0,h_ini(1:end-1)]);
        h_ini=h_ini';
        rho_ini=rho_ini';
        m_guess=[rho_ini;h_ini];

        %% Inversion
        [minv,misfit2{p,j},ite_num]=inversion_LM(fun,rhos,f,m_guess,50*ones(2*n_layer-1,1),ones(2*n_layer-1,1)*1e5,maxit,1e-4,3,0.25,0.75);
        showmsg(ite_num,maxit,misfit2{p,j});

        rho_inv2{p,j}=minv(1:n_layer);
        h_inv2{p,j}=minv(n_layer+1:end);

        %% Plotting
       
        ax1=subplot(121);
        loglog(f,rhos,shape{p},'MarkerSize',4,...
                     'MarkerEdgeColor',colorline{p},...
                     'MarkerFaceColor',colorline{p});
        
        set(gca,'XDir','reverse');
        hold on;
        loglog(f,fun(minv,f),'-','MarkerSize',4,'LineWidth',2,'Color',colorline{p});
        set(gca,'XDir','reverse');%将频率坐标设置为由大到小
        xlim([min(f),max(f)])
        
        % set(h_le,'box','off')
        xlabel('Frequency (Hz)');
        ylabel('Apparent resistivity (\Omega \cdot m)')
        set(gca,'fontsize',fsize);

        if j==1
            ylim([100,1500]);
        elseif j==2
            ylim([120,700]);
        elseif j==3
            ylim([150,1120]);
        end
        
        axis normal
        % 绘实际的和反演得到的地层电阻率结构
        Depth_inv1=[1,cumsum(h_inv2{p,j})'];
        Depth_inv2=[cumsum(h_inv2{p,j})',100000];
        Depth_inv=[Depth_inv1;Depth_inv2];
        rho_plot_inv=[rho_inv2{p,j}';rho_inv2{p,j}'];

        ax2=subplot(122);
        loglog(Depth_inv(:),rho_plot_inv(:),'--','LineWidth',2,'Color',colorline{p});
        if j==1
            ylim([100,3500]);
        elseif j==2
            ylim([45,2000]);
        elseif j==3
            ylim([50,2000]);
        end
        xlim([min(Depth_inv(:)),max(Depth_inv(:))]);
        title(['Inversion models from data with ',num2str(100*err(p)),'% Gassian noise']);
        set(h_le,'fontsize',14);
        % set(h_le,'box','off')
        xlabel('Depth (m)');
        ylabel('Resistivity (\Omega \cdot m)')
        
        set(gca,'fontsize',fsize);
        axis normal
        
    end
    h_e1=legend(ax1,'Data with 2% noise','Predicted data','Data with 5% noise','Predicted data','Data with 10% noise','Predicted data','location','best');
    set(h_e1,'fontsize',12);
    title(ax1,'Data fitting');
    h_e2=legend(ax2,'True model','Inversion model from data with 2% noise','Inversion model from data with 5% noise','Inversion model from data with 10% noise','location','best');
    set(h_e2,'fontsize',12);
    title(ax2,'Inversion models');
    
end