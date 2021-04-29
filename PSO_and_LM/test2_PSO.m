clear
close all
fun=@(para,freq)(Analytic_MT1D(para(1:ceil(length(para)/2)),para(ceil(length(para)/2)+1:end),freq));
f=logspace(-4,4,40);
f=f';
nf=length(f);

rho={500;[500;200];[500;200;800]};
h={[];600;[600;1000]};

rho_inv=cell(1,3);
misfit=cell(1,3);
h_inv=cell(1,3);

rho_inv2=cell(3,3);
h_inv2=cell(3,3);
misfit2=cell(3,3);
convergence=cell(4,3);
maxit=200;


n_layer=5;
N_size=50;
fsize=12.5;
for j=1:3
    [rhos,phase_z]=fun([rho{j};h{j}],f);
    %rhos=rhos+0.05*rhos.*randn(size(rhos,1),size(rhos,2));
%     n_layer=length(rho{j});    

    %% 反演
    %[minv,misfit,ite_num]=inversion_LM(fun,rhos,f,m_guess,10*ones(2*n_layer-1,1),ones(2*n_layer-1,1)*1e5,maxit,1e-5,3,0.25,0.75);
    m_lower=100*ones(2*n_layer-1,1);
    m_upper=ones(2*n_layer-1,1)*1.5e3;
    vmax=(m_upper-m_lower)/20;
    w=[0.687,2.25,1.125];
%     w=[0.8,1.9,2.1];
    [minv,misfit{j},convergence{1,j},ite_num]=inversion_pso(fun,rhos,f,N_size,m_lower,m_upper,vmax,w,1e-5,maxit);
    showmsg(ite_num,maxit,misfit{j});
    
    rho_inv{j}=minv(1:n_layer);
    h_inv{j}=minv(n_layer+1:end);
    
    %% 绘图
    figure('Position',[80 100 1250 450]);
    subplot(121)
    loglog(f,rhos,'k*-','MarkerSize',4,'LineWidth',1.5);
    set(gca,'XDir','reverse');
    hold on;
    loglog(f,fun(minv,f),'r-','MarkerSize',4,'LineWidth',1.5);
    set(gca,'XDir','reverse');%将频率坐标设置为由大到小
    xlim([min(f),max(f)])
    h_le=legend('Synthetic data (without noise)','Predicted data from the PSO result','location','best');
    title('Data fitting');
    set(h_le,'fontsize',fsize);
    % set(h_le,'box','off')
    xlabel('Frequency (Hz)');
    ylabel('Apparent resistivity (\Omega \cdot m)')
    set(gca,'fontsize',fsize);
    if j==1
        ylim([200,1200]);
    elseif j==2
        ylim([150,700]);
    elseif j==3
        ylim([200,1120]);
    end
    
    
    % 绘地层实际电阻率
    Depth_inv1=[1,cumsum(h_inv{j})'];
    Depth_inv2=[cumsum(h_inv{j})',100000];
    Depth_inv=[Depth_inv1;Depth_inv2];
    rho_plot_inv=[rho_inv{j}';rho_inv{j}'];
    
    Depth1=[1,cumsum(h{j}')];
    Depth2=[cumsum(h{j}'),100000];
    Depth=[Depth1;Depth2];
    rho_plot=[rho{j}';rho{j}'];

    
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
    h_le=legend('True model','Inversion model from PSO method','location','best');
    set(h_le,'fontsize',fsize);
    title('Inversion model');
    % set(h_le,'box','off')
    xlabel('Depth (m)');
    ylabel('Resistivity (\Omega \cdot m)')
    set(gca,'fontsize',fsize);
end


err=[0.02,0.05,0.1];
colorline={[255,0,0]/255;[5,150,60]/255;[30,80,255]/255};
shape={'>';'s';'p'};
for j=1:3
    rhos0=fun([rho{j};h{j}],f);
    figure('Position',[80 100 1250 450]);
    Depth1=[1,cumsum(h{j}')];
    Depth2=[cumsum(h{j}'),100000];
    Depth=[Depth1;Depth2];
    rho_plot=[rho{j}';rho{j}'];    
    ax2=subplot(122);
    loglog(Depth(:),rho_plot(:),'k','LineWidth',1.5);
    xlim([min(Depth(:)),max(Depth(:))]);
    hold on;   
    for p=1:3
        %加噪
        rhos=rhos0+err(p)*rhos.*randn(size(rhos,1),size(rhos,2));
        %% 反演
        %[minv,misfit,ite_num]=inversion_LM(fun,rhos,f,m_guess,10*ones(2*n_layer-1,1),ones(2*n_layer-1,1)*1e5,maxit,1e-5,3,0.25,0.75);
        m_lower=100*ones(2*n_layer-1,1);
        m_upper=ones(2*n_layer-1,1)*1.5e3;
        vmax=(m_upper-m_lower)/20;
        w=[0.687,2.25,1.125];

        

        [minv,misfit2{p,j},convergence{p+1,j},ite_num]=inversion_pso(fun,rhos,f,N_size,m_lower,m_upper,vmax,w,1e-3,maxit);
        showmsg(ite_num,maxit,misfit2{p,j});

        rho_inv2{p,j}=minv(1:n_layer);
        h_inv2{p,j}=minv(n_layer+1:end);

        %% 绘图
        
        ax1=subplot(121);
        loglog(f,rhos,shape{p},'MarkerSize',4,...
                     'MarkerEdgeColor',colorline{p},...
                     'MarkerFaceColor',colorline{p});
        set(gca,'XDir','reverse');
        hold on;
        loglog(f,fun(minv,f),'-','MarkerSize',4,'LineWidth',2,'Color',colorline{p});
        hold on;
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
        

        % 绘地层实际电阻率
        Depth_inv1=[1,cumsum(h_inv2{p,j})'];
        Depth_inv2=[cumsum(h_inv2{p,j})',100000];
        Depth_inv=[Depth_inv1;Depth_inv2];
        rho_plot_inv=[rho_inv2{p,j}';rho_inv2{p,j}'];
        subplot(122)
        loglog(Depth_inv(:),rho_plot_inv(:),'--','LineWidth',2,'Color',colorline{p});
        if j==1
            ylim([100,3500]);
        elseif j==2
            ylim([45,2000]);
        elseif j==3
            ylim([50,2000]);
        end
        xlim([min(Depth_inv(:)),max(Depth_inv(:))]);
%         h_le=legend('True model','Inversio models from PSO methods','location','best');
%         set(h_le,'fontsize',14);
%         title(['Inversion models wi',num2str(100*err(p)),'%高斯噪声数据的反演地层模型']);
%         % set(h_le,'box','off')
%         xlabel('深度 (m)');
%         ylabel('地层电阻率 (\Omega \cdot m)')
%         set(gca,'fontsize',fsize);

%         figure
%         semilogy(convergence(1,:),convergence(2,:),'b.','linewidth',2);
%         xlim([min(convergence(1,:)),max(convergence(1,:))]);
%         xlabel('迭代次数');
%         ylabel('相对拟合差');
%         set(gca,'fontsize',14)
    end
    h_e1=legend(ax1,'Data with 2% noise','Predicted data','Data with 5% noise','Predicted data','Data with 10% noise','Predicted data','location','best');
    set(h_e1,'fontsize',12);
    title(ax1,'Data fitting');
    h_e2=legend(ax2,'True model','Inversion model from data with 2% noise','Inversion model from data with 5% noise','Inversion model from data with 10% noise','location','best');
    set(h_e2,'fontsize',12);
    title(ax2,'Inversion models');
end
colorline2={[120,90,30]/255;[255,0,0]/255;[5,150,60]/255;[30,80,255]/255};
figure('Position',[80 100 1350 450]);
for j=1:3
    subplot(1,3,j);
    for p=1:4
        semilogy(convergence{p,j}(1,:),convergence{p,j}(2,:),'.',...
            'Color',colorline2{p},...
            'linewidth',2);
        hold on;
        xlim([min(convergence{p,j}(1,:)),max(convergence{p,j}(1,:))]);   
    end
    xlabel('Iteration');
    ylabel('Relative error');
    set(gca,'fontsize',14)         
    legend('Withoud noise','2% noise','5% noise','10% noise','location','best');
end
subplot(131)
title('Homogeneous half space');
subplot(132)
title('Two layer model');
subplot(133)
title('Three layer model');