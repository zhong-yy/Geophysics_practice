function [m,err,convergence,ite]=inversion_pso(f,y_obs,x,N_size,m_lower,m_upper,vmax,w,tol,maxit)
% =========================================================================
% 粒子群算法
% Particle swarm algorithm
%
% 钟乙源Yiyuan Zhong
% 2017, 11
%
% 中南大学，地球科学与信息物理学院
% Central South University
% =========================================================================
% f: 正演函数句柄
% x: 正演的已知输入，比如MT的频率，f的调用形式必须为f(m,x)
%
% N_size: 例子数
% m_lower, m_upper: 参数的上下限
% w=[w0,ag,al] PSO的三个参数
%
% tol: 当拟合差小于这个值时，停止迭代
% maxit: 最大迭代次数
% =========================================================================
% m: 反演的解
% convergence: 包含两行，第一行为迭代数，第二行为每次迭代的误差
% ite: 最终迭代次数
% =========================================================================



if norm(y_obs,2)>eps
    %relative error
   misfit=@(mvalue)(norm(f(mvalue,x)-y_obs,2)/norm(y_obs,2));
%     misfit=@(mvalue)(sqrt(sum((f(mvalue,x)-y_obs).^2./(y_obs.^2)))/length(y_obs));
    %misfit=@(mvalue)(norm((f(mvalue,x)-y_obs)/length(y_obs),2));
else
    %absolute error
    misfit=@(mvalue)(norm((f(mvalue,x)-y_obs)/length(y_obs),2));
end

%% initialization
n=length(m_lower);%number of parameters

%first convert m_lower and m_upper into column vectors
m_lower=m_lower(:);
m_upper=m_upper(:);

%initialize particles
particle=m_lower*ones(1,N_size)+...
    (m_upper-m_lower)*ones(1,N_size).*rand(n,N_size);
%v=particle/100;
v=zeros(n,N_size);

%compute the fitness values of initial particles
p_bestfit=zeros(1,N_size);
for k=1:N_size
    p_bestfit(k)=misfit(particle(:,k));
end
%
[g_bestfit,best_index]=min(p_bestfit);

gbest=particle(:,best_index);%global best position
pbest=particle;%previous best position for each particle

M_UPPER=m_upper*ones(1,N_size);
M_LOWER=m_lower*ones(1,N_size);

vmax=vmax(:);
vmax=vmax*ones(1,N_size);

for ite=1:maxit
    r1=ones(n,1)*rand(1,N_size);
    r2=ones(n,1)*rand(1,N_size);

 %   w(1)=0.999^ite*rand()+0.2;
    v=w(1)*v+...
        w(2)*r1.*(gbest*ones(1,N_size)-particle)+...
        w(3)*r2.*(pbest-particle);
    
    v(v>vmax)=vmax(v>vmax);
    v(v<-vmax)=-vmax(v<-vmax);
    %v(v>vmax)=rand()
    
    
    particle=particle+v;
%     for par_i=1:N_size
%         [m1,rms]=inversion_1d(f,y_obs,x,particle(:,par_i),10,1e-8,3,0.25,0.75);
%         m1=m1(:);
%         particle(:,par_i)=m1;
%     end;
        
    
    
    idx=particle>M_UPPER;
    particle(idx)=M_UPPER(idx);
    
    idx=particle<M_LOWER;
    particle(idx)=M_LOWER(idx);
    
    
    for k=1:N_size
        fitness=misfit(particle(:,k));
        if fitness<p_bestfit(k)
            p_bestfit(k)=fitness;
            pbest(:,k)=particle(:,k);
            %更新群体极值
            if fitness<g_bestfit
                g_bestfit=fitness;
                gbest=particle(:,k);
            end
        end
    end
    
    convergence(1,ite)=ite;
    convergence(2,ite)=g_bestfit;
       
    err=g_bestfit;
    if err<tol
        break;
    end   
end
m=gbest;
