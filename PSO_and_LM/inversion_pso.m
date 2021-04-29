function [m,err,convergence,ite]=inversion_pso(f,y_obs,x,N_size,m_lower,m_upper,vmax,w,tol,maxit)
% =========================================================================
% ����Ⱥ�㷨
% Particle swarm algorithm
%
% ����ԴYiyuan Zhong
% 2017, 11
%
% ���ϴ�ѧ�������ѧ����Ϣ����ѧԺ
% Central South University
% =========================================================================
% f: ���ݺ������
% x: ���ݵ���֪���룬����MT��Ƶ�ʣ�f�ĵ�����ʽ����Ϊf(m,x)
%
% N_size: ������
% m_lower, m_upper: ������������
% w=[w0,ag,al] PSO����������
%
% tol: ����ϲ�С�����ֵʱ��ֹͣ����
% maxit: ����������
% =========================================================================
% m: ���ݵĽ�
% convergence: �������У���һ��Ϊ���������ڶ���Ϊÿ�ε��������
% ite: ���յ�������
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
            %����Ⱥ�弫ֵ
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
