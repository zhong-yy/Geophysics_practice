function [m,misfit,ite]=inversion_LM(f,y_obs,x_input,m0,m_lower,m_upper,maxit,tol,lambda,R_low,R_up)
% =========================================================================
% LM法，或阻尼非线性最小二乘法
% 钟乙源
% 2017, 11
% 中南大学，地球科学与信息物理学院
% ==============================参数=======================================
% f: 为正问题函数句柄，f的调用形式必须为f(m0,x_input)
% y_obs: 为观测值（已知的输出）
% x_input为f的已知输入，
% m为需要反演的参数
% 
% m0: 为给定的初值
% m_lower: 模型参数下界
% m_upper: 模型参数上界
% maxit: 为最大迭代数;
% tol: 为两次迭代得到的参数的相对误差，如果两次迭代所得结果很接近则终止迭代
% lambda: 初始阻尼因子
% R_low: 用于调整阻尼因子的比例因子的下界，一般为0.25
% R_up: 用于调整阻尼因子的比例因子的上界，一般为0.75
% ==============================返回值=====================================
% m: 达到最大迭代次数或达到给定精度后，得到的模型参数向量
% misfit: 相对偏差
% ite: 最后迭代次数
% =========================================================================
N=length(m0);
M=length(y_obs);

if M~=size(x_input,1)
    error('输入参数不匹配');
end

dy=y_obs-f(m0,x_input);
%必须确保y_obs和f(m0,x_input)的输出都是行向量或都是列向量，否则会出错

m=m0(:);
m_upper=m_upper(:);
m_lower=m_lower(:);
for ite=1:maxit
    J=zeros(M,N);
    %计算雅克比矩阵
    temp_m2=m;
    temp_m1=m;
    for k=1:N
        if abs(m(k))>eps
            delta=0.01*abs(m(k));
        else
            delta=0.0001;
        end
        temp_m1(k)=m(k)-delta;
        temp_m2(k)=m(k)+delta;
        J(:,k)=(f(temp_m2,x_input)-f(temp_m1,x_input))/(2*delta);
        temp_m1(k)=m(k);
        temp_m2(k)=m(k);
    end
    
    %解dm
%     [U,S,V]=svd(J);
%     ind=find(abs(S)>1e-14);
%     S(ind)=1./S(ind);
%     S=S';
%     J_plus=V*S*U';
%     dm=J_plus*(dy(:));

    %阻尼最小二乘
    G=J.'*J;

%     whos G
%     
%     whos J
%     whos dy
    dm=(G+lambda*eye(N))\(J.')*(dy(:));
    
    m_old=m;
    m=m+dm;
    %范围
    m(m>m_upper)=m_upper(m>m_upper);
    m(m<m_lower)=m_lower(m<m_lower);    
    
    %自适应调整阻尼因子
    y_predict_old=f(m_old,x_input);
    dy_old=y_obs-y_predict_old;
    dy_old=dy_old(:);
    TR=dy_old'*J*dm;
    DQ=dm'*(J')*J*dm;

    y_predict_new=f(m,x_input);
    dy_new=y_predict_new-y_obs;
    dy_new=dy_new(:);
    DS=(dy_old)'*(dy_old)-...
        (dy_new)'*(dy_new);
    
    R_factor=DS/(2*TR-DQ);
    v=2-DS/TR;
    if v>10
        v=10;
    elseif v<2
        v=2;
    end
    
    if R_factor<R_low
        lambda=lambda*v;
    elseif R_factor>R_up
        lambda=lambda/v;
    end
    if lambda<1e-3
        lambda=1e-3;
    end
           
    

    
    dy=y_obs-f(m,x_input);
    if norm(dm,2)/max(norm(m,2),norm(m_old,2))<tol
%     if norm(dy,2)/norm(y_obs,2)<tol    
        %y_predict=f(m,x_input);
        misfit=norm(y_obs-y_predict_new,2)/norm(y_obs);        
        %fprintf(1,'迭代收敛！  迭代次数：%d \r\nrms=%e\n',[ite,misfit]);
        break;
    end
end
if ite==maxit
    misfit=norm(y_obs-y_predict_new,2)/norm(y_obs);
end

% while (norm(dy,2)/norm(y_obs,2)>1e-5)&&(number<maxit)
%     number=number+1;%计数器
%     temp_m2=m0;
%     temp_m1=m0;
%     J=zeros(M,N);
%     for k=1:N
%         if abs(m0(k))>eps
%             delta=0.01*abs(m0(k));
%         else
%             delta=0.0001;
%         end
%         temp_m1(k)=m0(k)-delta;
%         temp_m2(k)=m0(k)+delta;
%         J(:,k)=(f(temp_m2,x_input)-f(temp_m1,x_input))/(2*delta);
%         temp_m1(k)=m0(k);
%         temp_m2(k)=m0(k);
%     end
%     dm=(J.'*J)\(J.')*(dy(:));
%     m0=m0+dm;
%     dy=y_obs-f(m0,x_input);
% end
% m=m0;



%fprintf(1,'残差 %f\n\r',e);
