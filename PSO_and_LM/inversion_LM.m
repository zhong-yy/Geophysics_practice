function [m,misfit,ite]=inversion_LM(f,y_obs,x_input,m0,m_lower,m_upper,maxit,tol,lambda,R_low,R_up)
% =========================================================================
% LM�����������������С���˷�
% ����Դ
% 2017, 11
% ���ϴ�ѧ�������ѧ����Ϣ����ѧԺ
% ==============================����=======================================
% f: Ϊ�����⺯�������f�ĵ�����ʽ����Ϊf(m0,x_input)
% y_obs: Ϊ�۲�ֵ����֪�������
% x_inputΪf����֪���룬
% mΪ��Ҫ���ݵĲ���
% 
% m0: Ϊ�����ĳ�ֵ
% m_lower: ģ�Ͳ����½�
% m_upper: ģ�Ͳ����Ͻ�
% maxit: Ϊ��������;
% tol: Ϊ���ε����õ��Ĳ����������������ε������ý���ܽӽ�����ֹ����
% lambda: ��ʼ��������
% R_low: ���ڵ����������ӵı������ӵ��½磬һ��Ϊ0.25
% R_up: ���ڵ����������ӵı������ӵ��Ͻ磬һ��Ϊ0.75
% ==============================����ֵ=====================================
% m: �ﵽ������������ﵽ�������Ⱥ󣬵õ���ģ�Ͳ�������
% misfit: ���ƫ��
% ite: ����������
% =========================================================================
N=length(m0);
M=length(y_obs);

if M~=size(x_input,1)
    error('���������ƥ��');
end

dy=y_obs-f(m0,x_input);
%����ȷ��y_obs��f(m0,x_input)�����������������������������������

m=m0(:);
m_upper=m_upper(:);
m_lower=m_lower(:);
for ite=1:maxit
    J=zeros(M,N);
    %�����ſ˱Ⱦ���
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
    
    %��dm
%     [U,S,V]=svd(J);
%     ind=find(abs(S)>1e-14);
%     S(ind)=1./S(ind);
%     S=S';
%     J_plus=V*S*U';
%     dm=J_plus*(dy(:));

    %������С����
    G=J.'*J;

%     whos G
%     
%     whos J
%     whos dy
    dm=(G+lambda*eye(N))\(J.')*(dy(:));
    
    m_old=m;
    m=m+dm;
    %��Χ
    m(m>m_upper)=m_upper(m>m_upper);
    m(m<m_lower)=m_lower(m<m_lower);    
    
    %����Ӧ������������
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
        %fprintf(1,'����������  ����������%d \r\nrms=%e\n',[ite,misfit]);
        break;
    end
end
if ite==maxit
    misfit=norm(y_obs-y_predict_new,2)/norm(y_obs);
end

% while (norm(dy,2)/norm(y_obs,2)>1e-5)&&(number<maxit)
%     number=number+1;%������
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



%fprintf(1,'�в� %f\n\r',e);
