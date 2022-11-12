%% ��������ʽ��Ⲩ������

%�������
K=80;      %ʱ��������
J=50;       %�ռ�������
h=1/J;      %�ռ䲽��
tau=1/K;    %ʱ�䲽��
r=tau/h;    %�����
%�������
U=zeros(J+1,K+1);   %��ʾ��ֵ��
U0=zeros(J+1,K+1);  %��ʾ��ȷ��
error=zeros(J+1,K+1);%��ʾ������
%% �ȴ����ֵ
%��ʼ����
for j=1:J+1
    U(j,1)=exp((j-1) * h);               %t=0����ʼ�±߽�
end
%�߽�����
for k=1:K+1
    U(1,k)=exp((k-1) * tau);    %��߽�
    U(J+1,k)=exp(1+((k-1) * tau));  %�ұ߽�
end
%% ����ڶ���
for j=1:J-1
    U(j+1,2)=r^2/2 * (exp((j-1)*h) + exp((j+1)*h)) + (1-r^2+tau) * exp(j*h) ;
end
%% ����Gauss����Ԫ��ȥ���ⷽ���飨��ʽ���㣩AU_k+1=F��k=2,3,...(�ӵ����㿪ʼ��)
F=zeros(J-1,1);   %��ʽ�Ҷ�����
A=zeros(J-1,K-1); %��ʽϵ������
u=zeros(J-1,1);   %��ʽ��

%������ʽϵ������
A=(1+r^2/2) * diag(ones(J-1,1)) - r^2/4 * diag(ones(J-2,1),1) - r^2/4 * diag(ones(J-2,1),-1);
for k=2:K   %�ӵ����㿪ʼ��
    %������ʽ�Ҷ�����
    for j=1:J-1
       F(j)= r^2/2*(U(j+2,k)+U(j,k)) + (2-r^2)*U(j+1,k) + r^2/4*(U(j+2,k-1)+U(j,k-1)) - (1+r^2/2)*U(j+1,k-1);
    end
    F(1)=F(1)+1/4 * r^2 * U(1,k+1);
    F(J-1)=F(J-1)+1/4 * r^2 * U(J+1,k+1);
    %Gauss��ȥ�����
    u=Gauss_solve(A,F);
    %����ʽ�������ֵ����
    for j=2:J
        U(j,k+1)=u(j-1);
    end
end
%% ��ͼ 
%��ȷ��
for j=1:J+1
    for k=1:K+1
        U0(j,k)=exp((j-1)*h + (k-1)*tau);
    end
end
%��ͼ
[X,Y]=meshgrid(0:tau:1,0:h:1);
figure;
surf(X,Y,U0);
title('��ȷ���ͼ��');
hold on;
figure;
mesh(X,Y,U);
hold on;
title('��������ʽͼ��');
%���
error=abs(U0-U);
figure;
surf(X,Y,error);
title('���ͼ��');