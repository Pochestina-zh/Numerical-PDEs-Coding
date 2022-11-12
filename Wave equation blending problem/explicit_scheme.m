%% �����Ը�ʽ��Ⲩ������

%�������
K=100;      %ʱ��������
J=50;      %�ռ�������
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
%% ����ÿ��
%����ڶ��㣨�����ʽ��
for j=1:J-1
    U(j+1,2)=r^2/2 * (exp((j-1)*h) + exp((j+1)*h)) + (1-r^2+tau) * exp(j*h) ;
end
%�ӵ����㿪ʼ
for k=2:K    
    for j=2:J
        U(j,k+1)=r^2*(U(j-1,k)+U(j+1,k)) + 2*(1-r^2)*U(j,k) - U(j,k-1);
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
title('�����Ը�ʽͼ��');
%���
error=abs(U0-U);
figure;
surf(X,Y,error);
title('���ͼ��');
