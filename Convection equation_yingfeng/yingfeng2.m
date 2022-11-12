%% ӭ���ʽ��⴫�䷽�̣�a=-1,��ʾ��ƫ�ģ�
clear;
clc;

%�������
J=1000;       %�ռ�������,��Χx=[0,2pi]
K=5*J;      %ʱ��������,ʱ�䷶ΧΪT=[0,pi]
h=2*pi/J;   %�ռ䲽��
tau=0.1*h;  %ʱ�䲽��
r=0.1;      %�����

%���󣨽⣩����
U0=zeros(J+1,K+1);  %��ʾ��ȷ��
U=zeros(5*J+1,K+1);   %��ʾ������ֵ�⣬��Χ[-8pi,2pi]
U1=zeros(J+1,K+1);   %��ʾ��ֵ��
error=zeros(J+1,K+1);%��ʾ������

%% �ȴ����ֵ
%��ʼ����
for j=1:J+1
    U(j,1)=sin(2*pi-(j-1)*h);               %t=0����ʼ�±߽�[2pi,0]
end
%�����Ա߽�����
for i=1:4
    for j=1:J
        U(j+i*J,1)=U(j,1);                  %t=0�ĳ�ʼ�±߽�[0,-8pi]
    end
end
U(5*J+1,1)=U(1,1);  %x=-8pi���ĳ�ʼ����

%% ����ÿ��

%�ӵڶ��㿪ʼ
for k=2:K+1
    for j=k:5*J+1
        U(j,k)=(1-r) * U(j,k-1) + r * U(j-1,k-1);
    end
end

%���������Ա߽�ѭ����������ֵ��[2pi,0]
for k=1:K+1        
    U(1,k)=U(5*J+1,k); %ƽ�Ʊ߽�x=-8pi��x=2pi
end

for k=2:K+1
    for j=2:J+1
        U(j,k)=(1-r) * U(j,k-1) + r * U(j-1,k-1);
    end
end

%�������õ���ֵ��
for j=1:J+1
    for k=1:K+1
        U1(j,k)=U(j,k);         
    end
end

%% ��ͼ
%��ȷ��
for j=1:J+1
    for k=1:K+1
        U0(j,k)=sin(2*pi-(j-1)*h+(k-1)*tau);
    end
end
%��ͼ
[X,Y]=meshgrid(0:tau:pi,0:h:2*pi);
figure;
mesh(X,Y,U0);
title('��ȷ���ͼ��a=-1��');
x1=xlabel('t��(ʱ��)');        %x�����
x2=ylabel('x��(�ռ�)');        %y�����
x3=zlabel('���̽�u(x,t)');        %z�����
set(x1,'Rotation',30);    %x��������ת
set(x2,'Rotation',-30);    %y��������ת
hold on;
figure;
mesh(X,Y,U1);
hold on;
title('ӭ���ʽ��a=-1����ֵ��ͼ��');
x1=xlabel('t��(ʱ��)');        %x�����
x2=ylabel('x��(�ռ�)');        %y�����
x3=zlabel('���̽�u(x,t)');        %z�����
set(x1,'Rotation',30);    %x��������ת
set(x2,'Rotation',-30);    %y��������ת
%���
error=abs(U0-U1);
err=max(max(U0-U1));
figure;
mesh(X,Y,error);
title('���ͼ��a=-1��');
x1=xlabel('t��(ʱ��)');        %x�����
x2=ylabel('x��(�ռ�)');        %y�����
x3=zlabel('���̽�u(x,t)');        %z�����
set(x1,'Rotation',30);    %x��������ת
set(x2,'Rotation',-30);    %y��������ת
