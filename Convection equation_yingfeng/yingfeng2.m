%% 迎风格式求解传输方程（a=-1,显示左偏心）
clear;
clc;

%网格参数
J=1000;       %空间网格数,范围x=[0,2pi]
K=5*J;      %时间网格数,时间范围为T=[0,pi]
h=2*pi/J;   %空间步长
tau=0.1*h;  %时间步长
r=0.1;      %网格比

%矩阵（解）参数
U0=zeros(J+1,K+1);  %表示精确解
U=zeros(5*J+1,K+1);   %表示增广数值解，范围[-8pi,2pi]
U1=zeros(J+1,K+1);   %表示数值解
error=zeros(J+1,K+1);%表示解的误差

%% 先存初边值
%初始条件
for j=1:J+1
    U(j,1)=sin(2*pi-(j-1)*h);               %t=0，初始下边界[2pi,0]
end
%周期性边界条件
for i=1:4
    for j=1:J
        U(j+i*J,1)=U(j,1);                  %t=0的初始下边界[0,-8pi]
    end
end
U(5*J+1,1)=U(1,1);  %x=-8pi处的初始条件

%% 计算每层

%从第二层开始
for k=2:K+1
    for j=k:5*J+1
        U(j,k)=(1-r) * U(j,k-1) + r * U(j-1,k-1);
    end
end

%利用周期性边界循环出增广数值解[2pi,0]
for k=1:K+1        
    U(1,k)=U(5*J+1,k); %平移边界x=-8pi到x=2pi
end

for k=2:K+1
    for j=2:J+1
        U(j,k)=(1-r) * U(j,k-1) + r * U(j-1,k-1);
    end
end

%由增广解得到数值解
for j=1:J+1
    for k=1:K+1
        U1(j,k)=U(j,k);         
    end
end

%% 画图
%精确解
for j=1:J+1
    for k=1:K+1
        U0(j,k)=sin(2*pi-(j-1)*h+(k-1)*tau);
    end
end
%画图
[X,Y]=meshgrid(0:tau:pi,0:h:2*pi);
figure;
mesh(X,Y,U0);
title('精确解的图像（a=-1）');
x1=xlabel('t轴(时间)');        %x轴标题
x2=ylabel('x轴(空间)');        %y轴标题
x3=zlabel('方程解u(x,t)');        %z轴标题
set(x1,'Rotation',30);    %x轴名称旋转
set(x2,'Rotation',-30);    %y轴名称旋转
hold on;
figure;
mesh(X,Y,U1);
hold on;
title('迎风格式（a=-1）数值解图像');
x1=xlabel('t轴(时间)');        %x轴标题
x2=ylabel('x轴(空间)');        %y轴标题
x3=zlabel('方程解u(x,t)');        %z轴标题
set(x1,'Rotation',30);    %x轴名称旋转
set(x2,'Rotation',-30);    %y轴名称旋转
%误差
error=abs(U0-U1);
err=max(max(U0-U1));
figure;
mesh(X,Y,error);
title('误差图像（a=-1）');
x1=xlabel('t轴(时间)');        %x轴标题
x2=ylabel('x轴(空间)');        %y轴标题
x3=zlabel('方程解u(x,t)');        %z轴标题
set(x1,'Rotation',30);    %x轴名称旋转
set(x2,'Rotation',-30);    %y轴名称旋转
