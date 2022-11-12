%% 经典隐格式求解波动方程

%网格参数
K=80;      %时间网格数
J=50;       %空间网格数
h=1/J;      %空间步长
tau=1/K;    %时间步长
r=tau/h;    %网格比
%矩阵参数
U=zeros(J+1,K+1);   %表示数值解
U0=zeros(J+1,K+1);  %表示精确解
error=zeros(J+1,K+1);%表示解的误差
%% 先存初边值
%初始条件
for j=1:J+1
    U(j,1)=exp((j-1) * h);               %t=0，初始下边界
end
%边界条件
for k=1:K+1
    U(1,k)=exp((k-1) * tau);    %左边界
    U(J+1,k)=exp(1+((k-1) * tau));  %右边界
end
%% 先算第二层
for j=1:J-1
    U(j+1,2)=r^2/2 * (exp((j-1)*h) + exp((j+1)*h)) + (1-r^2+tau) * exp(j*h) ;
end
%% 再用Gauss列主元消去法解方程组（隐式三层）AU_k+1=F，k=2,3,...(从第三层开始算)
F=zeros(J-1,1);   %隐式右端向量
A=zeros(J-1,K-1); %隐式系数矩阵
u=zeros(J-1,1);   %隐式解

%生成隐式系数矩阵
A=(1+r^2/2) * diag(ones(J-1,1)) - r^2/4 * diag(ones(J-2,1),1) - r^2/4 * diag(ones(J-2,1),-1);
for k=2:K   %从第三层开始算
    %生成隐式右端向量
    for j=1:J-1
       F(j)= r^2/2*(U(j+2,k)+U(j,k)) + (2-r^2)*U(j+1,k) + r^2/4*(U(j+2,k-1)+U(j,k-1)) - (1+r^2/2)*U(j+1,k-1);
    end
    F(1)=F(1)+1/4 * r^2 * U(1,k+1);
    F(J-1)=F(J-1)+1/4 * r^2 * U(J+1,k+1);
    %Gauss消去法求解
    u=Gauss_solve(A,F);
    %将隐式解放入数值解中
    for j=2:J
        U(j,k+1)=u(j-1);
    end
end
%% 画图 
%精确解
for j=1:J+1
    for k=1:K+1
        U0(j,k)=exp((j-1)*h + (k-1)*tau);
    end
end
%画图
[X,Y]=meshgrid(0:tau:1,0:h:1);
figure;
surf(X,Y,U0);
title('精确解的图像');
hold on;
figure;
mesh(X,Y,U);
hold on;
title('经典隐格式图像');
%误差
error=abs(U0-U);
figure;
surf(X,Y,error);
title('误差图像');