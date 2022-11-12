%% 经典显格式求解波动方程

%网格参数
K=100;      %时间网格数
J=50;      %空间网格数
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
%% 计算每层
%先算第二层（三层格式）
for j=1:J-1
    U(j+1,2)=r^2/2 * (exp((j-1)*h) + exp((j+1)*h)) + (1-r^2+tau) * exp(j*h) ;
end
%从第三层开始
for k=2:K    
    for j=2:J
        U(j,k+1)=r^2*(U(j-1,k)+U(j+1,k)) + 2*(1-r^2)*U(j,k) - U(j,k-1);
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
title('经典显格式图像');
%误差
error=abs(U0-U);
figure;
surf(X,Y,error);
title('误差图像');
