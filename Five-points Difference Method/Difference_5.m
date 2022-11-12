%% Gauss_Seidel法和列主元消去法求解五点差分方程KU=F
%将二维问题当成一维问题求解，即二维的(m-1)*(n-1)个节点按顺序分行放到(m-1)*(n-1)维向量U中

clc
n=32;       %y方向的网格数
m=64;       %x方向的网格数
h2=1/n;     %y方向步长
h1=2/m;     %x方向步长
F=zeros((m-1)*(n-1),1); %KU=F，载荷向量
K=zeros((m-1)*(n-1));   %K为(m-1)*(n-1)方阵，系数矩阵
U=ones((m-1)*(n-1),1);  %U为(m-1)*(n-1)维向量,初始迭代矩阵
alpha=2*(1/h1^2+1/h2^2);
beta=1/h1^2;    gamma=1/h2^2;

%% 下面对K赋值
B= (-beta) * diag(ones(1,m-2),1)+(-beta) * diag(ones(1,m-2,1),-1)+alpha * eye(m-1);
A=-gamma * eye(m-1);
K=kron(diag(ones(1,n-2),1)+diag(ones(1,n-2),-1),A)+kron(diag(ones(1,n-1)),B);

%% 下面对F赋值
%先对非边值的节点基本赋值
for i=1:m-1
    for j=1:n-1
        F(i + (j-1)*(m-1))=(pi^2 - 1) * exp(i * h1) * sin(pi * j * h2);
    end
end
%对特殊节点赋值
for j=1:n-1
   F(1 + (j-1)*(m-1))=F(1 + (j-1)*(m-1))+beta * sin(pi * j * h2);                   %左边界
   F(m-1 + (j-1)*(m-1))=F(m-1 + (j-1)*(m-1)) +beta * exp(1)^2 * sin(pi * j * h2);    %右边界
end

%% 使用高斯赛德迭代求解KU=F
epsn=0.5 * 1e-5;  %误差限，可根据需要更改，本次我们使用1/2 * 1e-5
N=5000;    %最大迭代次数
u=Gauss_Seidel(K,F,U,epsn,N);  %将求解结果存到向量u中
u2=Gauss_solve(K,F);

%% 下面将解u放到矩阵U1中，将u2放入矩阵U3中，并与真解U2作比较
U1=zeros(m+1,n+1);  %U1用于存放向量u中元素二维化后的数据,上下边界为0
U2=zeros(m+1,n+1);  %U2用于存放真解
U3=zeros(m+1,n+1);  %U3用于存放向量u2

for i=2:m     %将数值解u的一维数据放到二维矩阵U1中
    for j=2:n
        U1(i,j)=u((i-1)+(j-2) * (m-1));
    end
end
for j=1:n+1
   U1(1,j)=sin(pi * (j-1) * h2);                 %左边界
   U1(m+1,j)=exp(1)^2 * sin(pi * (j-1) * h2);    %右边界
end
for i=1:m+1%真解产生的节点函数值放到U2中
    for j=1:n+1
        U2(i,j)=exp((i-1) * h1) * sin(pi * (j-1) * h2);
    end
end
for i=2:m     %将数值解u的一维数据放到二维矩阵U3中
    for j=2:n
        U3(i,j)=u2((i-1)+(j-2) * (m-1));
    end
end
for j=1:n+1
   U3(1,j)=sin(pi * (j-1) * h2);                 %左边界
   U3(m+1,j)=exp(1)^2 * sin(pi * (j-1) * h2);    %右边界
end

%% 画图
[X,Y]=meshgrid(0:h2:1,0:h1:2);
figure;
surf(X,Y,U2);%当要画数值解的图时，U2改为U1即可
title('理论解图像');%当要画数值解的图时，'理论解图像'改为'数值解图像'即可
hold on;
figure;
mesh(X,Y,U1);
title('Guass-Seidel解图像');
figure;
surf(X,Y,abs(U1-U2));
title('G-S误差');
figure;
surf(X,Y,U3);
title('列主元消去法解图像');
figure;
surf(X,Y,abs(U3-U2));
title('消去法误差');