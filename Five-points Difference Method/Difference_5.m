%% Gauss_Seidel��������Ԫ��ȥ���������ַ���KU=F
%����ά���⵱��һά������⣬����ά��(m-1)*(n-1)���ڵ㰴˳����зŵ�(m-1)*(n-1)ά����U��

clc
n=32;       %y�����������
m=64;       %x�����������
h2=1/n;     %y���򲽳�
h1=2/m;     %x���򲽳�
F=zeros((m-1)*(n-1),1); %KU=F���غ�����
K=zeros((m-1)*(n-1));   %KΪ(m-1)*(n-1)����ϵ������
U=ones((m-1)*(n-1),1);  %UΪ(m-1)*(n-1)ά����,��ʼ��������
alpha=2*(1/h1^2+1/h2^2);
beta=1/h1^2;    gamma=1/h2^2;

%% �����K��ֵ
B= (-beta) * diag(ones(1,m-2),1)+(-beta) * diag(ones(1,m-2,1),-1)+alpha * eye(m-1);
A=-gamma * eye(m-1);
K=kron(diag(ones(1,n-2),1)+diag(ones(1,n-2),-1),A)+kron(diag(ones(1,n-1)),B);

%% �����F��ֵ
%�ȶԷǱ�ֵ�Ľڵ������ֵ
for i=1:m-1
    for j=1:n-1
        F(i + (j-1)*(m-1))=(pi^2 - 1) * exp(i * h1) * sin(pi * j * h2);
    end
end
%������ڵ㸳ֵ
for j=1:n-1
   F(1 + (j-1)*(m-1))=F(1 + (j-1)*(m-1))+beta * sin(pi * j * h2);                   %��߽�
   F(m-1 + (j-1)*(m-1))=F(m-1 + (j-1)*(m-1)) +beta * exp(1)^2 * sin(pi * j * h2);    %�ұ߽�
end

%% ʹ�ø�˹���µ������KU=F
epsn=0.5 * 1e-5;  %����ޣ��ɸ�����Ҫ���ģ���������ʹ��1/2 * 1e-5
N=5000;    %����������
u=Gauss_Seidel(K,F,U,epsn,N);  %��������浽����u��
u2=Gauss_solve(K,F);

%% ���潫��u�ŵ�����U1�У���u2�������U3�У��������U2���Ƚ�
U1=zeros(m+1,n+1);  %U1���ڴ������u��Ԫ�ض�ά���������,���±߽�Ϊ0
U2=zeros(m+1,n+1);  %U2���ڴ�����
U3=zeros(m+1,n+1);  %U3���ڴ������u2

for i=2:m     %����ֵ��u��һά���ݷŵ���ά����U1��
    for j=2:n
        U1(i,j)=u((i-1)+(j-2) * (m-1));
    end
end
for j=1:n+1
   U1(1,j)=sin(pi * (j-1) * h2);                 %��߽�
   U1(m+1,j)=exp(1)^2 * sin(pi * (j-1) * h2);    %�ұ߽�
end
for i=1:m+1%�������Ľڵ㺯��ֵ�ŵ�U2��
    for j=1:n+1
        U2(i,j)=exp((i-1) * h1) * sin(pi * (j-1) * h2);
    end
end
for i=2:m     %����ֵ��u��һά���ݷŵ���ά����U3��
    for j=2:n
        U3(i,j)=u2((i-1)+(j-2) * (m-1));
    end
end
for j=1:n+1
   U3(1,j)=sin(pi * (j-1) * h2);                 %��߽�
   U3(m+1,j)=exp(1)^2 * sin(pi * (j-1) * h2);    %�ұ߽�
end

%% ��ͼ
[X,Y]=meshgrid(0:h2:1,0:h1:2);
figure;
surf(X,Y,U2);%��Ҫ����ֵ���ͼʱ��U2��ΪU1����
title('���۽�ͼ��');%��Ҫ����ֵ���ͼʱ��'���۽�ͼ��'��Ϊ'��ֵ��ͼ��'����
hold on;
figure;
mesh(X,Y,U1);
title('Guass-Seidel��ͼ��');
figure;
surf(X,Y,abs(U1-U2));
title('G-S���');
figure;
surf(X,Y,U3);
title('����Ԫ��ȥ����ͼ��');
figure;
surf(X,Y,abs(U3-U2));
title('��ȥ�����');