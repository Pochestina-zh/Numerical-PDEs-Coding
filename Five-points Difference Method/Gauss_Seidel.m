function x=Gauss_Seidel(A,b,x0,epsn,N)

%用于Guass-seidel迭代法解线性方程组Ax=b
%epsn为精度1e-3,最大迭代次数（默认500），x返回数值解向量

n=length(b);
x=zeros(n,1);%给x赋值
k=1;
while k<N %最大迭代次数
    for i=1:n
        if i==1
            x(1)=(b(1)-A(1,2:n)*x0(2:n))/A(1,1);
        elseif i==n
            x(n)=(b(n)-A(n,1:n-1)*x(1:n-1))/A(n,n);
        else 
            x(i)=(b(i)-A(i,1:i-1)*x(1:i-1)-A(i,i+1:n)*x0(i+1:n))/A(i,i);
        end
    end
    if norm(x-x0,inf)<epsn   %满足∞范数误差限就停止循环
        break;
    end
    
    x0=x;
    k=k+1;
end
    if k==N
        disp('迭代次数已达上限!');
    end
    disp(['迭代次数k =',num2str(k)])
