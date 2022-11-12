function x=Gauss_Seidel(A,b,x0,epsn,N)

%����Guass-seidel�����������Է�����Ax=b
%epsnΪ����1e-3,������������Ĭ��500����x������ֵ������

n=length(b);
x=zeros(n,1);%��x��ֵ
k=1;
while k<N %����������
    for i=1:n
        if i==1
            x(1)=(b(1)-A(1,2:n)*x0(2:n))/A(1,1);
        elseif i==n
            x(n)=(b(n)-A(n,1:n-1)*x(1:n-1))/A(n,n);
        else 
            x(i)=(b(i)-A(i,1:i-1)*x(1:i-1)-A(i,i+1:n)*x0(i+1:n))/A(i,i);
        end
    end
    if norm(x-x0,inf)<epsn   %����޷�������޾�ֹͣѭ��
        break;
    end
    
    x0=x;
    k=k+1;
end
    if k==N
        disp('���������Ѵ�����!');
    end
    disp(['��������k =',num2str(k)])
