function [ x ] = Gauss_solve( A,b )
	n = size(A,1);
    x = zeros(n,1);
    %寻找当前列的绝对值最大的元素
	for k = 1:n-1
        Max = abs(A(k,k));
        MaxIndex = k;
		for u = k+1:n
            if(abs(A(u,k)) > Max)
                MaxIndex = u;
                Max = abs(A(u,k));
            end
        end
        
        %交换增广矩阵相应的行和列
        temp = A(MaxIndex,:);
		A(MaxIndex,:) = A(k,:);
		A(k,:) = temp;
        
        bt = b(MaxIndex);
        b(MaxIndex) = b(k);
        b(k) = bt;
        
        %判断顺序余子式是否为0，若是则方程无法用Gauss求解
        Det = A(1:k,1:k);
        if(Det==0)
            error('This matrix can''t be solved by Gauss algorithm');
        end
        
        %Gauss消元法的消元计算
		for i = k+1:n
			Mik = A(i,k)/A(k,k);
			b(i) = b(i) - Mik*b(k);
			for j = k+1:n
				A(i,j) = A(i,j) - Mik*A(k,j);
            end
        end
    end
    
    %Gauss回代运算
	x(n) = b(n) / A(n,n);
	for i = n-1:-1:1
		x(i) = ( b(i) - sum(A(i,i+1:n)*x(i+1:n)) ) / A(i,i);
    end
end