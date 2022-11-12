function [ x ] = Gauss_solve( A,b )
	n = size(A,1);
    x = zeros(n,1);
    %Ѱ�ҵ�ǰ�еľ���ֵ����Ԫ��
	for k = 1:n-1
        Max = abs(A(k,k));
        MaxIndex = k;
		for u = k+1:n
            if(abs(A(u,k)) > Max)
                MaxIndex = u;
                Max = abs(A(u,k));
            end
        end
        
        %�������������Ӧ���к���
        temp = A(MaxIndex,:);
		A(MaxIndex,:) = A(k,:);
		A(k,:) = temp;
        
        bt = b(MaxIndex);
        b(MaxIndex) = b(k);
        b(k) = bt;
        
        %�ж�˳������ʽ�Ƿ�Ϊ0�������򷽳��޷���Gauss���
        Det = A(1:k,1:k);
        if(Det==0)
            error('This matrix can''t be solved by Gauss algorithm');
        end
        
        %Gauss��Ԫ������Ԫ����
		for i = k+1:n
			Mik = A(i,k)/A(k,k);
			b(i) = b(i) - Mik*b(k);
			for j = k+1:n
				A(i,j) = A(i,j) - Mik*A(k,j);
            end
        end
    end
    
    %Gauss�ش�����
	x(n) = b(n) / A(n,n);
	for i = n-1:-1:1
		x(i) = ( b(i) - sum(A(i,i+1:n)*x(i+1:n)) ) / A(i,i);
    end
end