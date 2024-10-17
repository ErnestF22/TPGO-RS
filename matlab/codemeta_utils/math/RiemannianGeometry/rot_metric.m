function d=rot_metric(R,A,B)

if (isa(A, 'function_handle'))
    A = A(R);
end

if (isa(B, 'function_handle'))
    B = B(R);
end

N=size(A,3);
if(size(B,3)~=N)
    error('A and B must contain the same number of tangent vectors')
end
if(N==1)
    d=trace(A'*B)/2;
else
    d=zeros(1,N);
    for n=1:N
        d(n)=trace(A(:,:,n)'*B(:,:,n))/2;
    end
end
