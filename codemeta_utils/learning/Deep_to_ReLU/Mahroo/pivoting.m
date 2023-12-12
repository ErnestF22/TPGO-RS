function [A,basic] = pivoting(A,P,basic)
n = P(1);
m = P(2);
b = zeros(size(A,1)-1,1);
b(m,1) = 1;
basicout =[];
for i=1:size(A,2)-1
    if all(A(1:end-1,i)==b)
        basicout = i;
    end
end
for i = 1:size(A,1)-1
    if i~=m
        c = A(i,n)/A(m,n);
        A(i,:) = A(i,:)-c*A(m,:);
    end
end
A(m,:) = A(m,:)/A(m,n);
if all(A(:,end)>=0)
    basic(basic==basicout) = n;
else
    basic = [];
end
end