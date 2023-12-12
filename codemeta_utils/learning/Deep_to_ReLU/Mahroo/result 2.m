function v = result(A,basic)
v = zeros(1,size(A,2)-1);
for i=1:size(basic,2)
    v(basic(i)) = A(:,basic(i))'*A(:,end);
end
end