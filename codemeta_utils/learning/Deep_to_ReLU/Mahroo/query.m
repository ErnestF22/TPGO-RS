function Q = query(A,nonbasic)
Q=[];
for i=1:size(nonbasic,2)
    n = nonbasic(i);
    m = 1:size(A,1)-1;
    m = m(A(1:end-1,n)>0);
    for j = 1:length(m)
        Q = [Q;n m(j)];
    end
end
end