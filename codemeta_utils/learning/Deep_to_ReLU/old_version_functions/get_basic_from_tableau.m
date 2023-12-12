function basic = get_basic_from_tableau(A)
A = round(A,5);
[r,c] = find(A(1:end-1,1:end-1)==1);
[~,idx] = sort(r);
basic = c(idx)';
end