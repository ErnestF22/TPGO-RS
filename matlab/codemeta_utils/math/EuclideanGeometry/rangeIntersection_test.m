function rangeIntersection_test
A=randn(5,4);
B=randn(5,3);

n=null([A B]);
S1=A*n(1:size(A,2),:);
S2=rangeIntersection(A,B,1e-6);

disp(size(S1))
disp(size(S2))

subspace(S1,S2)
subspace(S1,A)
subspace(S1,B)
subspace(S2,A)
subspace(S2,B)
