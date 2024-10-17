function [C,s]=rangeIntersection(A,B,tol)
m=size(A,2);
n=size(B,2);
if ~exist('tol','var')
    tol= max(m,n) * eps(class(A));
end
A=orth(A);
B=orth(B);
[U,S]=svd(A'*B,'econ');
%angles=acos(min(diag(S),1));
if m > 1 && n>1
    s = diag(S);
elseif m == 1 || n==1
    s = S(1);
else
    s = 0;
end
dimIntersection=sum((1-s)<tol);
C=A*U(:,1:dimIntersection);
