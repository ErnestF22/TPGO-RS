function [R,S,flagIsRigid] = distanceRigidityMatrix(E,x)
[d,n]=size(x);
m=size(E,1);
R=zeros(m,d*n);
for Midx=1:m
    i=E(Midx,1);
    j=E(Midx,2);
    R(Midx,(i-1)*d+1:i*d)=x(:,i)-x(:,j);
    R(Midx,(j-1)*d+1:j*d)=x(:,j)-x(:,i);
end
% S is the maximum rank of R which happens when G is rigid
if n>=d+2
    S=n*d-d*(d+1)/2;
else
    S=n*(n-1)/2;
end

if rank(R)==S
    flagIsRigid=true;
else
    flagIsRigid=false;
end
end