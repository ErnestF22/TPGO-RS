function y=homogeneousProjectPoints(x,bOrth)
D=size(bOrth,1)-1;
DOrth=size(bOrth,2);

eDp1=[zeros(D,1);1];
eDOrthp1=[zeros(DOrth,1);1];

A=[bOrth eDp1];
PD=[eye(D) zeros(D,1)];
NX=size(x,2);
ut=ones(1,NX);
x=[x;ut];

y=PD*(x-A*((A'*A)\(A'*x-eDOrthp1*ut)));
