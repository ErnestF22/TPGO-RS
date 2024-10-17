function house_test
d=5;
x=randn(d,1);
[v,beta]=house(x);
P=eye(d)-beta*(v*v');
eig(P)