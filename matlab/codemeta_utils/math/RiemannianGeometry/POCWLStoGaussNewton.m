function POCWLStoGaussNewton
preG=randn(5);
G=preG*preG';
R=chol(G);

r=@(x) sqrt(x'*G*x);

f=@(x) r(x)^2;

dr=@(x) (G*x)/r(x);

[tx,dtx]=real_randGeodFun(randn(5,1));

tr=@(t) r(tx(t));
dr=@(t) dr(tx(t))'*dtx(t);

check_der(tr,dr)
