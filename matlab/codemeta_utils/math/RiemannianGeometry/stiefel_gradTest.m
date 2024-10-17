function stiefel_gradTest
P=randn(4,3);
Y0=eye(4,3);

f=@(x) trace(x'*P);
fGrad=@(x) P;

[xt,dxt]=real_randGeodFun(Y0);
ft=@(t) f(xt(t));
df=@(t) stiefel_metric(xt(t),fGrad(xt(t)),dxt(t));

funCheckDer(ft,df)