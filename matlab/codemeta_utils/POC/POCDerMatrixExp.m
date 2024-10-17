function POCDerMatrixExp
p=0.7;
[F,dF]=real_randGeodFun(randn(5,5));
A=@(t) F(t)*F(t)';
dA=@(t) dF(t)*F(t)'+F(t)*dF(t)';
% funCheckDer(A,dA)
f=@(t) trace(A(t)^p);
% this DOES NOT WORK
df=@(t) trace(dA(t)*A(t)^(p-1));
funCheckDer(f,df)


