function funCheckDer_test
A=[1 2;2 1];

f=@(x) x'*A*x;
fGrad=@(x) A*x;

%Create a parametric line x(t) in a random direction
[xt,xDott]=real_randGeodFun(randn(2,1));

ft=@(t) f(xt(t));
dft=@(t) fGrad(xt(t))'*xDott(t);

funCheckDer(ft,dft)

