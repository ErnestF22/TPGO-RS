function POCDerSingularVector
[A,dA]=real_randGeodFun(randn(5));

funCheckDer(@(t) funAndDer(A(t),dA(t),1))


function [l,dl]=funAndDer(A,dA,k)
[U,S,V]=svd(A);
l=S(k,k);
dl=U(:,k)'*dA*V(:,k);
