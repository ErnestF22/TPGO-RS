function POCReviewDerivativeSchattenpNorm
paper='iccv';
%paper='aaai';
p=0.7;

switch paper
    case 'iccv'
        X=randn(5,10);
        %[F,dF]=real_randGeodFun(randn(10,10));
        %G=@(t) F(t)*F(t)';
        %dG=@(t) dF(t)*F(t)'+F(t)*dF(t)';
        [d,dd]=real_geodFun(rand(10,1),rand(10,1));
        G=@(t) diag(d(t));
        dG=@(t) diag(dd(t));
        %funCheckDer(G,dG)
        f=@(G) schatten(X*G,p)^2;
        D=@(G) p*schatten(X*G,p)^2*(X*G^2*X')^((p-2)/2);
        gradf=@(G) 2*X'*D(G)*X*G;
        df=@(G,dG) trace(gradf(G)'*dG);
        %funCheckDer(@(t) f(G(t)), @(t) df(G(t),dG(t)));
        A=@(t) X*G(t)'*G(t)*X';
        dA=@(t) X*dG(t)'*G(t)*X'+X*G(t)'*dG(t)*X';
        funCheckDer(@(t) funAndDer(A(t),dA(t),p))        
        
    case 'aaai'
        [X,dX]=real_randGeodFun(randn(10,5));
        A=@(t) X(t)'*X(t);
        dA=@(t) dX(t)'*X(t)+X(t)'*dX(t);
        f=@(X) trace((X'*X)^(p/2));
        D=@(X) p/2*(X'*X)^((p-2)/2);
        gradf=@(X) 2*X*D(X);
        df=@(X,dX) trace(gradf(X)'*X);
        %funCheckDer(@(t) f(X(t)), @(t) df(X(t),dX(t)))
        %funCheckDer(@(t) sAndDs(A(t),dA(t)))
        funCompare(@(t) f(X(t)), @(t) funAndDer(A(t),dA(t),p))
        funCheckDer(@(t) funAndDer(A(t),dA(t),p))        
end

function [f,df]=funAndDer(A,dA,p)
[s,ds,U]=sAndDs(A,dA);
f=trace(U*diag(s.^(p/2))*U');
df=trace(U*diag(p/2*s.^((p-2)/2))*diag(ds)*U');


function [s,ds,U,V]=sAndDs(A,dA)
[U,S,V]=svd(A,'econ');
s=diag(S);
ds=diag(U'*dA*U);
