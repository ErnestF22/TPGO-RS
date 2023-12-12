function stiefel_min_test
randn('state',0);
n=5;
p=3;

%A=randn(n,p);
A=orth(randn(n,p));
Y=orth(randn(n,p));

f=@(Y) 0.5*sum((Y(:)-A(:)).^2);

epsilon=0.1;
Nit=100;

for it=1:Nit
    d(it)=f(Y);
    gradf=stiefel_tangentProj(Y,Y-A);
    Y=stiefel_exp(Y,-epsilon*gradf);
end

plot(d)
[A Y]
