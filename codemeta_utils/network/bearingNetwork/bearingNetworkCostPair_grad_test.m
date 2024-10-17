function bearingNetworkCostPair_grad_test
global funs
funs=bearingCostFunctions('angleSq');
[x,dx]=real_randGeodFun(randn(2));
E=[1 2];
yg=bearingNetworkComputeBearings(x(0),E);

%funCheckDer(@(t) funAndDer(x(t),E,yg,dx(t)))
%funCheckDer(@(t) derAndDder(x(t),E,yg,dx(t)))
funCheckDer(@(t) gradAndDgrad(x(t),E,yg,dx(t)))


function [f,df]=funAndDer(x,E,yg,dx)
global funs

[y,ny]=bearingNetworkComputeBearings(x,E);
f=bearingNetworkCostPair(y,yg,ny,funs);
gradf=bearingNetworkCostPair_grad(y,yg,funs,ny);
df=gradf(:)'*dx(:);


function [df,ddf]=derAndDder(x,E,yg,dx)
global funs

[y,ny]=bearingNetworkComputeBearings(x,E);
[gradf,DgradPhi]=bearingNetworkCostPair_grad(y,yg,funs,ny);
df=gradf(:)'*dx(:);
ddf=dx(:)'*DgradPhi*dx(:);

function [gradf,dgradf]=gradAndDgrad(x,E,yg,dx)
global funs

[y,ny]=bearingNetworkComputeBearings(x,E);
[gradf,DgradPhi]=bearingNetworkCostPair_grad(y,yg,funs,ny);
dgradf=reshape(DgradPhi*dx(:),size(gradf));



