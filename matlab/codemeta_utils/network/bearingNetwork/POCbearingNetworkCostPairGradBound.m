function POCbearingNetworkCostPairGradBound
global funs
funs=bearingCostFunctions('angleSq');
[x,dx]=real_randGeodFun(5*randn(2));
E=[1 2];
yg=bearingNetworkComputeBearings(x(0),E);

t=linspace(0,1);
%funCheckDer(@(t) funAndDer(x(t),E,yg,dx(t)),'function',t)
funPlot(@(t) funAndNormGradSq(x(t),E,yg),t)

function [f,df]=funAndDer(x,E,yg,dx)
global funs

[y,ny]=bearingNetworkComputeBearings(x,E);
f=bearingNetworkCostPair(y,yg,ny,funs);
gradf=bearingNetworkCostPair_grad(y,yg,funs,ny);
df=gradf(:)'*dx(:);

function h=funAndNormGradSq(x,E,yg)
global funs

[y,ny]=bearingNetworkComputeBearings(x,E);
f=bearingNetworkCostPair(y,yg,ny,funs);
gradf=bearingNetworkCostPair_grad(y,yg,funs,ny);

h=[f;norm(gradf)^2];

