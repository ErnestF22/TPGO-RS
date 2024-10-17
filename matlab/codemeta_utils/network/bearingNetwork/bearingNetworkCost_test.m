function bearingNetworkCost_test
resetRands()

global funs;
costName='cosine';
%costName='angleSq';

funs=bearingCostFunctions(costName);

t_node=bearingNetworkBuildTestNetwork();
[xt,~,x0,v]=real_randGeodFun(t_node.Titruth);
E=t_node.E;
Yijtruth=bearingNetworkComputeBearings(x0,E);

f=@(t) costAndDer(xt(t),E,Yijtruth,v);
t=linspace(-1,1,50);

figure(1)
check_der(f,'function',t)

figure(2)
d1=@(t) der1(xt(t),E,Yijtruth,v);
d2=@(t) der2(xt(t),E,Yijtruth,v);

plotfun(d1,t)
hold on
plotfun(d2,t,'rx')
hold off

function [c,dc]=costAndDer(x,E,Yijtruth,v)
global funs
[Yij,nYij]=bearingNetworkComputeBearings(x,E);
[c,gradc]=bearingNetworkCost(E,Yij,Yijtruth,nYij,funs);
dc=v(:)'*gradc(:);

function dc=der1(x,E,Yijtruth,v)
global funs
[Yij,nYij]=bearingNetworkComputeBearings(x,E);
[~,gradc]=bearingNetworkCost(E,Yij,Yijtruth,nYij,funs);
dc=v(:)'*gradc(:);

function dc=der2(x,E,Yijtruth,v)
global funs
Yij=bearingNetworkComputeBearings(x,E);
gradc=bearingNetworkCost_grad(E,Yij,Yijtruth,funs);
dc=v(:)'*gradc(:);

