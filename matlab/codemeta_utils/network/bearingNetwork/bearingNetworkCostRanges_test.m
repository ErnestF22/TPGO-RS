function bearingNetworkCostRanges_test
resetRands()

global funs;
costName='squared';

funs=bearingCostFunctions(costName);

t_node=bearingNetworkBuildTestNetwork();
[xt,~,x0,v]=real_randGeodFun(t_node.Titruth);
E=t_node.E;
[yg,nyg]=bearingNetworkComputeBearings(x0,E);

f=@(t) costAndDer(xt(t),E,yg,nyg,v);
t=linspace(0,1,50);

figure(1)
check_der(f,'function',t)

figure(2)
d1=@(t) der1(xt(t),E,yg,nyg,v);
d2=@(t) der2(xt(t),E,yg,nyg,v);

plotfun(d1,t)
hold on
plotfun(d2,t,'rx')
hold off

function [c,dc]=costAndDer(x,E,yg,nyg,v)
global funs
[y,ny]=bearingNetworkComputeBearings(x,E);
[c,gradc]=bearingNetworkCostRanges(E,y,yg,ny,nyg,funs);
dc=v(:)'*gradc(:);

function dc=der1(x,E,yg,nyg,v)
global funs
[y,ny]=bearingNetworkComputeBearings(x,E);
[~,gradc]=bearingNetworkCostRanges(E,y,yg,ny,nyg,funs);
dc=v(:)'*gradc(:);

function dc=der2(x,E,yg,nyg,v)
global funs
[y,ny]=bearingNetworkComputeBearings(x,E);
gradc=bearingNetworkCostRangesGradient(E,y,yg,ny,nyg,funs);
dc=v(:)'*gradc(:);

