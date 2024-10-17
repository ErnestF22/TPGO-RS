function POCOrthogonalityGradient
global funs

NNodes=11;
t_node=bearingNetworkBuildTestNetwork(NNodes);
funs=bearingCostFunctions('angleSq');

x0=t_node.Ti;
E=t_node.E;
yg=t_node.Yijtruth;

v=x0;
%v=randn(size(x0));
x=real_geodFun(x0,v);

funCheckDer(@(t) costAndDer(x(t),E,yg,v),'function',linspace(-0.5,0.5));


function [c,dc]=costAndDer(x,E,yg,v)
global funs
[y,ny]=bearingNetworkComputeBearings(x,E);
c=bearingNetworkCost(E,y,yg,ny,funs);
gradc=bearingNetworkCostGradient(E,y,yg,funs);
dc=v(:)'*gradc(:);
