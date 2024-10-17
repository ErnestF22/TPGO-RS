function bearingNetworkCostGradient_test
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
df=@(t) derAndDder(xt(t),E,Yijtruth,v);
t=linspace(-1,1,50);

figure(1)
%check_der(f,'function',t)
check_der(df,'function',t)

function [c,dc]=costAndDer(x,E,Yijtruth,v)
global funs
[Yij,nYij]=bearingNetworkComputeBearings(x,E);
c=bearingNetworkCost(E,Yij,Yijtruth,nYij,funs);
gradc=bearingNetworkCostGradient(E,Yij,Yijtruth,funs);
dc=v(:)'*gradc(:);

function [dc,ddc]=derAndDder(x,E,Yijtruth,v)
global funs
[Yij,nYij]=bearingNetworkComputeBearings(x,E);
[gradc,Dgradc]=bearingNetworkCostGradient(E,Yij,Yijtruth,funs,nYij);
dc=v(:)'*gradc(:);
ddc=v(:)'*Dgradc*v(:);
