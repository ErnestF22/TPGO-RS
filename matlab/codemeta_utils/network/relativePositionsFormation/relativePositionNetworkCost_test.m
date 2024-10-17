function relativePositionNetworkCost_test
t_node=relativePositionNetworkTestNetwork(7);
[Tt,dTt]=real_randGeodFun(t_node.TiTruth);
E=t_node.E;
TijTruth=t_node.TijTruth;

funCheckDer(@(t) costAndDer(Tt(t),E,TijTruth,dTt(t)));

function [cost,dCost]=costAndDer(Ti,E,TijTruth,v)
Tij=relativePositionNetworkCompute(Ti,E);
[cost,gradCost]=relativePositionNetworkCost(E,Tij,TijTruth);
dCost=gradCost(:)'*v(:);



