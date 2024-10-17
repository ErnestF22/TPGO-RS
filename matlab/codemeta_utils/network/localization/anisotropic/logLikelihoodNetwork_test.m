function logLikelihoodNetwork_test

t_node=buildTestTagNetwork();
[Ri0,Ti0]=testNetworkGetRotTransl(t_node);
[Rij,Tij]=testNetworkGetRelativeRotTranslScales(t_node);
Gammaij=testNetworkGetDispersionMatricesRotationTranslation(t_node);
E=testNetworkGetEdges(t_node);

[Rit,dRit,Ri0,vi0,vi]=rot_randGeodFun(Ri0);
[Tit,dTit,Ti0,dTi]=real_randGeodFun(Ti0);

v=[vi;dTi];

f=@(t) cost(Rit(t),Tit(t),Rij,Tij,Gammaij,E,v);

check_der(f,'function')

function [c,dc]=cost(Ri,Ti,Rij,Tij,Gammaij,E,v)
[c,Dc]=logLikelihoodNetwork(Ri,Ti,Rij,Tij,Gammaij,E);
dc=Dc(:)'*v(:);
