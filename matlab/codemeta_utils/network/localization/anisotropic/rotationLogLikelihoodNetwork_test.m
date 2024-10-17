function rotationLogLikelihoodNetwork_test
t_node=buildTestTagNetwork();
[Ri0,Ti0]=testNetworkGetRotTransl(t_node);
[Rij,Tij]=testNetworkGetRelativeRotTranslScales(t_node);
[GammaijR,GammaijT]=testNetworkGetDispersionMatricesRotationTranslation(t_node);
E=testNetworkGetEdges(t_node);

[Rit,dRit,Ri0,vi0,vi]=rot_randGeodFun(Ri0);

f=@(t) cost(Rit(t),Rij,GammaijR,E,vi);

check_der(f,'function')

function [c,dc]=cost(Ri,Rij,GammaijR,E,vi)
[c,Dc]=rotationLogLikelihoodNetwork(Ri,Rij,GammaijR,E);
dc=sum(sum(vi.*Dc));
