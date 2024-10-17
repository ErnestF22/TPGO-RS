function rotLocFrobCostPair_test
[Ri,vi]=real_randGeodFun(rot_randn());
[Rj,vj]=real_randGeodFun(rot_randn());
Rij=Ri(0)'*Rj(0);

funCheckDer(@(t) costAndDer(Ri(t),Rj(t),Rij,vi(t),vj(t)))

function [c,dc]=costAndDer(Ri,Rj,Rij,vi,vj)
[c,gradc]=rotLocFrobCostPair(Ri,Rj,Rij);
dc=trace(vi'*gradc(:,:,1))+trace(vj'*gradc(:,:,2));

