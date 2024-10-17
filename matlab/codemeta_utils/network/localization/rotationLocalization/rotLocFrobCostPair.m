function [c,gradc]=rotLocFrobCostPair(Ri,Rj,Rij)
flagComputeGrad=nargout>1;

E=Ri*Rij-Rj;
c=sum(E(:).^2)/2;

if flagComputeGrad
    gradc=rotLocFrobCostPair_grad(Ri,Rj,Rij,E);
end
