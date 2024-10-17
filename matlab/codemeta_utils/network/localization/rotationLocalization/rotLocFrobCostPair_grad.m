function gradc=rotLocFrobCostPair_grad(Ri,Rj,Rij,E)
if ~exist('E','var')
    E=Ri*Rij-Rj;
end

gradc=cat(3,E*Rij',-E);

