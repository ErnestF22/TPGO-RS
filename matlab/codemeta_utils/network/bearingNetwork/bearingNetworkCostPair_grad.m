function [gradPhi,DgradPhi]=bearingNetworkCostPair_grad(yij,ygij,funs,nygij)

flagComputeDgrad=nargout>1;

if ~flagComputeDgrad
    gradPhiij=bearingCostGeneral_gradient(yij,ygij,funs);
else
    [gradPhiij,Hij]=bearingCostGeneral_gradient(yij,ygij,funs,nygij);
    DgradPhi=[Hij -Hij; -Hij' Hij];
end
gradPhi=[gradPhiij -gradPhiij];
