function [phi,gradPhi]=bearingNetworkCostPair(yij,ygij,nyij,funs)
flagComputeGrad=false;
if nargout>1
    flagComputeGrad=true;
end
phi=bearingCostGeneral(yij,ygij,nyij,funs);
if flagComputeGrad
    gradPhi=bearingNetworkCostPair_grad(yij,ygij,funs);
end
