function [phi,gradPhi]=bearingNetworkCostRangesPair(y,yg,ny,nyg,funs)
flagComputeGrad=false;
if nargout>1
    flagComputeGrad=true;
end
phi=bearingCostGeneralRanges(y,yg,ny,nyg,funs);
if flagComputeGrad
    gradPhi=bearingNetworkCostRangesPairGradient(y,yg,ny,nyg,funs);
end
