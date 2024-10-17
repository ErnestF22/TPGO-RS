function [phi,gradPhiVec]=bearingNetworkCostRotationVisibilityPair(Ri,Rj,y,ny,y0,funs)
flagComputeGrad=false;
if nargout>1
    flagComputeGrad=true;
end
phi=bearingCostRotationVisibility(Ri,y,ny,y0,funs)...
    +bearingCostRotationVisibility(Rj,-y,ny,y0,funs);
if flagComputeGrad
    gradPhiVec=bearingNetworkCostRotationVisibilityPairGradient(Ri,Rj,y,ny,y0,funs);
end
