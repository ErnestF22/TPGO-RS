%Evaluate cost and gradient for the rotation visibility constraint cost
%function [phi,gradPhiVec]=bearingCostRotationVisibility(R,y,ny,y0,funs)
%Inputs
%   R       body to world rotation (i.e., reference interpretation)
%   y       normalized bearing vectors
%   ny      range information corresponding to y
%   y0      normalized vector in body coordinates pointing at center of
%           constraint
%   funs    constraint function structure
function [phi,gradPhiVec]=bearingCostRotationVisibility(R,y,ny,y0,funs)
flagComputeGrad=false;
if nargout>1
    flagComputeGrad=true;
end

c=bearingComputeCosine(y,R*y0);
phi=ny*funs.f(c)';
if flagComputeGrad
    gradPhiVec=bearingCostRotationVisibilityGradient(R,y,ny,y0,funs,c);
end