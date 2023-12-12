%Evaluate cost and gradient for the rotation visibility constraint cost (no ranges)
%function [phi,gradPhiR]=bearingCostRotationVisibilityROnly(R,y,y0,funs)
%The gradient is computed with respect to R only.
%Inputs
%   R       body to world rotation (i.e., reference interpretation)
%   y       normalized bearing vectors
%   y0      normalized vector in body coordinates pointing at center of
%           constraint
%   funs    constraint function structure
function [phi,gradPhiR]=bearingCostRotationVisibilityROnly(R,y,y0,funs)
flagComputeGrad=false;
if nargout>1
    flagComputeGrad=true;
end

c=bearingComputeCosine(y,R*y0);
phi=sum(funs.f(c));
if flagComputeGrad
    gradPhiR=bearingCostRotationVisibilityROnlyGradient(R,y,y0,funs,c);
end
