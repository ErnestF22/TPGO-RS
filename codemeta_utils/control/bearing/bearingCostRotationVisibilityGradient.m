%Evaluate gradient for the rotation visibility constraint cost
%function gradPhiVec=bearingCostRotationVisibilityGradient(R,y,y0,funs,c)
%Inputs
%   R       body to world rotation
%   y       normalized bearing vectors
%   ny      range information corresponding to y
%   y0      normalized vector in body coordinates pointing at center of
%           constraint
%   funs    constraint function structure
%   c       inner products between y and y0 (can be omitted)
function gradPhiVec=bearingCostRotationVisibilityGradient(R,y,y0,funs,c)
if ~exist('c','var')
    c=bearingComputeCosine(y,R*y0);
end
a=bearingCostRotationVisibilityGradientRTerms(R,y,y0,funs,c);
gradRPhiVec=sum(ny.*a);
gradxPhi=bearingCostGeneral_gradient(y,R*y0,funs);
gradPhiVec=[gradRPhiVec;gradxPhi];
