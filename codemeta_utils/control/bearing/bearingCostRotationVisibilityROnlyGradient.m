%Evaluate gradient for the rotation visibility constraint cost (no ranges)
%function gradPhiRVec=bearingCostRotationVisibilityROnlyGradient(R,y,y0,funs,c)
%The gradient is computed with respect to R only.
%Inputs
%   R       body to world rotation
%   y       normalized bearing vectors
%   y0      normalized vector in body coordinates pointing at center of
%           constraint
%   funs    constraint function structure
%   c       inner products between y and y0 (can be omitted)
function gradPhiRVec=bearingCostRotationVisibilityROnlyGradient(R,y,y0,funs,c)
if ~exist('c','var')
    c=bearingComputeCosine(y,R*y0);
end
NLandmarks=size(y,2);
S=[0 -1; 1 0];
gradPhiRVec=0;
for iLandmark=1:NLandmarks
    gradPhiRVec=gradPhiRVec+funs.df(c(iLandmark))*y(:,iLandmark)'*R*S*y0;
end
