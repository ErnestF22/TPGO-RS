%Compute bearing+range cost and its gradient
%function [phi,gradPhi]=bearingCostGeneralRanges(y,yg,ny,nyg,funs)
function [phi,gradPhi]=bearingCostGeneralRanges(y,yg,ny,nyg,funs)
flagComputeGrad=false;
if nargout>1
    flagComputeGrad=true;
end

f=funs.f;

c=bearingComputeCosine(y,yg);
e=ny.*c-nyg;

phi=sum(f(e));

if flagComputeGrad
    gradPhi=bearingCostGeneralRangesGradient(y,yg,ny,nyg,funs);
end