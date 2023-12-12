%Compute the total bearing cost, its gradient and its Hessian
%function [phi,gradPhi,DgradPhi,DygradPhi]=bearingCostGeneral(YEval,YGoal,nYEval,funs)
%Inputs
%   YEval   current measured bearing vectors
%   YGoal   bearing vectors at the goal position
%   nYEval  norm of measured bearing vectors (range information)
%   funs    reshaping function used for the cost
%Outputs
%   phi         value of the cost
%   gradPhi     gradient of the cost with respect to YEval
%   DgradPhi    differential of the gradient with respect to YEval
%   DygradPhi   differential of the gradient with respect to YGoal
%See also bearingCostGeneral_gradient, bearingCostGeneral_DgradientMeasurements
function [phi,gradPhi,DgradPhi,DygradPhi]=bearingCostGeneral(y,yg,ny,funs)
flagComputeGrad=false;
flagComputeDGrad=false;
flagComputeDyGrad=false;


if nargout>1
    flagComputeGrad=true;
end
if nargout>2
    flagComputeDGrad=true;
end
if nargout>2
    flagComputeDyGrad=true;
end

f=funs.f;

c=bearingComputeCosine(y,yg);
phi=sum(ny.*f(c));

if flagComputeGrad
    if ~flagComputeDGrad
        gradPhi=bearingCostGeneral_gradient(y,yg,funs);
    else
        [gradPhi,DgradPhi]=bearingCostGeneral_gradient(y,yg,funs,ny);
    end
    if flagComputeDyGrad
        DygradPhi=bearingCostGeneral_DgradientMeasurements(y,yg,funs);
    end
end
