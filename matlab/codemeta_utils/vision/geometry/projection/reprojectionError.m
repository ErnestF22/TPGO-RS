%Compute reprojection error and, optionally, the Jacobians and the RMSE
%function [err,Jx]=reprojectionError(xTruth,P,X)
%Inputs
%   xTruth  [d x NX x NP] array
function [err,JxProjected,rmse]=reprojectionError(xTruth,P,X)
flagComputeJacobian=false;
flagComputeRMSE=false;
if nargout>1
    flagComputeJacobian=true;
end
if nargout>2
    flagComputeRMSE=true;
end

if flagComputeJacobian
    [xProjected,JxProjected]=projectFromP(P,X);
else
    xProjected=projectFromP(P,X);
end

err=xProjected-xTruth;

if flagComputeRMSE
    rmse=squeeze(sqrt(mean(sum(err.^2))));
end
