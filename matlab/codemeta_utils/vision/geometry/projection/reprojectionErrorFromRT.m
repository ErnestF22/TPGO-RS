%Compute reprojection error, Jacobians and RMSE from rotation and translation
function [err,JRTProjected,rmse]=reprojectionErrorFromRT(xTruth,R,T,X,varargin)

flagComputeJacobian=false;
flagComputeRMSE=false;
if nargout>1
    flagComputeJacobian=true;
end
if nargout>2
    flagComputeRMSE=true;
end

if flagComputeJacobian
    [xProjected,JRTProjected]=projectFromRT(R,T,X);
else
    xProjected=projectFromRT(R,T,X);
end

err=xProjected-xTruth;

if flagComputeRMSE
    rmse=squeeze(sqrt(mean(sum(err.^2))));
end
