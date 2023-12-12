%Compute the matrix form of the differential of the normalization operation
%function [DxNormalized,xNormalized]=cnormalizeDiffMat(x)
%Compute the matrix DxNormalized such that
%d/dt cnormalize(x)=DxNormalized*d\dt x
function [DxNormalized,xNormalized]=cnormalizeDiffMat(x)
[xNormalized,nx]=cnormalize(x);
d=size(x,1);
if nx==0
    DxNormalized=zeros(d);
else
    DxNormalized=(eye(d)-xNormalized*xNormalized')/nx;
end
