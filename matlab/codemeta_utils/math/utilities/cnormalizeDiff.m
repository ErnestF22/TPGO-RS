%Compute the differential of the normalization operation
%function [dxNormalized,DxNormalized]=cnormalizeDiff(x,dx)
%Compute the d/dt cnormalize(x) from x and d/dt x
function [dxNormalized,DxNormalized]=cnormalizeDiff(x,dx)
DxNormalized=cnormalizeDiffMat(x);
dxNormalized=DxNormalized*dx;

