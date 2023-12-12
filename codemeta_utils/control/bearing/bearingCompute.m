%Compute bearing information from base point and landmarks
%function [y,ny]=bearingCompute(x,xLandmarks)
%Inputs
%   x           coordinates of base point
%   xLandmarks  coordinates of landmarks
%Outputs
%   y       bearings
%   ny      ranges corresponding to the bearings
%
function [y,ny]=bearingCompute(x,xLandmarks)
Nx=size(x,2);
NxLandmarks=size(xLandmarks,2);
if Nx==1 && NxLandmarks>1
    z=x*ones(1,NxLandmarks);
elseif NxLandmarks==1 && Nx>1
    xLandmarks=xLandmarks*ones(1,Nx);
end
[y,ny]=cnormalize(xLandmarks-z);
