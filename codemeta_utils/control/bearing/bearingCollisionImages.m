%Compute the bearings and angles from the image of a landmark
%function [y,a]=bearingCollisionImageAngle(x,xLandmarks,rLandmarks)
function [y,a]=bearingCollisionImages(x,xLandmarks,rLandmarks)
Nx=size(x,2);
NxLandmarks=size(xLandmarks,2);
if Nx==1 && NxLandmarks>1
    x=x*ones(1,NxLandmarks);
elseif NxLandmarks==1 && Nx>1
    xLandmarks=xLandmarks*ones(1,Nx);
    rLandmarks=rLandmarks*ones(1,Nx);
end
[y,ny]=cnormalize(xLandmarks-x);
a=asin(min(1,rLandmarks./ny));
