%Compute range information from base point and landmarks
%function [ny]=bearingComputeRanges(x,xLandmarks)
%Inputs
%   x           coordinates of base point
%   xLandmarks  coordinates of landmarks
%Outputs
%   ny      ranges corresponding to the bearings
%
function [ny]=bearingComputeRanges(x,xLandmarks)
Nx=size(xLandmarks,2);
uVec=ones(1,Nx);
ny=sqrt(sum((xLandmarks-x*uVec).^2));
