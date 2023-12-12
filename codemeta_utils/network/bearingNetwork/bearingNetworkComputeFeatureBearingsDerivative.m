%Compute derivative of bearings for a network
%function [dy,R,B]=bearingNetworkComputeFeatureBearingsDerivative(dx,E,y,ny)
%Inputs
%   dx   derivative of agents' locations (i.e., \dot{x})
%   y   current feature bearing
%   ny  ranges corresponding to the bearings
%   E   network edge list
%Output
%   dy  derivative of y along v
%   D,R matrices to build the differential Dy=inv(D)*R that maps dx to dy
function [dy,D,R]=bearingNetworkComputeFeatureBearingsDerivative(dx,E,y,ny)
[d,NNodes]=size(dx);
[D,R]=bearingNetworkComputeFeatureBearingsDerivativeMat(E,y,ny,NNodes);
dy=reshape(D\R*dx(:),d,[]);
