%Compute normalized bearing vectors for features in a network
%function [y,ny]=bearingNetworkComputeFeatureBearings(x,xFeatures,E)
%Inputs
%   x           coordinates of the nodes
%   xFeatures   corrdinates of the features
%   E           edge list of the computed bearings
%       First column    index of the agent in x
%       Second column   index of the feature in xFeatures
function [y,ny]=bearingNetworkComputeFeatureBearings(x,xFeatures,E)
[y,ny]=cnormalize(xFeatures(:,E(:,2))-x(:,E(:,1)));
