%Compute normalized bearing vectors in a network
%function [y,ny]=bearingNetworkComputeBearings(x,E)
function [y,ny]=bearingNetworkComputeBearings(x,E)
[y,ny]=cnormalize(x(:,E(:,2),:)-x(:,E(:,1),:));
