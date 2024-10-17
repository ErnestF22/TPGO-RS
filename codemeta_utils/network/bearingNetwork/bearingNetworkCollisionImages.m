%Compute bearing vectors and angles for collision images in a network
%function [y,a]=bearingNetworkCollisionImages(x,E,r)
function [y,a]=bearingNetworkCollisionImages(x,E,r)
[y,ny]=cnormalize(x(:,E(:,2),:)-x(:,E(:,1),:));
a=asin(min(1,r./ny));
