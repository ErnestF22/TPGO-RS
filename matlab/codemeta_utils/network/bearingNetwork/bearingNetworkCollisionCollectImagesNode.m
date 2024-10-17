function [yNode,aNode]=bearingNetworkCollisionCollectImagesNode(y,a,E,iNode)
yNode=y(:,E(:,1)==iNode);
aNode=a(:,E(:,1)==iNode);
