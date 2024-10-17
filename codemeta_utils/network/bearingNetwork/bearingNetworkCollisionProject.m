function uProj=bearingNetworkCollisionProject(u,y,a,E)
NNodes=size(u,2);
uProj=zeros(size(u));
for iNode=1:NNodes
    [yi,ai]=bearingNetworkCollisionCollectImagesNode(y,a,E,iNode);
    uProj(:,iNode)=bearingSectorsProject(u(:,iNode),yi,ai);
end
