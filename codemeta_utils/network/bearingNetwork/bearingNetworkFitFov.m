function R=bearingNetworkFitFov(x,E,y0)
y=bearingNetworkComputeBearings(x,E);
NNodes=size(x,2);

R=zeros(2,2,NNodes);
for iNode=1:NNodes
    yNode=y(:,E(:,1)==iNode);
    yHeading=bearingFitFov(yNode);
    R(:,:,iNode)=householderRotation(y0,yHeading);
end
