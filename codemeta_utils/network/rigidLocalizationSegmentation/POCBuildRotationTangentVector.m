function n=POCBuildRotationTangentVector(y,R,v,nodeList)
[dimData,NNodes]=size(y);
dimCoordinates=locSegDimData2DimCoordinates(dimData);
n=zeros(dimCoordinates*NNodes,1);

if dimData==2
    error('This function has not been implemented yet')
end

NNodeList=length(nodeList);
idxNodeList=ones(dimCoordinates,1)*((nodeList-1)*dimCoordinates)...
    +(1:dimCoordinates)'*ones(1,NNodeList);

for iNodeList=1:NNodeList
    iNode=nodeList(iNodeList);
    yi=y(:,iNode);
    Ri=R(:,:,iNode);
    n(idxNodeList(:,iNodeList))=[-hat(yi)*v; Ri'*v];
end
