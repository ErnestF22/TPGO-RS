function n=POCbuildScalingTangentVector(y,nodeList)
[dimData,NNodes]=size(y);
NNodeList=length(nodeList);
idxTranslation=locSegTangentVectorIdx(dimData,NNodes);
dimCoordinates=locSegDimData2DimCoordinates(dimData);
n=zeros(dimCoordinates*NNodes,1);
for iNodeList=1:NNodeList
    iNode=nodeList(iNodeList);
    n(idxTranslation(:,iNode))=y(:,iNode);
end
