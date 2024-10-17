function n=locSegBuildScalingTangentVector(y,C,members)
[dimData,NNodes]=size(y);
dimCoordinates=locSegDimData2DimCoordinates(dimData);
n=zeros(dimCoordinates*NNodes,1);
idxTranslations=locSegTangentVectorIdx(dimData,NNodes);

for iNode=1:NNodes
    yi=y(:,iNode);
    ci=C(find(members(:,iNode),1,'first'));
    n(idxTranslations(:,iNode))=ci*yi;
end
