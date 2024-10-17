function n=locSegBuildRotationTangentVector(y,R,v,members)
[dimData,NNodes]=size(y);
dimCoordinates=locSegDimData2DimCoordinates(dimData);
n=zeros(dimCoordinates*NNodes,1);
[~,~,idxCoordinates]=locSegTangentVectorIdx(dimData,NNodes);
S=[0 -1; 1 0];

for iNode=1:NNodes
    yi=y(:,iNode);
    Ri=R(:,:,iNode);
    iAxis=find(members(:,iNode),1,'first');
    vi=v(:,iAxis);
    switch dimData
        case 2
            n(idxCoordinates(:,iNode))=[S*yi*vi; vi];            
        case 3
            n(idxCoordinates(:,iNode))=[-hat(yi)*vi; Ri'*vi];
    end
end
