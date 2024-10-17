function JCheck=POCBuildTranslationEdgesCheck(E)
NNodes=max(E(:));
NEdges=size(E,1);

dimData=3;
dimCoordinates=6;
idxBase=(0:NNodes-1)*dimCoordinates;
idxTranslation=(1:dimData)'*ones(1,NNodes)+ones(dimData,1)*idxBase;
idxEdges=reshape(1:NEdges*3,3,[]);

JCheck=zeros(NEdges*3,dimCoordinates*NNodes);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    idxRows=idxEdges(:,iEdge);
    idxCols=[idxTranslation(:,iNode); idxTranslation(:,jNode)];
    
    JCheck(idxRows,idxCols)=[eye(3) -eye(3)];
end
