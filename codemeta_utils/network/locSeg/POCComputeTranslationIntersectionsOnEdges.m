function C=POCComputeTranslationIntersectionsOnEdges(E,bNTranslations)
NNodes=max(E(:));
NEdges=size(E,1);

dimData=3;
dimCoordinates=6;
idxBase=(0:NNodes-1)*dimCoordinates;
idxTranslation=(1:dimData)'*ones(1,NNodes)+ones(dimData,1)*idxBase;


C=cell(1,NEdges);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    idxINode=idxTranslation(:,iNode);
    idxJNode=idxTranslation(:,jNode);
    
    C{iEdge}=null(bNTranslations(idxINode,:)-bNTranslations(idxJNode,:));
end