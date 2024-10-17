function C=POCcomputeIntersectionEdgesScaling(E,y,bN)
[dimData,NNodes]=size(y);
NEdges=size(E,1);
idxTranslation=locSegTangentVectorIdx(dimData,NNodes);
C=cell(1,NEdges);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    yi=y(:,iNode);
    yj=y(:,jNode);
    Ni=bN(idxTranslation(:,iNode),:);
    Nj=bN(idxTranslation(:,jNode),:);
    C{iEdge}=(yi'*Ni/(yi'*yi)-yj'*Nj/(yj'*yj));
end
