%Build matrix from translation measurements such that AMeasurements*vec(W)=bMeasurements
%function [AMeasurements,bMeasurements]=lowRankLocalization_translations_constraints(WInfo,tij,E)
function [AMeasurements,bMeasurements]=lowRankLocalization_translations_constraints(WInfo,Tij,E)
nbEdges=size(E,1);
dim=WInfo.dimAmbient;
nbConstraints=dim*nbEdges;
idxA=reshape(1:nbConstraints,dim,nbEdges);
AMeasurements=zeros(nbConstraints,WInfo.numel);
bMeasurements=zeros(nbConstraints,1);
for iEdge=1:nbEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    TijEdge=Tij(:,iEdge);
    [Aij,bij]=lowRankLocalization_translations_singleConstraint(WInfo,TijEdge,iNode,jNode);
    AMeasurements(idxA(:,iEdge),:)=Aij;
    bMeasurements(idxA(:,iEdge))=bij;
end
