%Build matrix from bearing measurements such that AMeasurements*vec(W)=bMeasurements
%function [AMeasurements,bMeasurements]=lowRankLocalization_bearings_constraints(WInfo,tij,E)
function [AMeasurements,bMeasurements]=lowRankLocalization_bearings_constraints(WInfo,tij,E)
nbEdges=size(E,1);
dim=WInfo.dimAmbient;
nbConstraints=(dim-1)*nbEdges;
idxA=reshape(1:nbConstraints,dim-1,nbEdges);
AMeasurements=zeros(nbConstraints,WInfo.numel);
for iEdge=1:nbEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    tijEdge=tij(:,iEdge);
    PPij=lowRankLocalization_bearings_singleConstraint(WInfo,tijEdge,iNode,jNode);
    AMeasurements(idxA(:,iEdge),:)=PPij;
end
bMeasurements=zeros(nbConstraints,1);
