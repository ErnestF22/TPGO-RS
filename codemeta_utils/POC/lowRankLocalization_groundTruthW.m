%Create low-rank shape-of-motion matrix
%function W=lowRankLocalization_groundTruthW(Ri,Ti,WInfo)
function [W]=lowRankLocalization_groundTruthW(Ri,Ti,WInfo)
nbNodes=size(Ri,3);
dim=WInfo.dimAmbient;
if nbNodes ~= WInfo.nbNodes
    error('Size of Ri not compatible with WInfo')
end
W=zeros(dim*nbNodes,nbNodes);
idxRows=reshape(1:dim*nbNodes,dim,nbNodes);

for iNode=1:nbNodes
    W(idxRows(:,iNode),:)=Ri(:,:,iNode)'*Ti;
end
if WInfo.flagRotationAugmented
    W=[W matStack(multitransp(Ri))];
end
