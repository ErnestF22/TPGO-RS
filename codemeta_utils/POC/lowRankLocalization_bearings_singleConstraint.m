%Makes the constraint for a single bearing measurement
function PPij=lowRankLocalization_singleConstraintDirection(WInfo,tij,iNode,jNode)
idxMatrix=lowRankLocalization_idxMatrix(WInfo);
Pij=orthComplement(tij)';

PPij=zeros(size(Pij,1),WInfo.numel);
PPij(:,idxMatrix(:,iNode,iNode))=Pij;
PPij(:,idxMatrix(:,iNode,jNode))=-Pij;
