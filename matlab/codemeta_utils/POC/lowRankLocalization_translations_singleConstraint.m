%Makes the constraint for a single relative translation measurement
function [Aij,bij]=lowRankLocalization_singleConstraintTranslation(WInfo,Tij,iNode,jNode)
idxMatrix=lowRankLocalization_idxMatrix(WInfo);
dim=WInfo.dimAmbient;
I=eye(dim);
Aij=zeros(dim,WInfo.numel);
Aij(:,idxMatrix(:,iNode,iNode))=-I;
Aij(:,idxMatrix(:,iNode,jNode))=I;
bij=Tij;
