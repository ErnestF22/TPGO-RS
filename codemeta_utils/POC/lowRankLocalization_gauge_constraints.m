%Return constraints for fixing global gauge ambiguity
function [AFix,bFix]=lowRankLocalization_gauge_constraints(WInfo)
iNodeFix=WInfo.iNodeFix;
dim=WInfo.dimAmbient;
idxMatrix=WInfo.idxMatrix;

%Since the first translation is at the origin, the entire first column of
%the matrix W must be zero
nbConstraintsTranslation=WInfo.size(1);
AFix=zeros(nbConstraintsTranslation,WInfo.numel);
AFix(sub2ind(size(AFix),(1:nbConstraintsTranslation)',vec(idxMatrix(:,:,iNodeFix))))=1;
bFix=zeros(nbConstraintsTranslation,1);

%we can fix the rotation only if the rotation columns are present
if WInfo.flagRotationAugmented
    nbConstraintsRotation=dim^2;
    AFixRotation=zeros(nbConstraintsRotation,WInfo.numel);
    AFixRotation(sub2ind(size(AFixRotation),(1:nbConstraintsRotation)',vec(idxMatrix(:,iNodeFix,WInfo.jNodeR))))=1;
    bFixRotation=vec(eye(3));
    %append constraints
    AFix=[AFix;AFixRotation];
    bFix=[bFix;bFixRotation];
end
