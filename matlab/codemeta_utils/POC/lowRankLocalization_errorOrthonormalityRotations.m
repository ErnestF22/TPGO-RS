function e=lowRankLocalization_errorOrthonormalityRotations(WInfo,W)
idxMatrix=WInfo.idxMatrix;
WRot=W(permute(idxMatrix(:,:,WInfo.jNodeR),[1 3 2]));
e=shiftdim(sqrt(sum(sum((multiprod(WRot,multitransp(WRot))-eye(3)).^2),2)));
