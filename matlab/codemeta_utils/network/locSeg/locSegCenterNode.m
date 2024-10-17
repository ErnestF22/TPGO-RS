function MCentered=locSegCenterNode(M,dimCoordinates,nodeCenter)
bNN=null(M(dimCoordinates*(nodeCenter-1)+(1:dimCoordinates),:));
MCentered=M*bNN;
