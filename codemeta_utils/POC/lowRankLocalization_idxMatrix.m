%Make matrix to index blocks in the shape-of-motion matrix
function idxMatrix=lowRankLocalization_idxMatrix(WInfo)
nbNodes=WInfo.nbNodes;
dim=WInfo.dimAmbient;
idxMatrix=reshape(1:WInfo.numel,dim,nbNodes,WInfo.size(2));
