function BReduced=bearingCluster_reducedMatrix(B,idxAnchorNode)
if ~exist('idxAnchorNode','var')
    idxAnchorNode=1;
end

NNodes=size(B,2);
idxNonAnchorNodes=setdiff(1:NNodes,idxAnchorNode);

BReduced=B(:,idxNonAnchorNodes);

