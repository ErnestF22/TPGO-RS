function [idxList,kList]=nearestNeighborGraphFromAdjacency(A)
NPoints=size(A,2);
kList=full(sum(A)');
idxList=zeros(NPoints,max(kList));
for iPoint=1:NPoints
    idxList(iPoint,1:kList(iPoint))=find(A(:,iPoint));
end
