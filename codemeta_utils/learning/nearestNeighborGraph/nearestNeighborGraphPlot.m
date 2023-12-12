function nearestNeighborGraphPlot(x,idxList,kList,varargin)
NPoints=size(x,2);
if ~exist('kList','var') || isempty(kList)
    kList=size(idxList,2)*ones(NPoints,1);
end

plotPoints(x,varargin{:})
NEdges=sum(kList);
E=zeros(NEdges,2);
idxKList=[0;cumsum(kList)];
for iPoint=1:NPoints
    k=kList(iPoint);
    E(idxKList(iPoint)+1:idxKList(iPoint+1),:)=...
        [iPoint*ones(k,1) idxList(iPoint,1:k)'];
end

plotEdges(x,E,varargin{:})
