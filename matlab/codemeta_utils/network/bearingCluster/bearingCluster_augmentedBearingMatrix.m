function U=bearingCluster_augmentedBearingMatrix(u)
[d,NEdges]=size(u);
U=zeros(d*NEdges,NEdges);
idxEdge=reshape(1:d*NEdges,d,NEdges);
for iEdge=1:NEdges
    U(idxEdge(:,iEdge),iEdge)=u(:,iEdge);
end
