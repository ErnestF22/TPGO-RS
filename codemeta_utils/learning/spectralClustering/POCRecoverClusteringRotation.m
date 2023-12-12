function POCRecoverClusteringRotation
NPointsCluster=[5 3 3];
NPoints=sum(NPointsCluster);
NClusters=length(NPointsCluster);

idxClusters=[0 cumsum(NPointsCluster)];
VTruth=zeros(NPoints,NClusters);
for iCluster=1:NClusters
    VTruth([idxClusters(iCluster)+1:idxClusters(iCluster+1)],iCluster)=1;
end

RTruth=rot_randn(NClusters);
V=VTruth*RTruth;

R=zeros(NClusters);
idxR=floor(rand*NPoints)+1;
for iCluster=1:NClusters
    
