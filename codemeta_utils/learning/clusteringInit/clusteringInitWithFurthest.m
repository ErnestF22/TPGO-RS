function [X,idxCluster]=clusteringInitWithFurthest(A,K)

%select first point as the point which is further from the median
m=median(A,2);
[~,idxCluster]=max(euclideanDistMatrix(m,A));

%select other points as the furthest ones from those already selected
distSq=euclideanDistMatrix(A,A);
NA=size(A,2);
idxCluster(2:K)=NaN;

for k=2:K
    idxNotCluster=setdiff(1:NA,idxCluster(1:k-1));
    [~,idxIdxClusterNext]=max(min(distSq(idxCluster(1:k-1),idxNotCluster),[],1));
    idxCluster(k)=idxNotCluster(idxIdxClusterNext);
end

X=A(:,idxCluster);
