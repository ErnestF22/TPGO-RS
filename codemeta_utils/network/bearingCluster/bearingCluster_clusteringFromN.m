function [membership,X]=bearingCluster_clusteringFromN(N)
K=size(N,1);
XInit=clusteringInitWithFurthest(N,K);
[membership,X]=kmeans(N',K,'start',XInit','emptyAction','singleton');
X=X';
membership=membership';


