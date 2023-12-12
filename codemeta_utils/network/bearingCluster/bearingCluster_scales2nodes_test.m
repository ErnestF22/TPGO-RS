function bearingCluster_scales2nodes_test

%datasetName='butterfly-rigid';
datasetName='3-connected-skew';
sigmaNoise=0.1;

[A,x]=bearingCluster_generateTest(datasetName);
E=adj2edges(A,'oriented');
u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);
M=bearingCluster_measurementMatrixFromE(E,u);

[U,S,V]=svd(M);
l=V(:,end);
if median(sign(l))<0
    l=-l;
end

x=bearingCluster_scales2nodes(E,u,l);

bearingClustering_plot(x,E)

