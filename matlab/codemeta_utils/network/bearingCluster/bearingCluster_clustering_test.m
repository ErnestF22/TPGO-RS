function bearingCluster_clustering_test
sigmaNoise=0.0;

d=2;
datasetName='complex-loop-5-2';
%datasetName='loop-4';
%datasetName='trapezoid';
%datasetName='degenerate-square';
%datasetName='butterfly';
%datasetName='3-connected';
%datasetName='3-connected-skew';

[A,x]=bearingCluster_generateTest(datasetName,'d',d);
E=adj2edges(A,'oriented');

%generate measurements
u=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);

%cluster
[membership,L]=bearingCluster_clustering(E,u);

fprintf('Framework contains:\n')
fprintf('\t%d shakes\n',size(L,2));
fprintf('\t%d rigid components\n',length(unique(membership)))
bearingClustering_plot(x,E,membership)
