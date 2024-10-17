function x=bearingCluster_embed(E,u)
C=grOrientedCycleBasis(E)';
M=bearingCluster_measurementMatrixFromC(C,u);
L=null(M);

l=L*cnormalize(randn(size(L,2),1));

x=bearingCluster_scales2nodes(E,u,l);
