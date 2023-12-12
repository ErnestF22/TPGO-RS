function [membership,info]=bearingCluster(x,E,sigmaNoise,tol)
%Cluster a framework (or subframework) into interdependent edges @ARMAN
[u,lambda]=bearingCluster_getBearingsScalesFromE(x,E,'noisy',sigmaNoise);

CUnsigned=grCycleBasisMultiComp(E);
C=grOrientCycleBasis(CUnsigned,E)';
M=bearingCluster_measurementMatrixFromC(C,u);

[L,s]=bearingCluster_nullSpaceBasis(M,tol);

N=bearingCluster_clusterVectors(L);

[membership,info]=projective_quickshift(N,'threshold',0.5);
end